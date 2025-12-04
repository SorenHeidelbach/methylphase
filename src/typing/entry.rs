use crate::typing::cli::{Commands, RunArgs};
use crate::typing::data::{self, load_config, load_dataset, Dataset};
use crate::typing::em::{self, fit_em, EmSettings};
use crate::typing::fastq_split::split_fastq;
use crate::typing::floria::parse_floria;
use crate::typing::impute::impute_dataset;
use crate::typing::methylation::{add_methylation_features, load_methylation_features};
use crate::typing::model_selection::{Criterion, ModelSelectionStrategy, SelectionStrategy};
use crate::typing::priors::DirichletPriors;
use anyhow::{anyhow, Result};
use crate::commands::split_reads;
use crate::cli::ClusterAlgorithm;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

/// Dispatch typing commands previously housed under the methyltyping binary.
pub fn run(command: Commands) -> Result<()> {
    match command {
        Commands::Fit(args) => {
            log_step(format!(
                "Loading config {} and dataset {}",
                args.config.display(),
                args.data.display()
            ));
            let config = load_config(&args.config)?;
            let dataset = load_dataset(&args.data, &config, args.delimiter_byte())?;
            log_step(describe_dataset(&dataset));
            let priors = DirichletPriors {
                alpha_pi: args.alpha_pi,
                alpha_phi: args.alpha_phi,
            };
            let settings = EmSettings {
                max_iters: args.max_iter,
                tol: args.tol,
                min_class_weight: args.min_class_weight,
            };
            log_step(format!(
                "Fitting EM with C = {} (max_iters = {}, tol = {})",
                args.classes, args.max_iter, args.tol
            ));
            let fit_start = Instant::now();
            let result = fit_em(&dataset, args.classes, &priors, settings)?;
            log_step(format!(
                "Fit complete in {:.2?}: log-likelihood = {:.4}, iterations = {}",
                fit_start.elapsed(),
                result.log_likelihood,
                result.iterations
            ));
            let model = em::ModelFile::from_parts(&config, &result.params);
            log_step(format!("Writing model to {}", args.output_model.display()));
            write_model(&args.output_model, &model)?;
            if let Some(path) = args.output_responsibilities {
                log_step(format!(
                    "Writing responsibilities to {}",
                    path.display()
                ));
                em::write_responsibilities(&path, &dataset, &result)?;
            }
        }
        Commands::SelectBestC(args) => {
            log_step(format!(
                "Loading config {} and dataset {}",
                args.config.display(),
                args.data.display()
            ));
            let config = load_config(&args.config)?;
            let dataset = load_dataset(&args.data, &config, args.delimiter_byte())?;
            log_step(describe_dataset(&dataset));
            let priors = DirichletPriors {
                alpha_pi: args.alpha_pi,
                alpha_phi: args.alpha_phi,
            };
            let settings = EmSettings {
                max_iters: args.max_iter,
                tol: args.tol,
                min_class_weight: args.min_class_weight,
            };
            let strategy = SelectionStrategy {
                min_classes: args.min_classes,
                max_classes: args.max_classes,
                priors: priors.clone(),
                settings,
                criterion: match args.criterion.as_str() {
                    "icl" => Criterion::Icl,
                    "cv" => Criterion::CrossValidated,
                    _ => Criterion::Bic,
                },
                penalty_multiplier: args.penalty_multiplier,
                cv_folds: args.cv_folds,
            };
            let criterion_label = match strategy.criterion {
                Criterion::Bic => "BIC",
                Criterion::Icl => "ICL",
                Criterion::CrossValidated => "CV-NLL",
            };
            log_step(format!(
                "Selecting best C in [{}..={}] using {}",
                args.min_classes, args.max_classes, criterion_label
            ));
            let (best_c, result, best_score) = match strategy.criterion {
                Criterion::CrossValidated => {
                    let mut best: Option<(usize, f64)> = None;
                    for c in args.min_classes..=args.max_classes {
                        let start = Instant::now();
                        let score =
                            strategy.cross_validated_score(&dataset, c, &|ds, cls| {
                                fit_em(ds, cls, &priors, settings)
                            })?;
                        log_step(format!(
                            "C = {} cross-val in {:.2?}: {} = {:.4}",
                            c,
                            start.elapsed(),
                            criterion_label,
                            score
                        ));
                        match best {
                            None => best = Some((c, score)),
                            Some((_, best_score)) if score < best_score => {
                                best = Some((c, score))
                            }
                            _ => {}
                        }
                    }
                    let (best_c, best_score) = best
                        .ok_or_else(|| anyhow!("no cross-validated fits ran"))?;
                    let fit_start = Instant::now();
                    let res = fit_em(&dataset, best_c, &priors, settings)?;
                    log_step(format!(
                        "Refit best C = {} on full data in {:.2?}: log-likelihood = {:.4}",
                        best_c,
                        fit_start.elapsed(),
                        res.log_likelihood
                    ));
                    (best_c, res, best_score)
                }
                _ => strategy.select_best_model(&dataset, |ds, c| {
                    let start = Instant::now();
                    let res = fit_em(ds, c, &priors, settings)?;
                    let score = strategy.score(ds, &res, c);
                    log_step(format!(
                        "C = {} finished in {:.2?}: log-likelihood = {:.4}, {} = {:.4}, iterations = {}",
                        c,
                        start.elapsed(),
                        res.log_likelihood,
                        criterion_label,
                        score,
                        res.iterations
                    ));
                    Ok(res)
                })?,
            };
            let model = em::ModelFile::from_parts(&config, &result.params);
            log_step(format!(
                "Best C = {} ({} = {:.4}, log-likelihood = {:.4})",
                best_c, criterion_label, best_score, result.log_likelihood
            ));
            log_step(format!("Writing best model to {}", args.output_model.display()));
            write_model(&args.output_model, &model)?;
            println!(
                "Selected C = {} with log-likelihood {:.4}",
                best_c,
                result.log_likelihood
            );
        }
        Commands::Impute(args) => {
            log_step(format!(
                "Loading config {}, dataset {}, and model {}",
                args.config.display(),
                args.data.display(),
                args.model.display()
            ));
            let config = load_config(&args.config)?;
            let dataset = load_dataset(&args.data, &config, args.delimiter_byte())?;
            log_step(describe_dataset(&dataset));
            let model: em::ModelFile = read_model(&args.model)?;
            log_step("Imputing missing values using provided model");
            let imputed = impute_dataset(&dataset, &model.params)?;
            log_step(format!("Writing imputed dataset to {}", args.output.display()));
            data::write_dataset(&args.output, &imputed, args.delimiter_byte())?;
        }
        Commands::ConvertFloria(args) => {
            log_step(format!("Parsing Floria haploset {}", args.floria.display()));
            let parsed = parse_floria(&args.floria)?;
            log_step(describe_dataset(&parsed.dataset));
            let (dataset, config) = if let Some(meth_path) = &args.methylation {
                log_step(format!(
                    "Loading methylation features {}",
                    meth_path.display()
                ));
                let (vals, motifs) =
                    load_methylation_features(meth_path, args.delimiter_byte())?;
                log_step("Combining haploset with methylation features");
                add_methylation_features(&parsed.dataset, &parsed.config, &vals, &motifs)?
            } else {
                (parsed.dataset, parsed.config)
            };
            log_step(describe_dataset(&dataset));
            data::write_dataset(&args.output_data, &dataset, args.delimiter_byte())?;
            let toml_str = toml::to_string_pretty(&config)?;
            log_step(format!(
                "Wrote dataset to {} and config to {}",
                args.output_data.display(),
                args.output_config.display()
            ));
            std::fs::write(&args.output_config, toml_str)?;
        }
        Commands::Run(args) => {
            let methylation = prepare_methylation(&args)?;
            run_pipeline(args, methylation)?;
        }
        Commands::SplitFastq(args) => {
            log_step(format!(
                "Splitting FASTQ {} using assignments {}",
                args.fastq.display(),
                args.assignments.display()
            ));
            split_fastq(&args.fastq, &args.assignments, args.delimiter_byte(), &args.out)?;
        }
    }
    Ok(())
}

fn write_model(path: &std::path::Path, model: &em::ModelFile) -> Result<()> {
    let mut file = File::create(path)?;
    let serialized = serde_json::to_string_pretty(model)?;
    file.write_all(serialized.as_bytes())?;
    Ok(())
}

fn read_model(path: &Path) -> Result<em::ModelFile> {
    let file = File::open(path)?;
    let model: em::ModelFile = serde_json::from_reader(file)?;
    Ok(model)
}

fn run_pipeline(args: RunArgs, methylation_path: PathBuf) -> Result<()> {
    let pipeline_start = Instant::now();
    log_step(format!(
        "Starting pipeline; writing outputs under {}",
        args.out.display()
    ));
    std::fs::create_dir_all(&args.out)?;

    log_step(format!("Parsing Floria haploset {}", args.floria.display()));
    let parse_start = Instant::now();
    let parsed = parse_floria(&args.floria)?;
    log_step(format!(
        "Parsed Floria in {:.2?}. {}",
        parse_start.elapsed(),
        describe_dataset(&parsed.dataset)
    ));

    log_step(format!(
        "Loading methylation features {}",
        methylation_path.display()
    ));
    let meth_start = Instant::now();
    let (meth_vals, motifs) = load_methylation_features(&methylation_path, args.delimiter_byte())?;
    log_step(format!(
        "Loaded methylation features in {:.2?}",
        meth_start.elapsed()
    ));

    log_step("Combining haploset and methylation features");
    let merge_start = Instant::now();
    let (dataset, config) =
        add_methylation_features(&parsed.dataset, &parsed.config, &meth_vals, &motifs)?;
    log_step(format!(
        "Merged dataset in {:.2?}. {}",
        merge_start.elapsed(),
        describe_dataset(&dataset)
    ));

    // Subcategory names (hap only; methylation is continuous).
    let mut subcat_names = Vec::with_capacity(parsed.subcat_names.len());
    subcat_names.extend(parsed.subcat_names.clone());

    // Save dataset and config.
    let dataset_path = args.out.join("dataset.tsv");
    let config_path = args.out.join("categories.toml");
    log_step(format!(
        "Writing dataset to {} and config to {}",
        dataset_path.display(),
        config_path.display()
    ));
    data::write_dataset(&dataset_path, &dataset, args.delimiter_byte())?;
    let toml_str = toml::to_string_pretty(&config)?;
    std::fs::write(&config_path, toml_str)?;

    let priors = DirichletPriors {
        alpha_pi: args.alpha_pi,
        alpha_phi: args.alpha_phi,
    };
    let settings = EmSettings {
        max_iters: args.max_iter,
        tol: args.tol,
        min_class_weight: args.min_class_weight,
    };
    let strategy = SelectionStrategy {
        min_classes: args.min_classes,
        max_classes: args.max_classes,
        priors: priors.clone(),
        settings,
        criterion: match args.criterion.as_str() {
            "icl" => Criterion::Icl,
            "cv" => Criterion::CrossValidated,
            _ => Criterion::Bic,
        },
        penalty_multiplier: args.penalty_multiplier,
        cv_folds: args.cv_folds,
    };

    let fits_dir = args.out.join("fits");
    std::fs::create_dir_all(&fits_dir)?;
    let criterion_label = match strategy.criterion {
        Criterion::Bic => "BIC",
        Criterion::Icl => "ICL",
        Criterion::CrossValidated => "CV-NLL",
    };
    log_step(format!(
        "Fitting class counts {}..={} using {}",
        args.min_classes, args.max_classes, criterion_label
    ));

    let mut best: Option<(usize, em::EmResult, f64)> = None;
    let mut summary = Vec::new();
    for c in args.min_classes..=args.max_classes {
        let score_start = Instant::now();
        let score = match strategy.criterion {
            Criterion::CrossValidated => strategy.cross_validated_score(&dataset, c, &|ds, cls| {
                fit_em(ds, cls, &priors, settings)
            })?,
            _ => {
                // We still fit full data to emit artifacts; score uses same fit.
                // Fallback updated below for BIC/ICL.
                0.0
            }
        };

        let fit_start = Instant::now();
        let result = fit_em(&dataset, c, &priors, settings)?;
        let final_score = match strategy.criterion {
            Criterion::CrossValidated => score,
            _ => strategy.score(&dataset, &result, c),
        };
        let log_duration = if matches!(strategy.criterion, Criterion::CrossValidated) {
            format!(
                "C = {} CV score in {:.2?}: {} = {:.4}; fit full data in {:.2?}: log-likelihood = {:.4}, iterations = {}",
                c,
                score_start.elapsed(),
                criterion_label,
                final_score,
                fit_start.elapsed(),
                result.log_likelihood,
                result.iterations
            )
        } else {
            format!(
                "C = {} finished in {:.2?}: log-likelihood = {:.4}, {} = {:.4}, iterations = {}",
                c,
                fit_start.elapsed(),
                result.log_likelihood,
                criterion_label,
                final_score,
                result.iterations
            )
        };
        log_step(log_duration);

        let model = em::ModelFile::from_parts(&config, &result.params);
        let model_path = fits_dir.join(format!("model_C{}.json", c));
        log_step(format!("Writing model to {}", model_path.display()));
        write_model(&model_path, &model)?;
        let resp_path = fits_dir.join(format!("responsibilities_C{}.tsv", c));
        log_step(format!(
            "Writing responsibilities to {}",
            resp_path.display()
        ));
        em::write_responsibilities(&resp_path, &dataset, &result)?;
        summary.push(serde_json::json!({
            "classes": c,
            "log_likelihood": result.log_likelihood,
            "score": final_score,
            "model_path": model_path.file_name().unwrap().to_string_lossy(),
            "responsibilities_path": resp_path.file_name().unwrap().to_string_lossy()
        }));
        match &best {
            None => best = Some((c, result, final_score)),
            Some((_, _, best_score)) if final_score < *best_score => {
                best = Some((c, result, final_score))
            }
            _ => {}
        }
    }

    let (best_c, best_result, best_score) = best.ok_or_else(|| anyhow!("no fits ran"))?;
    log_step(format!(
        "Best C = {} ({} = {:.4}, log-likelihood = {:.4})",
        best_c, criterion_label, best_score, best_result.log_likelihood
    ));
    let best_model = em::ModelFile::from_parts(&config, &best_result.params);
    log_step(format!(
        "Writing best model and responsibilities to {}",
        args.out.display()
    ));
    write_model(&args.out.join("best_model.json"), &best_model)?;
    em::write_responsibilities(
        &args.out.join("best_responsibilities.tsv"),
        &dataset,
        &best_result,
    )?;

    // Impute with best model.
    log_step("Imputing with best model");
    let imputed = impute_dataset(&dataset, &best_result.params)?;
    let imputed_path = args.out.join("imputed.tsv");
    log_step(format!("Writing imputed dataset to {}", imputed_path.display()));
    data::write_dataset(&imputed_path, &imputed, args.delimiter_byte())?;
    write_imputed_labels(
        &args.out.join("imputed_labels.tsv"),
        &imputed,
        &subcat_names,
        args.delimiter_byte(),
        &best_result,
    )?;

    // Summary JSON.
    let summary_path = args.out.join("summary.json");
    log_step(format!(
        "Writing summary to {}",
        summary_path.display()
    ));
    let summary_obj = serde_json::json!({
        "min_classes": args.min_classes,
        "max_classes": args.max_classes,
        "criterion": args.criterion,
        "best_classes": best_c,
        "best_score": best_score,
        "fits": summary
    });
    std::fs::write(summary_path, serde_json::to_string_pretty(&summary_obj)?)?;
    log_step(format!(
        "Pipeline finished in {:.2?}",
        pipeline_start.elapsed()
    ));

    Ok(())
}

fn write_imputed_labels(
    path: &PathBuf,
    dataset: &Dataset,
    subcat_names: &[Vec<String>],
    delimiter: u8,
    result: &em::EmResult,
) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_path(path)?;
    let mut header: Vec<String> = vec!["id".to_string()];
    header.extend(dataset.category_names.iter().cloned());
    header.extend(dataset.methylation_names.iter().cloned());
    // responsibilities columns + argmax class
    for c in 0..result.params.pi.len() {
        header.push(format!("resp_class_{}", c));
    }
    header.push("assigned_class".to_string());
    writer.write_record(header)?;
    for (i, id) in dataset.ids.iter().enumerate() {
        let mut row = vec![id.clone()];
        for k in 0..dataset.n_categories() {
            let label = match dataset.data[(i, k)] {
                Some(idx) => subcat_names
                    .get(k)
                    .and_then(|v| v.get(idx))
                    .cloned()
                    .unwrap_or_else(|| idx.to_string()),
                None => "NA".to_string(),
            };
            row.push(label);
        }
        if let Some(meth) = &dataset.methylation {
            for m in 0..dataset.n_methylation() {
                let val = match meth[(i, m)] {
                    Some(v) => format!("{:.6}", v),
                    None => "NA".to_string(),
                };
                row.push(val);
            }
        }
        // responsibilities
        let mut max_c = 0;
        let mut max_resp = -1.0;
        for c in 0..result.params.pi.len() {
            let r = result.responsibilities[(i, c)];
            if r > max_resp {
                max_resp = r;
                max_c = c;
            }
            row.push(format!("{:.6}", r));
        }
        row.push(max_c.to_string());
        writer.write_record(&row)?;
    }
    writer.flush()?;
    Ok(())
}

fn describe_dataset(dataset: &Dataset) -> String {
    let total_cat = dataset.n_blocks() * dataset.n_categories();
    let missing_cat = dataset.data.iter().filter(|v| v.is_none()).count();
    let mut parts = vec![format!("blocks = {}", dataset.n_blocks())];
    parts.push(format!(
        "categorical dims = {} (levels = {:?}, missing {}/{})",
        dataset.n_categories(),
        dataset.n_levels,
        missing_cat,
        total_cat
    ));
    if dataset.n_methylation() > 0 {
        let missing_meth = dataset
            .methylation
            .as_ref()
            .map(|m| m.iter().filter(|v| v.is_none()).count())
            .unwrap_or(0);
        let total_meth = dataset.n_blocks() * dataset.n_methylation();
        parts.push(format!(
            "methylation dims = {} (missing {}/{})",
            dataset.n_methylation(),
            missing_meth,
            total_meth
        ));
    }
    format!("Dataset overview: {}", parts.join(", "))
}

fn log_step(msg: impl AsRef<str>) {
    eprintln!("[methylphase::typing] {}", msg.as_ref());
}

fn prepare_methylation(args: &RunArgs) -> Result<PathBuf> {
    let bam = args.bam.clone();
    if args.motifs.is_empty() && args.motif_file.is_none() {
        return Err(anyhow!(
            "you must specify --motif or --motif-file for split-reads"
        ));
    }
    let split_out = args.out.join("split_reads");
    std::fs::create_dir_all(&split_out)?;
    let clustering_path = split_out.join("read_clustering_raw.tsv");

    log_step(format!(
        "Generating methylation features from BAM {}; outputs under {}",
        bam.display(),
        split_out.display()
    ));
    let seq_fallback = args.sequence_fallback.clone();
    split_reads::run(
        bam,
        args.motifs.clone(),
        args.motif_file.clone(),
        seq_fallback,
        None,
        split_out.clone(),
        5,
        None,
        false,
        1,
        ClusterAlgorithm::Gmm,
        None,
        Vec::new(),
    )?;
    Ok(clustering_path)
}
