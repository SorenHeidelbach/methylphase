use crate::typing::data::Dataset;
use crate::typing::error::MethylError;
use crate::typing::priors::DirichletPriors;
use ndarray::Array2;
use rand::distributions::Distribution;
use rand::thread_rng;
use rand_distr::Dirichlet;
use serde::{Deserialize, Serialize};
use std::f64::NEG_INFINITY;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Mixture weights and categorical emission parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmParams {
    pub pi: Vec<f64>,                // size C
    pub phi: Vec<Vec<Vec<f64>>>,     // [c][k][j] categorical
    /// Gaussian means per class per methylation dimension (if any).
    pub meth_mean: Vec<Vec<f64>>,    // [c][m]
    /// Gaussian variances per class per methylation dimension (if any).
    pub meth_var: Vec<Vec<f64>>,     // [c][m]
}

/// EM algorithm settings.
#[derive(Debug, Clone, Copy)]
pub struct EmSettings {
    pub max_iters: usize,
    pub tol: f64,
    /// Minimum effective mixture weight (responsibility mean) to keep a class.
    pub min_class_weight: f64,
}

/// Result of EM fitting.
#[derive(Debug, Clone)]
pub struct EmResult {
    pub params: EmParams,
    /// responsibilities tau[i][c]
    pub responsibilities: Array2<f64>,
    pub log_likelihood: f64,
    pub iterations: usize,
}

/// Model serialization container, capturing parameters and metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ModelFile {
    pub params: EmParams,
    pub n_levels: Vec<usize>,
    pub category_names: Vec<String>,
    pub methylation_names: Vec<String>,
}

impl ModelFile {
    pub fn from_parts(config: &crate::typing::data::CategoryConfig, params: &EmParams) -> Self {
        Self {
            params: params.clone(),
            n_levels: config.categories.levels.clone(),
            category_names: config.categories.names.clone(),
            methylation_names: config.methylation_names.clone(),
        }
    }
}

/// Fit EM for a fixed number of latent classes.
pub fn fit_em(
    dataset: &Dataset,
    n_classes: usize,
    priors: &DirichletPriors,
    settings: EmSettings,
) -> Result<EmResult, MethylError> {
    priors
        .validate()
        .map_err(|e| MethylError::Config(e.to_string()))?;
    if n_classes == 0 {
        return Err(MethylError::Config("n_classes must be > 0".into()));
    }
    if settings.min_class_weight < 0.0 || settings.min_class_weight >= 1.0 {
        return Err(MethylError::Config(
            "min_class_weight must be in [0, 1)".into(),
        ));
    }
    let mut params = initialize_params(dataset, n_classes, priors)?;
    let mut responsibilities = Array2::<f64>::zeros((dataset.n_blocks(), n_classes));
    let mut prev_ll = NEG_INFINITY;
    let mut current_ll = NEG_INFINITY;
    let mut iter = 0;

    while iter < settings.max_iters {
        e_step(dataset, &params, &mut responsibilities)?;
        m_step(dataset, &responsibilities, priors, &mut params)?;
        current_ll = log_likelihood(dataset, &params)?;
        if (current_ll - prev_ll).abs() < settings.tol {
            break;
        }
        prev_ll = current_ll;
        iter += 1;
    }

    if settings.min_class_weight > 0.0 {
        prune_small_classes(
            dataset,
            &mut params,
            &mut responsibilities,
            settings.min_class_weight,
        )?;
        current_ll = log_likelihood(dataset, &params)?;
    }

    Ok(EmResult {
        params,
        responsibilities,
        log_likelihood: current_ll,
        iterations: iter + 1,
    })
}

fn initialize_params(
    dataset: &Dataset,
    n_classes: usize,
    priors: &DirichletPriors,
) -> Result<EmParams, MethylError> {
    let mut rng = thread_rng();
    let pi = vec![1.0 / n_classes as f64; n_classes];

    let mut phi = Vec::with_capacity(n_classes);
    for _c in 0..n_classes {
        let mut class_phi = Vec::with_capacity(dataset.n_categories());
        for &levels in &dataset.n_levels {
            let alpha = vec![priors.alpha_phi; levels];
            let dirichlet = Dirichlet::new(&alpha).map_err(|e| {
                MethylError::Math(format!("invalid Dirichlet parameters: {}", e))
            })?;
            let sample = dirichlet.sample(&mut rng);
            class_phi.push(sample);
        }
        phi.push(class_phi);
    }

    // Initialize methylation Gaussian params.
    let m_dims = dataset.n_methylation();
    let mut meth_mean = vec![vec![0.0; m_dims]; n_classes];
    let mut meth_var = vec![vec![1.0; m_dims]; n_classes];
    if let Some(meth) = &dataset.methylation {
        for m in 0..m_dims {
            let mut vals = Vec::new();
            for i in 0..dataset.n_blocks() {
                if let Some(v) = meth[(i, m)] {
                    vals.push(v);
                }
            }
            let mean = if vals.is_empty() {
                0.0
            } else {
                vals.iter().copied().sum::<f64>() / vals.len() as f64
            };
            let var = if vals.len() <= 1 {
                1.0
            } else {
                let ss: f64 = vals.iter().map(|v| (v - mean).powi(2)).sum();
                (ss / (vals.len() as f64 - 1.0)).max(1e-6)
            };
            for c in 0..n_classes {
                meth_mean[c][m] = mean;
                meth_var[c][m] = var;
            }
        }
    }

    Ok(EmParams { pi, phi, meth_mean, meth_var })
}

fn e_step(
    dataset: &Dataset,
    params: &EmParams,
    responsibilities: &mut Array2<f64>,
) -> Result<(), MethylError> {
    let n_classes = params.pi.len();
    let n_blocks = dataset.n_blocks();
    for i in 0..n_blocks {
        let mut log_scores = vec![0.0; n_classes];
        for c in 0..n_classes {
            let mut log_score = params.pi[c].ln();
            for k in 0..dataset.n_categories() {
                if let Some(j) = dataset.data[(i, k)] {
                    let prob = *params
                        .phi
                        .get(c)
                        .and_then(|cats| cats.get(k))
                        .and_then(|levels| levels.get(j))
                        .ok_or_else(|| {
                            MethylError::Shape(format!(
                                "missing phi for class {}, category {}, level {}",
                                c, k, j
                            ))
                        })?;
                    if prob <= 0.0 {
                        return Err(MethylError::Math(
                            "encountered zero probability in phi".into(),
                        ));
                    }
                    log_score += prob.ln();
                }
            }
            if let Some(meth) = &dataset.methylation {
                for m in 0..dataset.n_methylation() {
                    if let Some(v) = meth[(i, m)] {
                        let mean = params.meth_mean[c][m];
                        let var = params.meth_var[c][m].max(1e-12);
                        let log_pdf =
                            -0.5 * ((2.0 * std::f64::consts::PI * var).ln())
                                - (v - mean).powi(2) / (2.0 * var);
                        log_score += log_pdf;
                    }
                }
            }
            log_scores[c] = log_score;
        }
        let max_log = log_scores
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let mut sum = 0.0;
        for c in 0..n_classes {
            let val = (log_scores[c] - max_log).exp();
            responsibilities[(i, c)] = val;
            sum += val;
        }
        if sum == 0.0 {
            return Err(MethylError::Math(
                "responsibility normalization underflow".into(),
            ));
        }
        for c in 0..n_classes {
            responsibilities[(i, c)] /= sum;
        }
    }
    Ok(())
}

fn m_step(
    dataset: &Dataset,
    responsibilities: &Array2<f64>,
    priors: &DirichletPriors,
    params: &mut EmParams,
) -> Result<(), MethylError> {
    let n_classes = params.pi.len();
    let n_blocks = dataset.n_blocks();
    // Update pi
    let mut class_counts = vec![0.0; n_classes];
    for c in 0..n_classes {
        class_counts[c] = responsibilities.column(c).sum();
    }
    let total = n_blocks as f64 + priors.alpha_pi * n_classes as f64;
    for c in 0..n_classes {
        params.pi[c] = (class_counts[c] + priors.alpha_pi) / total;
    }

    // Update phi
    for c in 0..n_classes {
        for k in 0..dataset.n_categories() {
            let levels = dataset.n_levels[k];
            let mut counts = vec![0.0; levels];
            for i in 0..n_blocks {
                if let Some(j) = dataset.data[(i, k)] {
                    counts[j] += responsibilities[(i, c)];
                }
            }
            let norm: f64 = counts.iter().sum::<f64>() + priors.alpha_phi * levels as f64;
            for j in 0..levels {
                params.phi[c][k][j] = (counts[j] + priors.alpha_phi) / norm;
            }
        }
    }

    // Update methylation Gaussian parameters
    if let Some(meth) = &dataset.methylation {
        for c in 0..n_classes {
            for m in 0..dataset.n_methylation() {
                let mut weight_sum = 0.0;
                let mut mean_num = 0.0;
                for i in 0..n_blocks {
                    if let Some(v) = meth[(i, m)] {
                        let w = responsibilities[(i, c)];
                        weight_sum += w;
                        mean_num += w * v;
                    }
                }
                let mean = if weight_sum > 0.0 {
                    mean_num / weight_sum
                } else {
                    0.0
                };
                let mut var_num = 0.0;
                for i in 0..n_blocks {
                    if let Some(v) = meth[(i, m)] {
                        let w = responsibilities[(i, c)];
                        var_num += w * (v - mean).powi(2);
                    }
                }
                let var = if weight_sum > 1e-9 {
                    (var_num / weight_sum).max(1e-6)
                } else {
                    1.0
                };
                params.meth_mean[c][m] = mean;
                params.meth_var[c][m] = var;
            }
        }
    }
    Ok(())
}

fn prune_small_classes(
    dataset: &Dataset,
    params: &mut EmParams,
    responsibilities: &mut Array2<f64>,
    min_weight: f64,
) -> Result<(), MethylError> {
    let n_classes = params.pi.len();
    let n_blocks = dataset.n_blocks();
    if n_classes == 0 {
        return Err(MethylError::Shape("no classes to prune".into()));
    }
    let effective: Vec<f64> = (0..n_classes)
        .map(|c| responsibilities.column(c).sum() / n_blocks as f64)
        .collect();
    let mut keep: Vec<usize> = effective
        .iter()
        .enumerate()
        .filter(|(_, w)| **w >= min_weight)
        .map(|(idx, _)| idx)
        .collect();
    if keep.is_empty() {
        if let Some((idx, _)) = effective
            .iter()
            .enumerate()
            .max_by(|a, b| a
                .1
                .partial_cmp(b.1)
                .unwrap_or(std::cmp::Ordering::Equal))
        {
            keep.push(idx);
        } else {
            return Err(MethylError::Other(
                "unable to compute effective weights during pruning".into(),
            ));
        }
    }
    if keep.len() == n_classes {
        return Ok(());
    }

    let mut new_params = EmParams {
        pi: Vec::with_capacity(keep.len()),
        phi: Vec::with_capacity(keep.len()),
        meth_mean: Vec::with_capacity(keep.len()),
        meth_var: Vec::with_capacity(keep.len()),
    };
    for &idx in &keep {
        new_params.phi.push(
            params
                .phi
                .get(idx)
                .cloned()
                .ok_or_else(|| MethylError::Shape("phi missing during prune".into()))?,
        );
        new_params.meth_mean.push(
            params
                .meth_mean
                .get(idx)
                .cloned()
                .ok_or_else(|| MethylError::Shape("meth_mean missing during prune".into()))?,
        );
        new_params.meth_var.push(
            params
                .meth_var
                .get(idx)
                .cloned()
                .ok_or_else(|| MethylError::Shape("meth_var missing during prune".into()))?,
        );
    }

    let mut new_pi: Vec<f64> = keep
        .iter()
        .map(|idx| effective[*idx].max(1e-12))
        .collect();
    let pi_sum: f64 = new_pi.iter().sum();
    for val in &mut new_pi {
        *val /= pi_sum;
    }
    new_params.pi = new_pi;

    let mut new_resp = Array2::<f64>::zeros((n_blocks, keep.len()));
    for i in 0..n_blocks {
        let mut row_sum = 0.0;
        for (new_idx, old_idx) in keep.iter().enumerate() {
            let val = responsibilities[(i, *old_idx)];
            new_resp[(i, new_idx)] = val;
            row_sum += val;
        }
        if row_sum == 0.0 {
            let uniform = 1.0 / keep.len() as f64;
            for new_idx in 0..keep.len() {
                new_resp[(i, new_idx)] = uniform;
            }
        } else {
            for new_idx in 0..keep.len() {
                new_resp[(i, new_idx)] /= row_sum;
            }
        }
    }

    *params = new_params;
    *responsibilities = new_resp;
    Ok(())
}

pub(crate) fn log_likelihood(dataset: &Dataset, params: &EmParams) -> Result<f64, MethylError> {
    let n_classes = params.pi.len();
    let mut total = 0.0;
    for i in 0..dataset.n_blocks() {
        let mut log_scores = vec![0.0; n_classes];
        for c in 0..n_classes {
            let mut log_score = params.pi[c].ln();
            for k in 0..dataset.n_categories() {
                if let Some(j) = dataset.data[(i, k)] {
                    let prob = params.phi[c][k][j];
                    if prob <= 0.0 {
                        return Err(MethylError::Math(
                            "encountered zero probability in phi".into(),
                        ));
                    }
                    log_score += prob.ln();
                }
            }
            if let Some(meth) = &dataset.methylation {
                for m in 0..dataset.n_methylation() {
                    if let Some(v) = meth[(i, m)] {
                        let mean = params.meth_mean[c][m];
                        let var = params.meth_var[c][m].max(1e-12);
                        let log_pdf =
                            -0.5 * ((2.0 * std::f64::consts::PI * var).ln())
                                - (v - mean).powi(2) / (2.0 * var);
                        log_score += log_pdf;
                    }
                }
            }
            log_scores[c] = log_score;
        }
        let max_log = log_scores
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        let sum: f64 = log_scores
            .iter()
            .map(|ls| (ls - max_log).exp())
            .sum();
        total += max_log + sum.ln();
    }
    Ok(total)
}

/// Write responsibilities to a CSV/TSV file (id + per-class probabilities).
pub fn write_responsibilities(
    path: &Path,
    dataset: &Dataset,
    result: &EmResult,
) -> Result<(), MethylError> {
    let mut writer = File::create(path)?;
    let header: Vec<String> = std::iter::once("id".to_string())
        .chain((0..result.params.pi.len()).map(|c| format!("class_{}", c)))
        .collect();
    writer.write_all(header.join("\t").as_bytes())?;
    writer.write_all(b"\n")?;
    for (i, id) in dataset.ids.iter().enumerate() {
        let mut row = vec![id.clone()];
        for c in 0..result.params.pi.len() {
            row.push(format!("{:.6}", result.responsibilities[(i, c)]));
        }
        writer.write_all(row.join("\t").as_bytes())?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

/// Number of free parameters for BIC/AIC style metrics.
pub fn parameter_count(n_classes: usize, n_levels: &[usize], n_methylation: usize) -> usize {
    let pi_params = if n_classes > 0 { n_classes - 1 } else { 0 };
    let phi_params: usize = n_classes
        * n_levels
            .iter()
            .map(|&m_k| if m_k > 0 { m_k - 1 } else { 0 })
            .sum::<usize>();
    let meth_params = n_classes * n_methylation * 2; // mean and variance per dim
    pi_params + phi_params + meth_params
}
