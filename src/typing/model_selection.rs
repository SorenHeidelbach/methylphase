use crate::typing::data::Dataset;
use crate::typing::em::{parameter_count, EmResult, EmSettings};
use crate::typing::error::MethylError;
use crate::typing::priors::DirichletPriors;

/// Strategy trait to choose the best latent class count.
pub trait ModelSelectionStrategy {
    fn select_best_model<F>(
        &self,
        dataset: &Dataset,
        fit_fn: F,
    ) -> Result<(usize, EmResult, f64), MethylError>
    where
        F: Fn(&Dataset, usize) -> Result<EmResult, MethylError>;
}

/// Available criteria for model selection.
#[derive(Debug, Clone, Copy)]
pub enum Criterion {
    Bic,
    Icl,
    CrossValidated,
}

/// Selection over a range of class counts using a criterion.
#[derive(Debug, Clone)]
pub struct SelectionStrategy {
    pub min_classes: usize,
    pub max_classes: usize,
    pub priors: DirichletPriors,
    pub settings: EmSettings,
    pub criterion: Criterion,
    pub penalty_multiplier: f64,
    pub cv_folds: usize,
}

impl SelectionStrategy {
    pub fn score(&self, dataset: &Dataset, result: &EmResult, n_classes: usize) -> f64 {
        let n = dataset.n_blocks() as f64;
        let p = parameter_count(n_classes, &dataset.n_levels, dataset.n_methylation()) as f64;
        let bic = -2.0 * result.log_likelihood + self.penalty_multiplier * p * n.ln();
        match self.criterion {
            Criterion::Bic => bic,
            Criterion::Icl => {
                let entropy = -result
                    .responsibilities
                    .iter()
                    .filter(|v| **v > 0.0)
                    .map(|v| v * v.ln())
                    .sum::<f64>();
                bic + 2.0 * entropy
            }
            Criterion::CrossValidated => {
                panic!("cross-validated score requires fold evaluation")
            }
        }
    }

    pub fn cross_validated_score<F>(
        &self,
        dataset: &Dataset,
        n_classes: usize,
        fit_fn: &F,
    ) -> Result<f64, MethylError>
    where
        F: Fn(&Dataset, usize) -> Result<EmResult, MethylError>,
    {
        let n = dataset.n_blocks();
        if n < 2 {
            return Err(MethylError::Config(
                "cross-validation requires at least 2 rows".into(),
            ));
        }
        let folds = self.cv_folds.max(2).min(n);
        let mut fold_indices = vec![Vec::new(); folds];
        for idx in 0..n {
            fold_indices[idx % folds].push(idx);
        }

        let mut total_ll = 0.0;
        let mut total_rows = 0usize;
        for fold in &fold_indices {
            if fold.is_empty() {
                continue;
            }
            let mut train_indices = Vec::with_capacity(n - fold.len());
            for i in 0..n {
                if !fold.contains(&i) {
                    train_indices.push(i);
                }
            }
            let train = dataset.select_rows(&train_indices)?;
            let test = dataset.select_rows(fold)?;
            let fit = fit_fn(&train, n_classes)?;
            let ll = crate::typing::em::log_likelihood(&test, &fit.params)?;
            total_ll += ll;
            total_rows += test.n_blocks();
        }
        if total_rows == 0 {
            return Err(MethylError::Other(
                "cross-validation produced zero test rows".into(),
            ));
        }
        Ok(-total_ll / total_rows as f64)
    }
}

impl ModelSelectionStrategy for SelectionStrategy {
    fn select_best_model<F>(
        &self,
        dataset: &Dataset,
        fit_fn: F,
    ) -> Result<(usize, EmResult, f64), MethylError>
    where
        F: Fn(&Dataset, usize) -> Result<EmResult, MethylError>,
    {
        if self.min_classes == 0 || self.max_classes < self.min_classes {
            return Err(MethylError::Config(
                "invalid class range for selection".into(),
            ));
        }
        match self.criterion {
            Criterion::CrossValidated => {
                let mut best_cv: Option<(usize, f64)> = None;
                for c in self.min_classes..=self.max_classes {
                    let score = self.cross_validated_score(dataset, c, &fit_fn)?;
                    match &best_cv {
                        None => best_cv = Some((c, score)),
                        Some((_, best_score)) if score < *best_score => {
                            best_cv = Some((c, score))
                        }
                        _ => {}
                    }
                }
                let (best_c, best_score) = best_cv
                    .ok_or_else(|| MethylError::Other("no model fitted".into()))?;
                let best_fit = fit_fn(dataset, best_c)?;
                Ok((best_c, best_fit, best_score))
            }
            _ => {
                let mut best: Option<(usize, EmResult, f64)> = None;
                for c in self.min_classes..=self.max_classes {
                    let result = fit_fn(dataset, c)?;
                    let score = self.score(dataset, &result, c);
                    match &best {
                        None => best = Some((c, result, score)),
                        Some((_, _, best_score)) if score < *best_score => {
                            best = Some((c, result, score))
                        }
                        _ => {}
                    }
                }
                best.ok_or_else(|| MethylError::Other("no model fitted".into()))
            }
        }
    }
}
