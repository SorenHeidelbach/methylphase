use crate::typing::data::Dataset;
use crate::typing::em::EmParams;
use crate::typing::error::MethylError;
use ndarray::Array2;
use std::f64::NEG_INFINITY;

/// Impute missing categorical values using MAP with fitted parameters.
pub fn impute_dataset(dataset: &Dataset, params: &EmParams) -> Result<Dataset, MethylError> {
    let responsibilities = compute_responsibilities(dataset, params)?;
    let mut imputed = dataset.clone();
    for i in 0..dataset.n_blocks() {
        for k in 0..dataset.n_categories() {
            if dataset.data[(i, k)].is_none() {
                let m_k = dataset.n_levels[k];
                let mut best_j = 0;
                let mut best_p = -1.0;
                for j in 0..m_k {
                    let mut prob = 0.0;
                    for c in 0..params.pi.len() {
                        prob += responsibilities[(i, c)] * params.phi[c][k][j];
                    }
                    if prob > best_p {
                        best_p = prob;
                        best_j = j;
                    }
                }
                imputed.data[(i, k)] = Some(best_j);
            }
        }
    }
    // Impute continuous methylation using per-target responsibilities that condition on
    // all observed features except the missing methyl dimension.
    if let Some(meth) = &dataset.methylation {
        let mut new_meth = meth.clone();
        for i in 0..dataset.n_blocks() {
            for m in 0..dataset.n_methylation() {
                if meth[(i, m)].is_none() {
                    let resp_m = compute_row_responsibilities(dataset, params, i, Some(m))?;
                    let mut num = 0.0;
                    for (c, r) in resp_m.iter().enumerate() {
                        num += *r * params.meth_mean[c][m];
                    }
                    new_meth[(i, m)] = Some(num);
                }
            }
        }
        imputed.methylation = Some(new_meth);
    }
    Ok(imputed)
}

fn compute_responsibilities(
    dataset: &Dataset,
    params: &EmParams,
) -> Result<Array2<f64>, MethylError> {
    let n_classes = params.pi.len();
    let mut responsibilities = Array2::<f64>::zeros((dataset.n_blocks(), n_classes));
    for i in 0..dataset.n_blocks() {
        let resp_row = compute_row_responsibilities(dataset, params, i, None)?;
        for (c, r) in resp_row.iter().enumerate() {
            responsibilities[(i, c)] = *r;
        }
    }
    Ok(responsibilities)
}

fn compute_row_responsibilities(
    dataset: &Dataset,
    params: &EmParams,
    row: usize,
    exclude_methyl: Option<usize>,
) -> Result<Vec<f64>, MethylError> {
    let n_classes = params.pi.len();
    let mut log_scores = vec![0.0; n_classes];
    for c in 0..n_classes {
        let mut log_score = params.pi[c].ln();
        for k in 0..dataset.n_categories() {
            if let Some(j) = dataset.data[(row, k)] {
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
                if Some(m) == exclude_methyl {
                    continue;
                }
                if let Some(v) = meth[(row, m)] {
                    let mean = params.meth_mean[c][m];
                    let var = params.meth_var[c][m].max(1e-12);
                    let log_pdf = -0.5 * ((2.0 * std::f64::consts::PI * var).ln())
                        - (v - mean).powi(2) / (2.0 * var);
                    log_score += log_pdf;
                }
            }
        }
        log_scores[c] = log_score;
    }
    let max_log = log_scores.iter().cloned().fold(NEG_INFINITY, f64::max);
    let mut sum = 0.0;
    let mut resp = vec![0.0; n_classes];
    for c in 0..n_classes {
        let val = (log_scores[c] - max_log).exp();
        resp[c] = val;
        sum += val;
    }
    if sum == 0.0 {
        return Err(MethylError::Math(
            "responsibility normalization underflow".into(),
        ));
    }
    for r in &mut resp {
        *r /= sum;
    }
    Ok(resp)
}
