use crate::typing::data::{CategoryConfig, Dataset};
use crate::typing::error::MethylError;
use csv::ReaderBuilder;
use ndarray::Array2;
use std::collections::HashMap;
use std::path::Path;

/// Load methylation file with columns: read_id, cluster_id, motif1, motif2, ...
/// cluster_id -1 is treated as unassigned; motif values are continuous per motif.
pub fn load_methylation_features(
    path: &Path,
    delimiter: u8,
) -> Result<(HashMap<String, Vec<Option<f64>>>, Vec<String>), MethylError> {
    let mut reader = ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(true)
        .from_path(path)?;

    let headers = reader.headers()?.clone();
    if headers.len() < 3 {
        return Err(MethylError::Parse(
            "methylation file must have at least read_id, cluster_id, and one motif column".into(),
        ));
    }
    let motif_names: Vec<String> = headers.iter().skip(2).map(|s| s.to_string()).collect();

    let mut map: HashMap<String, Vec<Option<f64>>> = HashMap::new();
    for record in reader.records() {
        let record = record?;
        let id = record
            .get(0)
            .ok_or_else(|| MethylError::Parse("missing read_id".into()))?
            .to_string();
        let mut vals: Vec<Option<f64>> = Vec::with_capacity(motif_names.len());
        for field in record.iter().skip(2) {
            if field == "NA" || field.is_empty() {
                vals.push(None);
            } else {
                let v: f64 = field.parse().map_err(|e| {
                    MethylError::Parse(format!("invalid methylation value {}: {}", field, e))
                })?;
                vals.push(Some(v));
            }
        }
        map.insert(id, vals);
    }
    Ok((map, motif_names))
}

/// Attach continuous methylation features: one column per motif, value = methylation score.
/// Missing reads or missing values become NA.
pub fn add_methylation_features(
    dataset: &Dataset,
    config: &CategoryConfig,
    values: &HashMap<String, Vec<Option<f64>>>,
    motif_names: &[String],
) -> Result<(Dataset, CategoryConfig), MethylError> {
    let existing: std::collections::HashSet<&String> = dataset.ids.iter().collect();
    let mut extra_ids: Vec<String> = values
        .keys()
        .filter(|id| !existing.contains(id))
        .cloned()
        .collect();
    extra_ids.sort_unstable();

    let mut all_ids = dataset.ids.clone();
    all_ids.extend(extra_ids.into_iter());

    let n_rows = all_ids.len();
    let n_categories = dataset.n_categories();

    let mut data = Array2::<Option<usize>>::from_elem((n_rows, n_categories), None);
    // Copy existing categorical observations; extra rows remain None (no hap evidence).
    for ((i, k), v) in dataset.data.indexed_iter() {
        data[(i, k)] = *v;
    }

    let mut methylation = Array2::<Option<f64>>::from_elem((n_rows, motif_names.len()), None);

    for (row, id) in all_ids.iter().enumerate() {
        if let Some(vals) = values.get(id) {
            if vals.len() != motif_names.len() {
                return Err(MethylError::Parse(format!(
                    "motif count mismatch for id {}: {} vs {}",
                    id,
                    vals.len(),
                    motif_names.len()
                )));
            }
            for (m, v) in vals.iter().enumerate() {
                methylation[(row, m)] = *v;
            }
        }
    }

    let new_dataset = Dataset {
        ids: all_ids,
        data,
        n_levels: config.categories.levels.clone(),
        category_names: config.categories.names.clone(),
        methylation: Some(methylation),
        methylation_names: motif_names.to_vec(),
    };
    new_dataset.validate()?;

    let mut new_config = config.clone();
    new_config.methylation_names = motif_names.to_vec();

    Ok((new_dataset, new_config))
}
