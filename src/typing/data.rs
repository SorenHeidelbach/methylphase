use crate::typing::error::MethylError;
use csv::ReaderBuilder;
use ndarray::Array2;
use serde::{Deserialize, Serialize};
use std::path::Path;

/// Configuration describing category names and levels.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct CategoryConfig {
    pub categories: Categories,
    /// Optional methylation category names (continuous features).
    #[serde(default)]
    pub methylation_names: Vec<String>,
}

/// Inner struct for TOML parsing.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Categories {
    pub names: Vec<String>,
    pub levels: Vec<usize>,
}

/// Dataset representation holding categorical observations with missing data.
#[derive(Debug, Clone)]
pub struct Dataset {
    pub ids: Vec<String>,
    /// Observations: shape [n_blocks][n_categories], values are Option<usize>.
    pub data: Array2<Option<usize>>,
    /// Number of subcategories per category.
    pub n_levels: Vec<usize>,
    /// Optional category names.
    pub category_names: Vec<String>,
    /// Continuous methylation values: shape [n_blocks][n_meth], Option<f64>.
    pub methylation: Option<Array2<Option<f64>>>,
    /// Names for methylation dimensions.
    pub methylation_names: Vec<String>,
}

impl Dataset {
    pub fn n_blocks(&self) -> usize {
        self.ids.len()
    }

    pub fn n_categories(&self) -> usize {
        self.n_levels.len()
    }

    pub fn n_methylation(&self) -> usize {
        self.methylation
            .as_ref()
            .map(|a| a.shape()[1])
            .unwrap_or(0)
    }

    /// Validate observations are consistent with category levels.
    pub fn validate(&self) -> Result<(), MethylError> {
        for ((i, k), value) in self.data.indexed_iter() {
            if let Some(v) = value {
                let max = self
                    .n_levels
                    .get(k)
                    .ok_or_else(|| MethylError::InconsistentCategory {
                        category: k,
                        message: "category index out of bounds".to_string(),
                    })?;
                if v >= max {
                    return Err(MethylError::InconsistentCategory {
                        category: k,
                        message: format!(
                            "value {} at row {} exceeds max level {}",
                            v, i, max
                        ),
                    });
                }
            }
        }
        if let Some(meth) = &self.methylation {
            if meth.shape()[0] != self.n_blocks() {
                return Err(MethylError::Shape(format!(
                    "methylation rows {} != ids {}",
                    meth.shape()[0],
                    self.n_blocks()
                )));
            }
            if meth.shape()[1] != self.methylation_names.len() {
                return Err(MethylError::Shape(format!(
                    "methylation cols {} != methylation_names {}",
                    meth.shape()[1],
                    self.methylation_names.len()
                )));
            }
        }
        Ok(())
    }

    /// Create a new dataset with only the given row indices (order preserved).
    pub fn select_rows(&self, indices: &[usize]) -> Result<Dataset, MethylError> {
        let mut ids = Vec::with_capacity(indices.len());
        let mut data = Array2::<Option<usize>>::from_elem((indices.len(), self.n_categories()), None);
        for (row_idx, &src) in indices.iter().enumerate() {
            let id = self
                .ids
                .get(src)
                .ok_or_else(|| MethylError::Shape(format!("row {} out of bounds", src)))?;
            ids.push(id.clone());
            for k in 0..self.n_categories() {
                data[(row_idx, k)] = *self
                    .data
                    .get((src, k))
                    .ok_or_else(|| MethylError::Shape(format!("row {} col {} out of bounds", src, k)))?;
            }
        }

        let methylation = if let Some(meth) = &self.methylation {
            let mut out = Array2::<Option<f64>>::from_elem((indices.len(), self.n_methylation()), None);
            for (row_idx, &src) in indices.iter().enumerate() {
                for m in 0..self.n_methylation() {
                    out[(row_idx, m)] = *meth
                        .get((src, m))
                        .ok_or_else(|| MethylError::Shape(format!("row {} meth {} out of bounds", src, m)))?;
                }
            }
            Some(out)
        } else {
            None
        };

        Ok(Dataset {
            ids,
            data,
            n_levels: self.n_levels.clone(),
            category_names: self.category_names.clone(),
            methylation,
            methylation_names: self.methylation_names.clone(),
        })
    }
}

/// Load TOML configuration describing categories.
pub fn load_config(path: &Path) -> Result<CategoryConfig, MethylError> {
    let contents = std::fs::read_to_string(path)?;
    let config: CategoryConfig = toml::from_str(&contents)?;
    if config.categories.names.len() != config.categories.levels.len() {
        return Err(MethylError::Config(format!(
            "names and levels length mismatch ({} vs {})",
            config.categories.names.len(),
            config.categories.levels.len()
        )));
    }
    Ok(config)
}

/// Load dataset from TSV/CSV with NA representing missing values.
pub fn load_dataset(
    path: &Path,
    config: &CategoryConfig,
    delimiter: u8,
) -> Result<Dataset, MethylError> {
    let mut reader = ReaderBuilder::new()
        .delimiter(delimiter)
        .flexible(true)
        .from_path(path)?;

    let mut ids = Vec::new();
    let mut rows: Vec<Vec<Option<usize>>> = Vec::new();
    for record in reader.records() {
        let record = record?;
        if record.is_empty() {
            continue;
        }
        let id = record
            .get(0)
            .ok_or_else(|| MethylError::Parse("missing id column".to_string()))?
            .to_string();
        ids.push(id);
        let mut row = Vec::new();
        for field in record.iter().skip(1) {
            if field == "NA" {
                row.push(None);
            } else {
                let val: usize = field.parse().map_err(|e| {
                    MethylError::Parse(format!("failed to parse value {}: {}", field, e))
                })?;
                row.push(Some(val));
            }
        }
        rows.push(row);
    }

    let n_categories = config.categories.levels.len();
    let n_rows = rows.len();
    let mut data = Array2::<Option<usize>>::from_elem((n_rows, n_categories), None);
    for (i, row) in rows.iter().enumerate() {
        if row.len() != n_categories {
            return Err(MethylError::Parse(format!(
                "row {} has {} categories, expected {}",
                i,
                row.len(),
                n_categories
            )));
        }
        for (k, value) in row.iter().enumerate() {
            data[(i, k)] = *value;
        }
    }

    let dataset = Dataset {
        ids,
        data,
        n_levels: config.categories.levels.clone(),
        category_names: config.categories.names.clone(),
        methylation: None,
        methylation_names: config.methylation_names.clone(),
    };
    dataset.validate()?;
    Ok(dataset)
}

/// Write dataset to TSV/CSV, converting None to "NA".
pub fn write_dataset(
    path: &Path,
    dataset: &Dataset,
    delimiter: u8,
) -> Result<(), MethylError> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_path(path)?;

    let mut header: Vec<String> = vec!["id".to_string()];
    header.extend(dataset.category_names.iter().cloned());
    header.extend(dataset.methylation_names.iter().cloned());
    writer.write_record(header)?;

    for (i, id) in dataset.ids.iter().enumerate() {
        let mut record = vec![id.clone()];
        for k in 0..dataset.n_categories() {
            let val = match dataset.data[(i, k)] {
                Some(v) => v.to_string(),
                None => "NA".to_string(),
            };
            record.push(val);
        }
        if let Some(meth) = &dataset.methylation {
            for m in 0..dataset.n_methylation() {
                let val = match meth[(i, m)] {
                    Some(v) => format!("{:.6}", v),
                    None => "NA".to_string(),
                };
                record.push(val);
            }
        }
        writer.write_record(&record)?;
    }
    writer.flush()?;
    Ok(())
}
