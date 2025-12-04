use serde::{Deserialize, Serialize};

/// Dirichlet hyperparameters applied to mixture weights and categorical emissions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DirichletPriors {
    pub alpha_pi: f64,
    pub alpha_phi: f64,
}

impl DirichletPriors {
    pub fn validate(&self) -> Result<(), String> {
        if self.alpha_pi <= 0.0 || self.alpha_phi <= 0.0 {
            return Err("Dirichlet alphas must be positive".to_string());
        }
        Ok(())
    }
}
