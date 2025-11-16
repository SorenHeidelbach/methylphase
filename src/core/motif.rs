use anyhow::{bail, Context, Result};

#[derive(Debug, Clone)]
pub struct MotifQuery {
    raw: String,
    pattern: Vec<u8>,
    masks: Vec<u8>,
    mod_label: String,
    mod_index: usize,
}

impl MotifQuery {
    pub fn parse(spec: &str) -> Result<Self> {
        let parts: Vec<_> = spec.split('_').collect();
        if parts.len() != 3 {
            bail!("motif specification `{spec}` must be motif_modtype_modposition");
        }

        let motif = parts[0].to_uppercase();
        if motif.is_empty() {
            bail!("motif portion may not be empty");
        }

        let mod_label = parts[1].to_string();
        let mod_position: usize = parts[2]
            .parse()
            .with_context(|| format!("invalid modification index in `{spec}`"))?;
        if mod_position == 0 || mod_position > motif.len() {
            bail!(
                "modification index must be within motif bounds (1..={}), got {}",
                motif.len(),
                mod_position
            );
        }

        let masks = motif
            .chars()
            .map(|ch| {
                base_mask_from_symbol(ch).with_context(|| format!("invalid motif base `{ch}`"))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Self {
            raw: spec.to_string(),
            pattern: motif.into_bytes(),
            masks,
            mod_label,
            mod_index: mod_position - 1,
        })
    }

    pub fn raw(&self) -> &str {
        &self.raw
    }

    pub fn motif(&self) -> &[u8] {
        &self.pattern
    }

    pub fn mod_label(&self) -> &str {
        &self.mod_label
    }

    pub fn mod_position(&self) -> usize {
        self.mod_index + 1
    }

    pub fn matches<'a>(&self, sequence: &'a [u8]) -> Vec<usize> {
        let motif_len = self.pattern.len();
        if motif_len == 0 || sequence.len() < motif_len {
            return Vec::new();
        }

        sequence
            .windows(motif_len)
            .enumerate()
            .filter_map(|(idx, window)| {
                if window
                    .iter()
                    .zip(&self.masks)
                    .all(|(base, mask)| base_matches(*base, *mask))
                {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn modification_position_for_match(&self, match_start: usize) -> usize {
        match_start + self.mod_index
    }
}

const A_BIT: u8 = 0b0001;
const C_BIT: u8 = 0b0010;
const G_BIT: u8 = 0b0100;
const T_BIT: u8 = 0b1000;

fn base_mask_from_symbol(symbol: char) -> Result<u8> {
    use anyhow::bail;
    let mask = match symbol {
        'A' => A_BIT,
        'C' => C_BIT,
        'G' => G_BIT,
        'T' | 'U' => T_BIT,
        'R' => A_BIT | G_BIT,
        'Y' => C_BIT | T_BIT,
        'S' => G_BIT | C_BIT,
        'W' => A_BIT | T_BIT,
        'K' => G_BIT | T_BIT,
        'M' => A_BIT | C_BIT,
        'B' => C_BIT | G_BIT | T_BIT,
        'D' => A_BIT | G_BIT | T_BIT,
        'H' => A_BIT | C_BIT | T_BIT,
        'V' => A_BIT | C_BIT | G_BIT,
        'N' => A_BIT | C_BIT | G_BIT | T_BIT,
        _ => bail!("unsupported motif character `{symbol}`"),
    };

    Ok(mask)
}

fn base_matches(base: u8, mask: u8) -> bool {
    base_mask_from_base(base)
        .map(|b| b & mask != 0)
        .unwrap_or(false)
}

fn base_mask_from_base(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(A_BIT),
        b'C' => Some(C_BIT),
        b'G' => Some(G_BIT),
        b'T' | b'U' => Some(T_BIT),
        b'N' => Some(A_BIT | C_BIT | G_BIT | T_BIT),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_degenerate_motif_matches() {
        let motif = MotifQuery::parse("RGATCY_6mA_2").unwrap();
        let hits = motif.matches(b"AGATCT");
        assert_eq!(hits, vec![0]);

        assert!(motif.matches(b"TGATCA").is_empty());
    }
}
