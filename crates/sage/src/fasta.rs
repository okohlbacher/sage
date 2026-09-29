use crate::enzyme::{Digest, EnzymeParameters};
use rayon::prelude::*;
use std::sync::Arc;

#[derive(Clone)]
pub struct Fasta {
    pub targets: Vec<(Arc<str>, String)>,
    decoy_tag: String,
    // Should we ignore decoys in the fasta database
    // and generate them internally?
    generate_decoys: bool,
}

/// First word of a header line. An empty header gets a unique placeholder: several
/// anonymous records must not merge into one protein.
fn accession(header: &str, record: usize) -> String {
    match header.split_ascii_whitespace().next() {
        Some(word) => word.to_string(),
        None => {
            log::warn!("FASTA record #{} has an empty header", record + 1);
            format!("{}{}", PLACEHOLDER, record + 1)
        }
    }
}

const PLACEHOLDER: &str = "unnamed_protein_";

/// Upper-case residues (lower case is used for soft masking) and drop a terminal stop
/// codon (`*`, e.g. translated/Ensembl FASTAs). Anything else that is not a residue
/// makes the peptides containing it invalid; they are dropped when digested.
fn clean(sequence: String) -> String {
    let mut sequence = sequence.to_ascii_uppercase();
    while sequence.ends_with('*') {
        sequence.pop();
    }
    sequence
}

impl Fasta {
    // Parse a string into a fasta database
    pub fn parse<S: Into<String>>(contents: String, decoy_tag: S, generate_decoys: bool) -> Fasta {
        let decoy_tag = decoy_tag.into();

        let mut targets = Vec::new();
        let mut last_id = "";
        let mut records = 0usize;
        // indices into `targets` of records named by a placeholder (empty header)
        let mut placeholders = Vec::new();
        let mut s = String::new();

        for line in contents.as_str().lines() {
            if line.is_empty() {
                continue;
            }
            let line = line.trim();
            if let Some(id) = line.strip_prefix('>') {
                if !s.is_empty() {
                    let acc: Arc<str> = Arc::from(accession(last_id, records));
                    records += 1;
                    let seq = clean(std::mem::take(&mut s));
                    if !acc.contains(&decoy_tag) || !generate_decoys {
                        if acc.starts_with(PLACEHOLDER) && last_id.trim().is_empty() {
                            placeholders.push(targets.len());
                        }
                        targets.push((acc, seq));
                    }
                }
                last_id = id;
            } else {
                s.push_str(line);
            }
        }

        if !s.is_empty() {
            let acc: Arc<str> = Arc::from(accession(last_id, records));
            if !acc.contains(&decoy_tag) || !generate_decoys {
                if acc.starts_with(PLACEHOLDER) && last_id.trim().is_empty() {
                    placeholders.push(targets.len());
                }
                targets.push((acc, clean(s)));
            }
        }

        // placeholder names must not collide with explicit accessions (a record may be
        // called `unnamed_protein_2`); rename colliding placeholders
        // ponytail: a decoy tag occurring inside the placeholder text (e.g. "unnamed_")
        // is not handled; such a tag is not realistic.
        if !placeholders.is_empty() {
            let explicit = targets
                .iter()
                .enumerate()
                .filter(|(ix, _)| !placeholders.contains(ix))
                .map(|(_, (acc, _))| acc.clone())
                .collect::<std::collections::HashSet<_>>();
            for &ix in &placeholders {
                let mut name = targets[ix].0.to_string();
                while explicit.contains(name.as_str()) {
                    name.push('_');
                }
                targets[ix].0 = Arc::from(name);
            }
        }

        Fasta {
            targets,
            decoy_tag,
            generate_decoys,
        }
    }

    pub fn digest(&self, enzyme: &EnzymeParameters) -> Vec<Digest> {
        self.targets
            .par_iter()
            .flat_map_iter(|(protein, sequence)| {
                enzyme
                    .digest(sequence, protein.clone())
                    .into_iter()
                    .filter_map(|mut digest| {
                        if protein.contains(&self.decoy_tag) {
                            if !self.generate_decoys {
                                digest.decoy = true;
                                Some(digest)
                            } else {
                                None
                            }
                        } else {
                            Some(digest)
                        }
                    })
            })
            .collect()
    }

    pub fn iter_chunks(&self, chunk_size: usize) -> impl Iterator<Item = Self> + '_ {
        self.targets
            .chunks(chunk_size)
            .map(move |target_chunk| Self {
                targets: target_chunk.to_vec(),
                decoy_tag: self.decoy_tag.clone(),
                generate_decoys: self.generate_decoys,
            })
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn sequences_are_cleaned() {
        let fasta = Fasta::parse(
            ">sp|P1|A desc\npeptidek\nAAK*\n>\nMKR\n".into(),
            "rev_",
            false,
        );
        assert_eq!(fasta.targets[0].0.as_ref(), "sp|P1|A");
        assert_eq!(fasta.targets[0].1, "PEPTIDEKAAK");
        // an empty header does not panic and gets a unique name
        assert_eq!(fasta.targets[1].0.as_ref(), "unnamed_protein_2");

        // ... also when an explicit accession already uses that name
        let clash = Fasta::parse(
            ">unnamed_protein_2\nPEPTIDEK\n>\nAAAAAAAK\n".into(),
            "rev_",
            false,
        );
        assert_eq!(clash.targets.len(), 2);
        assert_ne!(clash.targets[0].0, clash.targets[1].0);
        assert_eq!(fasta.targets[1].1, "MKR");
    }
}
