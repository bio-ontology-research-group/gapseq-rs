//! Loaders for gapseq's reference-data tables in `dat/`.
//!
//! Every loader is self-contained: pass in a path (or the `DataRoot` for the
//! batch loader), get back a plain Rust value. No implicit paths, no globals.
//!
//! The parsers are tolerant of the exact minor variations seen in real
//! gapseq `dat/` files (some tables have 13 columns, others 14; MNXref
//! prefixes blocks of `#`-comments). Row validation errors carry line
//! numbers so the culprit can be found quickly.

pub mod biomass;
pub mod common;
pub mod complex;
pub mod exception;
pub mod medium_rules;
pub mod mnxref;
pub mod pathway;
pub mod seed;
pub mod stoich_parse;
pub mod subex;
pub mod tcdb;

pub use biomass::{BiomassTemplate, BiomassComponent, BiomassGroup, BiomassError};
pub use common::DbError;
pub use complex::{ComplexSubunitEntry, ComplexSubunitTable};
pub use exception::ExceptionRow;
pub use medium_rules::MediumRule;
pub use mnxref::{MnxrefSeed, MnxrefSeedOther};
pub use pathway::{PathwayRow, PathwayTable, PwySource};
pub use seed::{SeedCpdRow, SeedRxnRow, load_seed_metabolites, load_seed_reactions};
pub use stoich_parse::{parse_stoichiometry, StoichTerm, StoichParseError};
pub use subex::SubexRow;
pub use tcdb::TcdbSubstrateRow;

use std::path::{Path, PathBuf};

/// The full set of reference tables gapseq consults. Fields are populated
/// lazily as the relevant loaders are called; [`DataRoot::load`] loads
/// everything at once.
#[derive(Default)]
pub struct DataRoot {
    pub root: PathBuf,
    pub seed_rxns: Vec<SeedRxnRow>,
    pub seed_cpds: Vec<SeedCpdRow>,
    pub mnxref_seed: Vec<MnxrefSeed>,
    pub mnxref_seed_other: Vec<MnxrefSeedOther>,
    pub meta_pwy: PathwayTable,
    pub kegg_pwy: PathwayTable,
    pub seed_pwy: PathwayTable,
    pub custom_pwy: PathwayTable,
    pub medium_rules: Vec<MediumRule>,
    pub complex_subunit: ComplexSubunitTable,
    pub subex: Vec<SubexRow>,
    pub tcdb_substrates: Vec<TcdbSubstrateRow>,
    pub exception: Vec<ExceptionRow>,
    pub biomass_gram_pos: Option<BiomassTemplate>,
    pub biomass_gram_neg: Option<BiomassTemplate>,
    pub biomass_archaea: Option<BiomassTemplate>,
}

impl DataRoot {
    /// Load every reference table from a single `dat/` root. Individual
    /// missing files (e.g. if the user only has the core bundle) are logged
    /// at `warn!` and left empty.
    pub fn load(root: impl AsRef<Path>) -> Result<Self, DbError> {
        let root = root.as_ref().to_path_buf();
        let mut out = DataRoot { root: root.clone(), ..Default::default() };

        out.seed_rxns = load_seed_reactions(root.join("seed_reactions_corrected.tsv"))?;
        out.seed_cpds = load_seed_metabolites(root.join("seed_metabolites_edited.tsv"))?;
        out.mnxref_seed = mnxref::load_mnxref_seed(root.join("mnxref_seed.tsv"))?;
        out.mnxref_seed_other = mnxref::load_mnxref_seed_other(root.join("mnxref_seed-other.tsv"))?;
        out.meta_pwy = PathwayTable::load(root.join("meta_pwy.tbl"), PwySource::MetaCyc)?;
        out.kegg_pwy = PathwayTable::load(root.join("kegg_pwy.tbl"), PwySource::Kegg)?;
        out.seed_pwy = PathwayTable::load(root.join("seed_pwy.tbl"), PwySource::Seed)?;
        out.custom_pwy = PathwayTable::load(root.join("custom_pwy.tbl"), PwySource::Custom)?;
        out.medium_rules = medium_rules::load(root.join("medium_prediction_rules.tsv"))?;
        out.complex_subunit = ComplexSubunitTable::load(root.join("complex_subunit_dict.tsv"))?;
        out.subex = subex::load(root.join("subex.tbl"))?;
        out.tcdb_substrates = tcdb::load_substrates(root.join("tcdb_substrates.tbl"))?;
        out.exception = exception::load(root.join("exception.tbl"))?;

        let bm_dir = root.join("biomass");
        out.biomass_gram_pos = BiomassTemplate::load_opt(bm_dir.join("biomass_Gram_pos.json"))?;
        out.biomass_gram_neg = BiomassTemplate::load_opt(bm_dir.join("biomass_Gram_neg.json"))?;
        out.biomass_archaea = BiomassTemplate::load_opt(bm_dir.join("biomass_archaea.json"))?;

        Ok(out)
    }
}
