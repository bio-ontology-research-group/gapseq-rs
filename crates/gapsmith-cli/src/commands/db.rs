//! `gapsmith db <subcommand>` — introspection of the reference-data tables.
//!
//! Milestone M2 provides `inspect`, which loads every table under the
//! resolved `--data-dir` and prints row counts plus a few summary stats. This
//! is the smoke test that the loaders work end-to-end against the real
//! `dat/` bundle.

use clap::{Parser, Subcommand};
use gapsmith_db::DataRoot;
use gapsmith_io::resolve_data_dir;
use gapsmith_core::SeedStatus;
use std::path::Path;

#[derive(Debug, Parser)]
pub struct Args {
    #[command(subcommand)]
    pub cmd: DbCmd,
}

#[derive(Debug, Subcommand)]
pub enum DbCmd {
    /// Load every table and print row counts.
    Inspect,
}

pub fn run(args: Args, data_dir_override: Option<&Path>) -> anyhow::Result<()> {
    match args.cmd {
        DbCmd::Inspect => inspect(data_dir_override),
    }
}

fn inspect(data_dir_override: Option<&Path>) -> anyhow::Result<()> {
    let root = resolve_data_dir(data_dir_override)?;
    eprintln!("loading tables from `{}` ...", root.display());
    let d = DataRoot::load(&root)?;

    println!("data root: {}", d.root.display());

    let approved = d.seed_rxns.iter().filter(|r| r.gapseq_status == SeedStatus::Approved).count();
    let corrected = d.seed_rxns.iter().filter(|r| r.gapseq_status == SeedStatus::Corrected).count();
    let not_assessed = d.seed_rxns.iter().filter(|r| r.gapseq_status == SeedStatus::NotAssessed).count();
    let removed = d.seed_rxns.iter().filter(|r| r.gapseq_status == SeedStatus::Removed).count();

    println!();
    println!("SEED reactions:      {:>7}  (approved {approved}, corrected {corrected}, not_assessed {not_assessed}, removed {removed})", d.seed_rxns.len());
    println!("SEED metabolites:    {:>7}", d.seed_cpds.len());
    println!("mnxref_seed:         {:>7}", d.mnxref_seed.len());
    println!("mnxref_seed-other:   {:>7}", d.mnxref_seed_other.len());
    println!("meta_pwy:            {:>7}", d.meta_pwy.len());
    println!("kegg_pwy:            {:>7}", d.kegg_pwy.len());
    println!("seed_pwy:            {:>7}", d.seed_pwy.len());
    println!("custom_pwy:          {:>7}", d.custom_pwy.len());
    println!("medium_rules:        {:>7}", d.medium_rules.len());
    println!("complex_subunit:     {:>7}", d.complex_subunit.len());
    println!("subex:               {:>7}", d.subex.len());
    println!("tcdb_substrates:     {:>7}", d.tcdb_substrates.len());
    println!("exception:           {:>7}", d.exception.len());

    println!();
    print_biomass("Gram+", d.biomass_gram_pos.as_ref());
    print_biomass("Gram−", d.biomass_gram_neg.as_ref());
    print_biomass("Archaea", d.biomass_archaea.as_ref());

    // Sanity: the first approved reaction should parse its stoichiometry.
    if let Some(r) = d.seed_rxns.iter().find(|r| r.gapseq_status == SeedStatus::Approved) {
        match r.parse_stoich() {
            Ok(terms) => {
                println!();
                println!(
                    "sample approved reaction: {} ({} terms, reversibility {:?})",
                    r.id,
                    terms.len(),
                    r.reversibility()
                );
            }
            Err(e) => {
                eprintln!("warning: failed to parse stoichiometry for {}: {}", r.id, e);
            }
        }
    }
    Ok(())
}

fn print_biomass(label: &str, t: Option<&gapsmith_db::BiomassTemplate>) {
    match t {
        Some(bm) => {
            let n_components: usize = bm.met_groups.iter().map(|g| g.components.len()).sum();
            println!(
                "biomass {:<8} id={:<10} domain={:<8} {} groups / {} components",
                label,
                bm.id,
                bm.domain,
                bm.met_groups.len(),
                n_components
            );
        }
        None => println!("biomass {:<8} <missing>", label),
    }
}
