#!/usr/bin/env python3
"""Validate an SBML L3V1+FBC2+groups file emitted by gapsmith-sbml.

Run under the uv venv at gapsmith/tools/.sbml-validate/.

Two checks:

1. libSBML consistency check — reports every error/warning the canonical
   validator finds. Exit 2 on structural errors, 1 on warnings-only.
2. COBRApy round-trip — loads the document, prints a summary, runs FBA on
   the model's objective (if any), and verifies reaction/metabolite counts
   match a user-supplied pair when given.

Usage:
    validate_sbml.py <file.xml>            # validation + cobra load
    validate_sbml.py <file.xml> N_RXN N_MET  # additionally check counts
"""

from __future__ import annotations

import sys
from pathlib import Path


def check_libsbml(path: Path) -> tuple[int, int]:
    """Run libSBML consistency checks. Returns (errors, warnings)."""
    import libsbml

    reader = libsbml.SBMLReader()
    doc = reader.readSBML(str(path))

    print(f"== libSBML ==  level={doc.getLevel()}  version={doc.getVersion()}")
    print(f"  plugins: fbc={doc.getPlugin('fbc') is not None}  "
          f"groups={doc.getPlugin('groups') is not None}")

    # Strict L3V1 + FBC + groups check.
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MODELING_PRACTICE, False)
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, True)
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_GENERAL_CONSISTENCY, True)
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_IDENTIFIER_CONSISTENCY, True)
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_MATHML_CONSISTENCY, True)
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_SBO_CONSISTENCY, True)
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_OVERDETERMINED_MODEL, True)
    n = doc.checkConsistency()
    errs = 0
    warns = 0
    for i in range(n):
        err = doc.getError(i)
        sev = err.getSeverity()
        label = {
            libsbml.LIBSBML_SEV_INFO: "info",
            libsbml.LIBSBML_SEV_WARNING: "warn",
            libsbml.LIBSBML_SEV_ERROR: "ERROR",
            libsbml.LIBSBML_SEV_FATAL: "FATAL",
        }.get(sev, "?")
        if sev >= libsbml.LIBSBML_SEV_ERROR:
            errs += 1
        elif sev == libsbml.LIBSBML_SEV_WARNING:
            warns += 1
        line = err.getLine()
        short = err.getShortMessage().replace("\n", " ")
        msg = err.getMessage().replace("\n", " ").strip()
        print(f"  [{label}] L{line} {short}: {msg[:180]}")
    if errs == 0 and warns == 0:
        print("  libSBML: 0 errors, 0 warnings")
    else:
        print(f"  libSBML: {errs} errors, {warns} warnings")

    model = doc.getModel()
    if model is not None:
        fbc = model.getPlugin("fbc")
        groups = model.getPlugin("groups")
        n_obj = fbc.getNumObjectives() if fbc else 0
        n_gp = fbc.getNumGeneProducts() if fbc else 0
        n_groups = groups.getNumGroups() if groups else 0
        print(
            f"  model: {model.getNumSpecies()} species, "
            f"{model.getNumReactions()} reactions, "
            f"{model.getNumCompartments()} compartments, "
            f"{model.getNumParameters()} parameters, "
            f"fbc:objectives={n_obj}, fbc:geneProducts={n_gp}, "
            f"groups:groups={n_groups}"
        )

    return errs, warns


def check_cobra(path: Path, expect_rxns: int | None, expect_mets: int | None) -> bool:
    """Load with COBRApy and report the summary."""
    try:
        import cobra
    except Exception as e:
        print(f"== COBRApy ==  import failed: {e}")
        return False

    try:
        m = cobra.io.read_sbml_model(str(path))
    except Exception as e:
        print(f"== COBRApy ==  load failed: {e.__class__.__name__}: {e}")
        return False

    print(f"== COBRApy == loaded model id={m.id!r}")
    print(
        f"  {len(m.metabolites)} metabolites, {len(m.reactions)} reactions, "
        f"{len(m.genes)} genes, {len(m.groups)} groups"
    )
    obj_rxns = [r.id for r in m.reactions if r.objective_coefficient != 0]
    print(f"  objective reactions: {obj_rxns}")

    ok = True
    if expect_rxns is not None and len(m.reactions) != expect_rxns:
        print(f"  FAIL: expected {expect_rxns} reactions, got {len(m.reactions)}")
        ok = False
    if expect_mets is not None and len(m.metabolites) != expect_mets:
        print(f"  FAIL: expected {expect_mets} metabolites, got {len(m.metabolites)}")
        ok = False

    # Try FBA if there's an objective with any stoichiometry.
    if obj_rxns and any(r.metabolites for r in m.reactions):
        try:
            sol = m.optimize()
            print(f"  FBA optimum: {sol.objective_value:.6g}  status={sol.status}")
        except Exception as e:
            print(f"  FBA failed: {e}")
    return ok


def main() -> int:
    argv = sys.argv[1:]
    if not argv:
        print(__doc__)
        return 64
    path = Path(argv[0])
    expect_rxns = int(argv[1]) if len(argv) >= 2 else None
    expect_mets = int(argv[2]) if len(argv) >= 3 else None

    errs, warns = check_libsbml(path)
    print()
    cobra_ok = check_cobra(path, expect_rxns, expect_mets)

    if errs > 0 or not cobra_ok:
        return 2
    if warns > 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
