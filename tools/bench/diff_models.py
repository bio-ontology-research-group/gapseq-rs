"""Semantic comparison of two gapsmith-filled SBML models."""
import sys
import cobra

def summarize(path):
    m = cobra.io.read_sbml_model(path)
    rxn_ids = {r.id for r in m.reactions}
    met_ids = {me.id for me in m.metabolites}
    try:
        m.objective.direction = "max"
        sol = m.optimize()
        bio = sol.objective_value
        status = sol.status
    except Exception as e:
        bio = None
        status = f"err: {e}"
    return {
        "rxns": len(rxn_ids),
        "mets": len(met_ids),
        "rxn_ids": rxn_ids,
        "met_ids": met_ids,
        "biomass": bio,
        "fba_status": status,
    }

def compare(a, b):
    ja_rxn = len(a["rxn_ids"] & b["rxn_ids"]) / max(1, len(a["rxn_ids"] | b["rxn_ids"]))
    ja_met = len(a["met_ids"] & b["met_ids"]) / max(1, len(a["met_ids"] | b["met_ids"]))
    only_a = a["rxn_ids"] - b["rxn_ids"]
    only_b = b["rxn_ids"] - a["rxn_ids"]
    return {
        "rxn_jaccard": ja_rxn,
        "met_jaccard": ja_met,
        "rxns_only_A": sorted(only_a),
        "rxns_only_B": sorted(only_b),
        "n_only_A": len(only_a),
        "n_only_B": len(only_b),
    }

a_path, b_path = sys.argv[1], sys.argv[2]
a = summarize(a_path)
b = summarize(b_path)
print(f"A ({a_path}): {a['rxns']} rxns, {a['mets']} mets, biomass={a['biomass']}, status={a['fba_status']}")
print(f"B ({b_path}): {b['rxns']} rxns, {b['mets']} mets, biomass={b['biomass']}, status={b['fba_status']}")
cmp = compare(a, b)
print(f"reaction Jaccard: {cmp['rxn_jaccard']:.4f}")
print(f"metabolite Jaccard: {cmp['met_jaccard']:.4f}")
print(f"reactions only in A: {cmp['n_only_A']}")
print(f"reactions only in B: {cmp['n_only_B']}")
if cmp["n_only_A"] <= 30:
    print(f"  A-only: {cmp['rxns_only_A']}")
if cmp["n_only_B"] <= 30:
    print(f"  B-only: {cmp['rxns_only_B']}")
