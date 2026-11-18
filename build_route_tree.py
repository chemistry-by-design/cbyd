"""Write data/route_tree.jsonl: one nested retrosynthetic tree per route in
data/routes (format in the README). A route with several or no targets is
written as ``{"id": ..., "trees": [...]}`` and reported.

    conda run -n cbyd python build_route_tree.py
"""

import json
import os
import re
from collections import Counter
from functools import lru_cache
from multiprocessing import Pool

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

# Neutralizer for protonation-invariant matching (an acid and its conjugate base —
# HCl/[Cl-], TFA/TFA- — map to the same neutral form). Optional: falls back to
# exact matching if the standardizer isn't available.
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
    _UNCHARGER = rdMolStandardize.Uncharger()
except Exception:  # pragma: no cover - depends on RDKit build
    _UNCHARGER = None


# SMILES utilities----------------------------------------------
def split_rxn_smiles(rxn_smiles):
    """Split a reaction SMILES on '>' that is not part of '->' or '<-'."""
    return re.split(r"(?<!-)>(?!-)", rxn_smiles)


def canonicalize_smiles(smiles, remove_map=True):
    """Canonicalize a SMILES string, optionally stripping atom-map numbers.
    Raises ValueError on an unparseable SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    if remove_map:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
    try:
        smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        mol = Chem.MolFromSmiles(smiles)
    except Exception:
        pass
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def safe_canonicalize(smiles, remove_map=True):
    """canonicalize_smiles that returns None instead of raising."""
    try:
        return canonicalize_smiles(smiles, remove_map=remove_map)
    except Exception:
        return None


def _protonation_key(smiles):
    """Canonical SMILES of the neutralized (uncharged) molecule, so a species and
    its protonation/charge variant share a key (e.g. 'Cl' and '[Cl-]' -> 'Cl';
    TFA and TFA- -> 'O=C(O)C(F)(F)F'). None if unparseable."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    if _UNCHARGER is not None:
        try:
            mol = _UNCHARGER.uncharge(mol)
        except Exception:
            pass
    try:
        return Chem.MolToSmiles(mol)
    except Exception:
        return None


# Reaction cleaning: which molecules of a reaction become tree nodes---------------------------------------------

# Spectator salt-forming acids and their conjugate anions: counter-ions that
# should not become tree nodes.
_HALOGEN_ACID_SMILES_RAW = [
    "Cl", "Cl-", "[Cl-]", "Cl[H]", "ClH", "Cl.H", "HCl", "[H][Cl]",
    "Br", "Br-", "[Br-]", "Br[H]", "BrH", "Br.H", "HBr", "[H][Br]",
    "I", "I-", "[I-]", "IH", "I.H", "HI", "[H][I]",
    "F", "F-", "[F-]", "FH", "F.H", "HF", "[H][F]",
]

# Common non-halide counter-ions (the acid and every anionic protonation state) —
# treated like halide counter-ions.
_OXOACID_SALT_SMILES_RAW = [
    "OP(=O)(O)O", "OP(=O)(O)[O-]", "OP(=O)([O-])[O-]", "[O-]P(=O)([O-])[O-]",  # H3PO4 / H2PO4- / HPO4-2 / PO4-3
    "OS(=O)(=O)O", "OS(=O)(=O)[O-]", "[O-]S(=O)(=O)[O-]",                       # H2SO4 / HSO4- / SO4-2
    "Cc1ccc(S(=O)(=O)O)cc1", "Cc1ccc(S(=O)(=O)[O-])cc1",                        # TsOH / TsO- (tosylate)
    "OS(C)(=O)=O", "[O-]S(C)(=O)=O",                                            # MsOH / MsO- (mesylate)
    "OS(=O)(=O)C(F)(F)F", "[O-]S(=O)(=O)C(F)(F)F",                              # TfOH / TfO- (triflate)
    "OC(=O)C(F)(F)F", "[O-]C(=O)C(F)(F)F",                                      # TFA / TFA- (trifluoroacetate)
]


@lru_cache(maxsize=None)
def get_halogen_acid_set():
    """Canonical SMILES of spectator counter-ions (halide acids/anions plus the
    phosphate, sulfate, sulfonate and TFA acids and their anions)."""
    acids = set()
    for salt in _HALOGEN_ACID_SMILES_RAW + _OXOACID_SALT_SMILES_RAW:
        can = safe_canonicalize(salt)
        if can:
            acids.add(can)
    return frozenset(acids)


def _comp_map_nums(comp):
    """Atom-map numbers used in a single SMILES component (empty if unparseable)."""
    mol = Chem.MolFromSmiles(comp)
    if mol is None:
        return set()
    return {a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum()}


def remove_spectator_halides(rxn_smiles):
    """Mapping-aware halide cleanup.

    A halide (HBr/HCl/[Cl-]/…) is KEPT only if it is *incorporated into an organic
    product* — i.e. one of its atom-map numbers appears in a product component that
    is itself not a halide salt (e.g. HBr -> [Br:1][CH:2]… keeps HBr). A halide that
    maps only to another halide/counter-ion (HCl -> [Cl-]) or isn't mapped in is a
    spectator and is stripped (when a >1-carbon organic remains on the side)."""
    halogen_acids = get_halogen_acid_set()

    parts = split_rxn_smiles(rxn_smiles)
    if len(parts) < 3:
        return rxn_smiles
    reactants_smi, agents, products_smi = parts[0], parts[1], parts[-1]

    # atom-map numbers that belong to a NON-halide (organic) product component
    organic_map_nums = set()
    for c in products_smi.split("."):
        if not c:
            continue
        can = safe_canonicalize(c)
        if can is not None and can in halogen_acids:
            continue                       # a standalone halide product doesn't count
        organic_map_nums |= _comp_map_nums(c)

    def count_carbons(smi):
        mol = Chem.MolFromSmiles(smi)
        return 0 if mol is None else sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)

    def strip_side(side_smi):
        comps = [c for c in side_smi.split(".") if c]
        non_halide, halide = [], []
        for c in comps:
            can = safe_canonicalize(c)
            is_halide = can is not None and can in halogen_acids
            incorporated = bool(_comp_map_nums(c) & organic_map_nums)
            if is_halide and not incorporated:
                halide.append(c)           # spectator halide -> candidate to strip
            else:
                non_halide.append(c)       # organic, or a halide fixed into the product
        if halide and any(count_carbons(c) > 1 for c in non_halide):
            return non_halide
        return comps

    reactants_new = strip_side(reactants_smi)
    products_new = strip_side(products_smi)
    return ".".join(reactants_new) + ">" + agents + ">" + ".".join(products_new)


def check_if_reactant_product_same(rxn_smiles):
    """True only if reactants == products (a genuine salt/complex no-op like
    ``A.HCl >> A.HCl``), ignoring spectator counter-ions when both sides carry them
    alongside other species. Comparison is charge-SENSITIVE, so a real protonation
    step (``A >> AH+``) is NOT dropped."""
    salts = get_halogen_acid_set()

    def split_salt(smi):
        comps = [c for c in smi.split(".") if c]
        canon = [(c, safe_canonicalize(c)) for c in comps]
        return ([c for c, ca in canon if ca not in salts],
                [c for c, ca in canon if ca in salts])

    parts = split_rxn_smiles(rxn_smiles)
    if len(parts) < 3:
        return False
    reactant_smiles, product_smiles = parts[0], parts[2]

    nh_react, h_react = split_salt(reactant_smiles)
    nh_prod, h_prod = split_salt(product_smiles)

    if h_react and h_prod and nh_react and nh_prod:
        rc, pc = ".".join(nh_react), ".".join(nh_prod)
    else:
        rc, pc = reactant_smiles, product_smiles
    a, b = safe_canonicalize(rc), safe_canonicalize(pc)
    return a is not None and a == b


def clean_rxn_duplicates(rxn_smiles):
    """Remove molecules unchanged between both sides (e.g. A.B>>C.B -> A>>C)."""
    parts = split_rxn_smiles(rxn_smiles)
    if len(parts) < 3:
        return rxn_smiles
    reactants_smi, agents, products_smi = parts[0], parts[1], parts[-1]
    reactants = reactants_smi.split(".")
    products = products_smi.split(".")

    reactants_canon = [safe_canonicalize(r) for r in reactants]
    products_canon = [safe_canonicalize(p) for p in products]
    common = set(reactants_canon) & set(products_canon)
    common.discard(None)

    reactants_new = [
        r for r, rc in zip(reactants, reactants_canon)
        if rc not in common or reactants_canon.count(rc) > products_canon.count(rc)
    ]
    products_new = [
        p for p, pc in zip(products, products_canon)
        if pc not in common or products_canon.count(pc) > reactants_canon.count(pc)
    ]

    # Also drop a spectator that appears on BOTH sides differing only by protonation
    # / charge (e.g. HCl / [Cl-], TFA / TFA-). Guarded: only when both sides keep >1
    # component, so a genuine sole reactant/product is never stripped.
    if len(reactants_new) > 1 and len(products_new) > 1:
        r_keys = [_protonation_key(r) for r in reactants_new]
        p_keys = [_protonation_key(p) for p in products_new]
        shared = {k for k in r_keys if k} & {k for k in p_keys if k}
        if shared:
            r2 = [r for r, k in zip(reactants_new, r_keys) if k not in shared]
            p2 = [p for p, k in zip(products_new, p_keys) if k not in shared]
            if r2 and p2:                 # never empty a side
                reactants_new, products_new = r2, p2

    return ".".join(reactants_new) + ">" + agents + ">" + ".".join(products_new)


def clean_reaction(rxn_smiles):
    """(reactants, products) as canonical SMILES of the molecules that become tree
    nodes, or None if the reaction is skipped (salt/complex formation) or invalid."""
    if check_if_reactant_product_same(rxn_smiles):
        return None
    cleaned = clean_rxn_duplicates(rxn_smiles)
    # mapping-aware: keep a halide that's incorporated into an organic product
    # (e.g. HBr adding across a diene), strip spectator counter-ions (HCl -> [Cl-]).
    cleaned = remove_spectator_halides(cleaned)
    parts = split_rxn_smiles(cleaned)
    if len(parts) < 3:
        return None
    reactants = [safe_canonicalize(r) for r in parts[0].split(".") if r]
    products = [safe_canonicalize(p) for p in parts[-1].split(".") if p]
    reactants = [r for r in reactants if r]
    products = [p for p in products if p]
    if not reactants or not products:
        return None
    return reactants, products


# Route -> tree -------------------------------------------------------
NODE_CAP = 2500           # guard against pathologically convergent blow-up


def mapped_smiles(rxn):
    """The atom-mapped reaction SMILES of one route reaction."""
    return (rxn.get("reaction_smiles_mapped") or "").strip()


def build_trees(synthesis):
    """Nested retrosynthetic tree(s) for a route's ``synthesis`` list.

    Returns ``(trees, warnings, connected, n_targets, orphan_reactions)``.

    The tree is built from the atom-mapped SMILES, so reagents and solvents (its
    middle section) are not molecule nodes. Starting at each final product, it walks
    backwards (product <- reaction <- reactants), making a fresh copy of every node:
    a building block used by several reactions gets its whole upstream sub-tree
    duplicated, so each path is an independent chain. A reactant repeated within one
    reaction (a dimerization A.A>>B, or 4x MeI) appears that many times as a child.
    """
    warnings = []
    agents = set()      # canonical SMILES that ever appear as a reagent/agent
    per_rxn = []        # index -> (reactants, products) | None
    source = []         # index -> the reaction SMILES the nodes come from
    produced, consumed = set(), set()
    produced_by = {}    # canonical SMILES -> indices of the reactions that make it

    for index, rxn in enumerate(synthesis):
        smiles = mapped_smiles(rxn) or (rxn.get("reaction_smiles") or "").strip()
        source.append(smiles)
        cleaned = None
        if smiles:
            sparts = split_rxn_smiles(smiles)
            if len(sparts) >= 3:
                agents.update(a for a in map(safe_canonicalize, sparts[1].split(".")) if a)
            try:
                cleaned = clean_reaction(smiles)
            except Exception as e:  # pragma: no cover - defensive
                warnings.append(f"Rxn {index + 1}: failed to parse ({e}).")
        per_rxn.append(cleaned)
        if cleaned is not None:
            reactants, products = cleaned
            consumed.update(reactants)
            produced.update(products)
            for p in products:
                produced_by.setdefault(p, []).append(index)

    reactants_of = [cl[0] if cl else [] for cl in per_rxn]
    used_rxn = set()
    n_nodes = [0]
    capped = [False]

    def expand(smiles, path, is_root):
        n_nodes[0] += 1
        node = {"target": smiles} if is_root else {"smiles": smiles}
        node["type"] = "mol"
        prods = [j for j in produced_by.get(smiles, []) if j not in path]
        if not prods or n_nodes[0] > NODE_CAP:
            capped[0] = capped[0] or n_nodes[0] > NODE_CAP
            return node
        producer = max(prods)                         # the reaction that made it
        n_nodes[0] += 1
        used_rxn.add(producer)
        children = []
        for r, mult in Counter(reactants_of[producer]).items():
            child = expand(r, path | {producer}, False)
            children.extend([child] * mult)
        node["children"] = [{
            "type": "reaction",
            "reaction_id": str(synthesis[producer].get("reaction_id", producer)),
            "reaction_smiles_mapped": mapped_smiles(synthesis[producer]),
            "children": children,
        }]
        return node

    trees = []
    for root in sorted(p for p in produced if p not in consumed):
        tree = expand(root, set(), True)
        if root not in agents:                        # a reagent-derived root is a byproduct
            trees.append(tree)
    if capped[0]:
        warnings.append(f"Tree truncated at {NODE_CAP} nodes (highly convergent route).")

    # reactions never reached above (salt/complex, unparseable, or orphaned)
    orphan_reactions = []
    for index in range(len(synthesis)):
        if index in used_rxn:
            continue
        number = index + 1
        if per_rxn[index] is None:
            if source[index] and check_if_reactant_product_same(source[index]):
                warnings.append(f"Rxn {number}: reactants == products (salt/complex); "
                                "shown with no molecule edges.")
            elif not source[index]:
                warnings.append(f"Rxn {number}: no reaction SMILES available.")
            else:
                warnings.append(f"Rxn {number}: SMILES could not be parsed — shown as invalid (red).")
        else:
            orphan_reactions.append(number)   # valid reaction, product/reactant mismatch

    n_targets = len(trees)
    connected = n_targets == 1 and not orphan_reactions
    return trees, warnings, connected, n_targets, orphan_reactions


def route_tree(route_path):
    """(tree, warnings) for one route JSON. A single-target route gives one tree
    dict; otherwise a list of trees."""
    with open(route_path) as f:
        synthesis = json.load(f)["synthesis"]
    trees, warnings, connected, n_targets, orphans = build_trees(synthesis)
    if not connected:
        warnings.append(f"route not a single connected tree: {n_targets} target(s), "
                        f"orphans={orphans}")
    return (trees[0] if len(trees) == 1 else trees), warnings


def _one(route_path):
    name = os.path.splitext(os.path.basename(route_path))[0]
    try:
        tree, warnings = route_tree(route_path)
    except Exception as e:                 # malformed route file
        return name, None, [f"failed: {e}"]
    return name, tree, warnings


def build_route_trees(routes_dir, out_path, processes=None):
    files = sorted(f for f in os.listdir(routes_dir) if f.endswith(".json"))
    n_ok = n_multi = n_failed = 0
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with Pool(processes or min(32, os.cpu_count() or 1)) as pool, open(out_path, "w") as out:
        for name, tree, warnings in pool.imap(_one, [os.path.join(routes_dir, f) for f in files],
                                              chunksize=8):
            for w in warnings:
                print(f"[WARNING] {name}: {w}")
            if tree is None:
                n_failed += 1
                continue
            if isinstance(tree, list):
                n_multi += 1
                record = {"id": name, "trees": tree}
            else:
                n_ok += 1
                record = {"id": name, **tree}
            out.write(json.dumps(record) + "\n")
    print(f"[SUCCESS]: wrote {out_path} — {n_ok + n_multi} routes "
          f"({n_multi} not a single tree, {n_failed} failed to load)")
    return out_path


if __name__ == "__main__":
    repo = os.path.dirname(os.path.abspath(__file__))
    build_route_trees(os.path.join(repo, "data", "routes"),
                      os.path.join(repo, "data", "route_tree.jsonl"))
