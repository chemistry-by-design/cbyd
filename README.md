# Chemistry by Design (CbyD)

This repository contains the CbyD dataset. The code for the paper's analysis and for building the route viewer is in [cbyd_analysis](https://github.com/chemistry-by-design/cbyd_analysis).

## Route files

`data/routes/<Compound>_<Year>_<Author>.json`:

```jsonc
{
  "background": {
    "name": ..., "author": ..., "year": ...,
    "drug_name": ..., "manufacturer": ...,       
    "citation": ..., "doi": ...
  },
  "synthesis": [                                 // one entry per step
    {
      "reaction_id": ...,                        // <Compound>_<Year>_<Author>_<step>
      "reaction_smiles": ...,                    // reactants (incl. reagents) >> product
      "reaction_smiles_mapped": ...,             // reactants > reagents > product
      "mapping_confidence": ...,
      "conditions": {
        "raw_description": ...,
        "GENERAL_REAGENTS": [...], "GENERAL_REAGENTS_SMILES": [...],
        "CATALYSTS": [...],        "CATALYSTS_SMILES": [...],
        "LIGANDS": [...],          "LIGANDS_SMILES": [...],
        "SOLVENTS": [...],         "SOLVENTS_SMILES": [...],
        "MODIFIERS": [...],
        "TEMPERATURE": [...],
        "YIELD": [...]
      }
    }
  ]
}
```

## Route trees (for CASP)

Route files list steps without linking them. For retrosynthesis work, we recommend building trees with `build_route_tree.py`. It is tailored to CbyD: it uses the atom mapping to separate reactants from reagents, strips spectator counter-ions and salts, and links each step to the one that made its reactants. 

The output has one AiZynthFinder-style tree per line, matching the trees in the viewer.

```bash
conda env create -f environment.yml && conda activate cbyd
python build_route_tree.py   # writes data/route_tree.jsonl
```

## Reporting errors

In the viewer, click **Report a mistake** (toolbar for a whole route, or the bottom of a reaction/molecule panel). This lets you email us, open a pre-filled GitHub issue, or edit the route's JSON and submit a pull request. 

Without the viewer, email [cbc-cbyd@arizona.edu](mailto:cbc-cbyd@arizona.edu) with the route name, step, and the problem.

## Submitting a route

1. Draw the route in ChemDraw following the [formatting guide](submit_new_route/General%20Formatting.pdf). Name the file `product_name_YEAR_LastAuthor.cdx`.
2. Check it against the [submit_new_route](submit_new_route) (`GUIDE_EXAMPLE_nimbolide_2025_Sarlah.cdx` is annotated).
3. Email the `.cdx` file and the DOI to [cbc-cbyd@arizona.edu](mailto:cbc-cbyd@arizona.edu).

## Citation

Williams, R. E.; Tran, S. B. A.; Leatherwood, H. M.; Roh, J.; Wolff, A.; Federice, J. G.; Lauta, N. R.; Yu, K.; Nick, S. J.; Jončev, Z.; Schwaller, P.; Coley, C. W.; Njardarson, J. T. Chemistry by Design (CbyD): A Machine-Readable Dataset of Expert-Curated Total Synthesis Routes. *arXiv* **2026**, [arXiv:2611.24996](https://arxiv.org/abs/2611.24996).


```bibtex
@article{williams2026cbyd,
  title   = {Chemistry by Design ({CbyD}): A Machine-Readable Dataset of Expert-Curated Total Synthesis Routes},
  author  = {Williams, Ryan E. and Tran, Suong B. A. and Leatherwood, Hayden M. and Roh, Jihye and Wolff, Aquila and Federice, John G. and Lauta, Nicholas R. and Yu, Kevin and Nick, Steven J. and Jon{\v{c}}ev, Zlatko and Schwaller, Philippe and Coley, Connor W. and Njardarson, Jon T.},
  journal = {arXiv preprint arXiv:2611.24996},
  year    = {2026},
  doi     = {10.48550/arXiv.2611.24996},
  url     = {https://arxiv.org/abs/2611.24996}
}
```
