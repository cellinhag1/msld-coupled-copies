# msld-coupled-copies

Multisite λ-dynamics for a **homodimer that binds two identical copies of the
same ligand**, where both copies must change substituent together.

## The problem

A symmetric receptor with two equivalent sites binds two copies of one
molecule. In a standard MSλD setup each site gets its own λ, so copy A can
carry substituent 3 while copy B carries substituent 5. That state is not a
real molecule being tested: no such compound exists, and the free energy
collected there does not correspond to anything you could synthesise.

With 7 substituents and two independent sites there are 49 combinations, of
which only the 7 matched ones are physical. Most of the sampling is spent in
states that have to be discarded.

## The fix

Put substituent *N* of copy A and substituent *N* of copy B into **one
selection**, so they become **one BLOCK and one λ**:

```
substituent N  ->  { copy A atoms of sub N }  OR  { copy B atoms of sub N }
                   ------------- one block, one lambda -------------
```

Both copies therefore always carry the same substituent. The 49-state space
collapses to the 7 physical ones. `adex` excludes the cross terms
(A-sub*N* with B-sub*M*, *N* ≠ *M*), while A-sub*N* and B-sub*N* interact
inside block *N* and so are scaled by the same λ, which is exact at the
endpoints.

To ALF and BLOCK this looks like an ordinary one-site, seven-substituent
problem. The coupling lives entirely in how the selections are built, so
nothing downstream needs to know about it.

The selection construction is in `build/build_complex.py`; look for the
`COUPLED selections` comment.

## Layout

| File | Purpose |
|---|---|
| `prep/make_refs.py` | Build per-copy reference ligands from the crystal pose, taking bond orders from a mol2 and coordinates from the solvated build |
| `prep/prep_ligands.py` | Align the analog panel onto each copy, write mol2 files, run CGenFF, then MCS and charge renormalization. Fixed random seed and constrained embedding, so the core reproduces to ~0.05 Å |
| `prep/write_copy_pdbs.py` | Emit core and substituent fragment PDBs per copy, and verify the atom-name map by reproducing the build coordinates |
| `build/build_complex.py` | pyCHARMM system build: delete the crystal ligands, generate both ligand copies with their patches, build the coupled selections, set up BLOCK |
| `build/alf_info.py` | Site counts, temperature and bias settings for ALF |

## Free energy with one ligand per copy

The complex leg has two ligand copies and the solvent leg has one, so the
legs are not symmetric:

    ddG_bind(pair) = ddG_complex - 2 * ddG_solvent

The free ligands are independent in solution, hence the factor of two, and the
solvent error propagates doubled as well.

## Stereochemistry

Substituents that introduce a new stereocentre must be placed on the same face
in both copies, or the two copies are no longer the same molecule and the
coupling is not physically meaningful. `prep_ligands.py` enforces one
consistent face by constrained embedding against the reference, and the chosen
isomer is recorded per substituent. **This is a modelling choice.** If the
experimental stereochemistry is known, set it explicitly rather than accepting
the default.

## Configuration

Nothing is hard-coded to a machine. Set before running:

```bash
export PREP_DIR=/path/to/charmm-gui/output   # solvated build
export PANEL_CSV=/path/to/panel.csv          # "name,SMILES" per line
export MSLD_PY_PREP=/path/to/msld-py-prep    # Vilseck-Lab clone
export MSLDPREP_PATH=/path/to/msldprep       # parent dir of the msldprep package
export TOPPAR=/path/to/toppar                # CHARMM36 / CGenFF
```

No ligand panel is distributed here; supply your own CSV.

## Requirements

pyCHARMM with BLaDE, [ALF](https://github.com/RyanLeeHayes/ALF), CGenFF, RDKit,
and [Vilseck-Lab/msld-py-prep](https://github.com/Vilseck-Lab/msld-py-prep).

## License

MIT, see `LICENSE`.
