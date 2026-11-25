"""Build one reference QAD mol2/sdf per ligand copy (HETA, HETB).

Bond orders and atom order come from complex_prep/qad/QAD.mol2.
Coordinates come from the solvated CHARMM-GUI system (step3_pbcsetup.crd),
so both references sit in the exact frame of the complex we simulate.

Run:  module load anaconda/2024.10 ; python3 make_refs.py
"""
from pathlib import Path
import os
import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D

PREP = Path(os.environ.get("PREP_DIR", "prep"))  # CHARMM-GUI solvated build
OUT = Path(__file__).resolve().parent / "refs"
OUT.mkdir(exist_ok=True)

template = Chem.MolFromMol2File(str(PREP / "qad/QAD.mol2"), removeHs=False)
assert template is not None and template.GetNumAtoms() == 34

# rtf atom names, in rtf order (the crd uses the same order)
rtf_names = [l.split()[1] for l in open(PREP / "qad/qad.rtf") if l.startswith("ATOM")]

def crd_segment(segid):
    rows = []
    for line in open(PREP / "step3_pbcsetup.crd").readlines()[4:]:
        t = line.split()
        if len(t) >= 8 and t[7] == segid:
            rows.append((t[3], float(t[4]), float(t[5]), float(t[6])))
    return rows

for seg in ("HETA", "HETB"):
    rows = crd_segment(seg)
    assert len(rows) == 34, (seg, len(rows))
    assert [r[0] for r in rows] == rtf_names, "crd atom order differs from qad.rtf"
    m = Chem.Mol(template)
    conf = m.GetConformer()
    for i, (name, x, y, z) in enumerate(rows):
        el = m.GetAtomWithIdx(i).GetSymbol()
        # rtf names start with the element symbol (S, O1, C10, H3 ...)
        assert name.rstrip("0123456789") == el, (i, name, el)
        conf.SetAtomPosition(i, Point3D(x, y, z))
        m.GetAtomWithIdx(i).SetProp("_TriposAtomName", name)
    # sanity: every bond 0.9-1.9 A in the complex frame
    pos = conf.GetPositions()
    lens = [np.linalg.norm(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()]) for b in m.GetBonds()]
    assert 0.9 < min(lens) and max(lens) < 1.9, (seg, min(lens), max(lens))
    Chem.AssignStereochemistryFrom3D(m)
    m.SetProp("_Name", f"qad_{seg[-1].lower()}")
    with Chem.SDWriter(str(OUT / f"qad_{seg[-1].lower()}.sdf")) as w:
        w.write(m)
    print(f"{seg}: 34 atoms, bonds {min(lens):.3f}-{max(lens):.3f} A, "
          f"centroid {pos.mean(axis=0).round(2)} -> refs/qad_{seg[-1].lower()}.sdf")

a = Chem.MolFromMolFile(str(OUT / "qad_a.sdf"), removeHs=False).GetConformer().GetPositions()
b = Chem.MolFromMolFile(str(OUT / "qad_b.sdf"), removeHs=False).GetConformer().GetPositions()
print(f"closest heavy-atom contact between the two copies: "
      f"{min(np.linalg.norm(a[i]-b[j]) for i in range(20) for j in range(20)):.2f} A")
