"""STING ligand panel prep for coupled two-copy MSLD (adapted from Zoe's
journal_club/py_prep_crimm/zoe_msld_prep.ipynb).

Differences from Zoe's notebook:
  * TWO reference poses (QAD copy A = HETA, copy B = HETB, from the solvated
    complex, see make_refs.py). Every analog is embedded onto BOTH, with the
    same atom order, so copy A and copy B share one set of atom names.
  * Fixed random seed for ConstrainedEmbed (reproducible), and each analog's
    stereocentre is taken from its copy-A embedding and enforced on copy B, so
    both copies are the same stereoisomer.
  * msld-py-prep runs on the copy-A mol2s only (one hybrid topology). Copy-B
    coordinates are mapped onto the build afterwards by write_copy_b.py.

Stages:  python3 prep_ligands.py align     (align + mol2 + CGenFF)
         python3 prep_ligands.py mcs       (msld-py-prep pass 1 -> MCS_for_MSLD.txt)
         python3 prep_ligands.py crn       (msld-py-prep pass 2 -> build.sting/)
Run with: module load anaconda/2024.10
"""
import os
import json, re, shutil, subprocess, sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdFMCS import FindMCS
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

HERE = Path(__file__).resolve().parent
sys.path.insert(0, os.environ.get("MSLDPREP_PATH", "."))  # dir containing the msldprep package
from msldprep import sybyl_typer  # noqa: E402

WORK = HERE / "msld_build"
MOL2_DIR = WORK / "mol2s"
CSV = Path(os.environ.get("PANEL_CSV", "panel.csv"))  # one "name,SMILES" per line
SYSNAME = "sting"
REF_ID = "qad"
SEED = 20261019
MSLD_PY_PREP = Path(os.environ["MSLD_PY_PREP"])  # clone of Vilseck-Lab/msld-py-prep
TOPPAR_SRC = Path(os.environ["TOPPAR"])  # CHARMM36/CGenFF toppar directory


def bond_path(atom_ids, mol):
    s = set(atom_ids)
    return [b.GetIdx() for b in mol.GetBonds() if b.GetBeginAtomIdx() in s and b.GetEndAtomIdx() in s]


def embed_on(mol, ref):
    """Same as msldprep.align.align_common_core, but with a fixed seed.
    Returns (core-pinned atom count, core RMSD)."""
    mcs = FindMCS([Chem.RemoveHs(mol), Chem.RemoveHs(ref)], completeRingsOnly=True)
    core = Chem.PathToSubmol(ref, bond_path(ref.GetSubstructMatch(mcs.queryMol), ref), useQuery=True)
    AllChem.ConstrainedEmbed(mol, core, useTethers=True, randomseed=SEED)
    return core.GetNumAtoms(), float(mol.GetProp("EmbedRMS"))


def cip(mol):
    Chem.AssignStereochemistryFrom3D(mol)
    return {a.GetIdx(): a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")}


def stage_align():
    for d in (WORK, MOL2_DIR):
        d.mkdir(exist_ok=True)
    ref_a = Chem.MolFromMolFile(str(HERE / "refs/qad_a.sdf"), removeHs=False)
    ref_b = Chem.MolFromMolFile(str(HERE / "refs/qad_b.sdf"), removeHs=False)
    smiles = [l.strip() for l in open(CSV) if l.strip()]
    ligs = {REF_ID: (Chem.Mol(ref_a), Chem.Mol(ref_b))}
    report = []
    # Face rule: every substituent replaces the SAME prochiral QAD hydrogen on its
    # carbon -- H2 on C2 (acid side), H4 on C3 (ketone side). Both are the roomier
    # positions in the pocket (identical in HETA and HETB), and one face per carbon
    # keeps the remaining chain H's in the common core.
    ref_pos = ref_a.GetConformer().GetPositions()
    target_h = {3: 22, 4: 24}   # QAD atom idx C2 -> H2, C3 -> H4
    ref_heavy = Chem.RemoveHs(ref_a)
    for i, smi in enumerate(smiles):
        lid = f"mol_{i}"
        best = None
        for iso in EnumerateStereoisomers(Chem.MolFromSmiles(smi)):
            cand = AllChem.AddHs(iso)
            nA, rA = embed_on(cand, ref_a)
            X = cand.GetConformer().GetPositions()
            match = cand.GetSubstructMatch(ref_heavy)        # analog idx for each QAD heavy idx
            inv = {a: q for q, a in enumerate(match)}
            core = set(match)
            for a in cand.GetAtoms():
                if a.GetIdx() in core or a.GetSymbol() == "H":
                    continue
                anchor = [inv[n.GetIdx()] for n in a.GetNeighbors() if n.GetIdx() in core]
                if anchor:
                    d = np.linalg.norm(X[a.GetIdx()] - ref_pos[target_h[anchor[0]]])
            if best is None or d < best[0]:
                best = (d, cand, nA, rA)
        dface, mA, nA, rA = best
        assert dface < 1.0, f"{lid}: substituent not on the target face ({dface:.2f} A)"
        stereo = cip(mA)                      # stereocentre picked by the face rule
        mB = Chem.Mol(mA)                     # same atoms, now with chiral tags set
        mB.RemoveAllConformers()
        nB, rB = embed_on(mB, ref_b)
        assert cip(mB) == stereo, f"{lid}: copy B embedded as a different stereoisomer"
        ligs[lid] = (mA, mB)
        report.append((lid, Chem.MolToSmiles(Chem.RemoveHs(mA)), nA, rA, rB, stereo))

    print(f"{'ligand':7s} {'pinned':>6s} {'rmsA':>6s} {'rmsB':>6s}  stereo   SMILES")
    for lid, smi, n, ra, rb, st in report:
        print(f"{lid:7s} {n:6d} {ra:6.3f} {rb:6.3f}  {str(list(st.values())):8s} {smi}")

    names_all = {}
    with Chem.SDWriter(str(WORK / "aligned_a.sdf")) as wa, Chem.SDWriter(str(WORK / "aligned_b.sdf")) as wb:
        for lid, (mA, mB) in ligs.items():
            mA.SetProp("_Name", lid); mB.SetProp("_Name", lid)
            wa.write(mA); wb.write(mB)
            sybyl_typer.annotate_sybyl_types(mA)
            names = sybyl_typer.write_mol2_file(mA, str(MOL2_DIR / f"{lid}.mol2"), mol_name="LIG")
            assert len(set(names.values())) == len(names)
            names_all[lid] = {
                "names": {str(k): v for k, v in names.items()},
                "xyz_a": mA.GetConformer().GetPositions().round(4).tolist(),
                "xyz_b": mB.GetConformer().GetPositions().round(4).tolist(),
            }
    (WORK / "copy_coords.json").write_text(json.dumps(names_all))
    print("wrote aligned_a.sdf, aligned_b.sdf, mol2s/*.mol2, copy_coords.json")

    cmd = " ; ".join(f"cgenff {l}.mol2 --all > {l}.str 2> {l}.cgenff.log" for l in ligs)
    subprocess.run(["bash", "-lc", f"module load cgenff/2026_May && {cmd}"], cwd=MOL2_DIR, check=True)
    for l in ligs:
        txt = (MOL2_DIR / f"{l}.str").read_text()
        m = re.search(r"^RESI\s+\S+\s+(\S+)\s+!\s*param penalty=\s*([\d.]+)\s*;\s*charge penalty=\s*([\d.]+)", txt, re.M)
        print(f"CGenFF {l:6s} charge {m.group(1):>7s}  param penalty {m.group(2):>7s}  charge penalty {m.group(3):>7s}")


def stage_mcs():
    lig_ids = [REF_ID] + [f"mol_{i}" for i in range(6)]
    for f in ["msld_chk.py", "msld_mcs.py", "msld_crn.py", "msld_prm.py", "msld_wrt.py",
              "CGenFF_atomtypes.txt", "Lg_Solvate.sh", "vis_check.py"]:
        shutil.copy(MSLD_PY_PREP / f, WORK / f)
    src = (WORK / "msld_mcs.py").read_text()
    src, n = re.subn(r'refname\s*=\s*"[^"]*"', f'refname = "{REF_ID}"', src)
    assert n == 1
    (WORK / "msld_mcs.py").write_text(src)
    # same msld_crn.py patches as Zoe's notebook (debug print + charge-offset fix)
    src = (WORK / "msld_crn.py").read_text()
    src = re.sub(r"\n\s*print\('site = ', site, 'frag = '.*\)", "", src)
    a = "        sitesum=[]\n        Qfrag.append([])\n"
    b = ("            Qfrag[site].append(pd.Series(np.zeros(len(Fatoms[site][frag])),index=Fatoms[site][frag]))\n"
         "            for mol in range(len(mols)):\n                if mols[mol] == frags[site][frag]:\n                    break\n")
    c = "                QQ=np.rint(sitesum) # round to nearest integer with numpy\n                QQ=[int(q) for q in QQ]\n"
    assert src.count(a) == 1 and src.count(b) == 1 and src.count(c) == 1
    src = src.replace(a, a + "        fragmol=[]\n").replace(b, b + "            fragmol.append(mol)\n")
    src = src.replace(c, "                if nsites == 1:\n"
                         "                    QQ=[int(round(float(rtfinfo[fragmol[f]]['QNET']) - float(rtfinfo[refnum]['QNET']))) for f in range(len(sitesum))]\n"
                         "                else:\n"
                         "                    QQ=np.rint(sitesum)\n                    QQ=[int(q) for q in QQ]\n")
    (WORK / "msld_crn.py").write_text(src)
    tp = WORK / "toppar"; tp.mkdir(exist_ok=True)
    for f in ["top_all36_prot.rtf", "par_all36m_prot.prm", "top_all36_na.rtf", "par_all36_na.prm",
              "top_all36_carb.rtf", "par_all36_carb.prm", "top_all36_cgenff.rtf", "par_all36_cgenff.prm",
              "toppar_water_ions.str"]:
        link = tp / f
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(TOPPAR_SRC / f)
        assert link.exists()
    shutil.copy(MSLD_PY_PREP / "toppar" / "toppar.str", tp / "toppar.str")
    for l in lig_ids:
        for ext in (".mol2", ".str"):
            shutil.copy(MOL2_DIR / f"{l}{ext}", WORK / f"{l}{ext}")
    (WORK / "mol_list.txt").write_text("\n".join(lig_ids) + "\n")
    sys.path.insert(0, str(WORK))
    import os; os.chdir(WORK)
    import msld_chk, msld_mcs
    msld_chk.MsldCHK("mol_list.txt")
    print("reference:", msld_mcs.MsldMCS("mol_list.txt", "MCS_for_MSLD.txt", cutoff=0.8, debug=False))


def stage_crn():
    import os; os.chdir(WORK); sys.path.insert(0, str(WORK))
    import msld_crn, msld_prm, msld_wrt
    out = f"build.{SYSNAME}"
    nsubs = [int(x) for x in re.search(r"^NSUBS\s+(.*)$", Path("MCS_for_MSLD.txt").read_text(), re.M).group(1).split()]
    if Path(out).exists():
        shutil.rmtree(out)
    msld_crn.MsldCRN("MCS_for_MSLD.txt", out, [[] for _ in nsubs], [[] for _ in nsubs], ChkQChange=True, verbose=True, debug=False)
    msld_prm.MsldPRM(out, True, verbose=False, debug=False)
    msld_wrt.writeALF_Files(SYSNAME, out, True)
    print("nsubs", Path(out, "nsubs").read_text().split(), "nblocks", Path(out, "nblocks").read_text().split())


if __name__ == "__main__":
    {"align": stage_align, "mcs": stage_mcs, "crn": stage_crn}[sys.argv[1]]()
