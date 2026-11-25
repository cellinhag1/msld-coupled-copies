"""Write per-copy PDBs for the coupled two-copy complex leg.

copy A (segid LIGA): build.sting coordinates (aligned onto HETA).
copy B (segid LIGB): the same build atoms, with coordinates taken from each
  source ligand's copy-B embedding (aligned onto HETB). Atom correspondence:
  build name -> (source ligand, original mol2 name) via translation.txt ->
  RDKit atom index via the names saved in copy_coords.json.

Check: mapping copy A through the same route must reproduce build.sting's own
coordinates (to PDB precision) -- proves the name mapping is right.

Output: complex/prep/corea.pdb, coreb.pdb, site1_sub<N>_fraga.pdb, _fragb.pdb
"""
import json, re
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
B = HERE / "msld_build"
OUT = HERE / "complex/prep"
OUT.mkdir(parents=True, exist_ok=True)

cc = json.loads((B / "copy_coords.json").read_text())
name2idx = {lig: {v: int(k) for k, v in d["names"].items()} for lig, d in cc.items()}

# translation.txt -> {section: (source ligand, {new: old})}
sections, cur = {}, None
reflig = re.search(r"^REFLIG\s+(\S+)", (B / "MCS_for_MSLD.txt").read_text(), re.M).group(1)
for line in (B / "translation.txt").read_text().splitlines():
    if line.strip() == "CORE":
        cur = "core"; sections[cur] = (reflig, {})
    elif line.startswith("SITE"):
        m = re.match(r"SITE \d+ (\S+) \(from (\S+)\)", line)
        cur = m.group(1); sections[cur] = (m.group(2), {})
    elif cur and len(line.split()) == 2:
        old, new = line.split(); sections[cur][1][new] = old


def rewrite(src, key, dst, copy, segid):
    lig, new2old = sections[key]
    out, maxdev = [], 0.0
    for line in src.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        idx = name2idx[lig][new2old[name]]
        xa = np.array(cc[lig]["xyz_a"][idx])
        xb = np.array(cc[lig]["xyz_b"][idx])
        own = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        maxdev = max(maxdev, np.abs(own - xa).max())          # copy-A route must match the build
        x = own if copy == "a" else xb
        new = line[:30] + "%8.3f%8.3f%8.3f" % tuple(x) + line[54:72] + segid
        out.append(new)
    assert maxdev < 2e-3, f"{src.name}: copy-A mapping deviates {maxdev:.4f} A from build -> name map wrong"
    dst.write_text("\n".join(out) + "\nEND\n")
    return len(out), maxdev


for copy, seg in (("a", "LIGA"), ("b", "LIGB")):
    n, dev = rewrite(B / "build.sting/core.pdb", "core", OUT / f"core{copy}.pdb", copy, seg)
    print(f"core{copy}.pdb  {n} atoms  (copy-A check max dev {dev:.4f} A)")
    for f in sorted((B / "build.sting").glob("site1_sub*_frag.pdb")):
        key = f.name.replace("_frag.pdb", "")
        n, dev = rewrite(f, key, OUT / f"{key}_frag{copy}.pdb", copy, seg)
        print(f"{key}_frag{copy}.pdb  {n} atoms  from {sections[key][0]}  (check {dev:.4f} A)")
