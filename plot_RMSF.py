#!/usr/bin/env python3

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import sys
import csv
from collections import defaultdict

# -------- USER INPUT --------
TRAJ_FILE = "production.nc"
TOP_FILE  = "solvated.pdb"
ATOM_SELECTION = "protein and backbone"

CSV_FILE = "rmsf.csv"
PNG_FILE = "rmsf.png"
# ----------------------------

print("Loading trajectory...")
traj = md.load(TRAJ_FILE, top=TOP_FILE)

print("Selecting atoms...")
atom_indices = traj.topology.select(ATOM_SELECTION)

if len(atom_indices) == 0:
    print("ERROR: atom selection returned 0 atoms!")
    sys.exit(1)

print(f"Selected {len(atom_indices)} atoms")

print("Using first frame as reference")
ref = traj[0]

print("Aligning trajectory...")
traj.superpose(ref, atom_indices=atom_indices)

print("Computing RMSF...")
rmsf_nm = md.rmsf(traj, ref, atom_indices=atom_indices)
rmsf_angstrom = rmsf_nm * 10

# Residue numbers
top = traj.topology
# Group RMSF by residue
res_rmsf = defaultdict(list)
for atom_idx, r in zip(atom_indices, rmsf_angstrom):
    resid = traj.topology.atom(atom_idx).residue.resSeq
    res_rmsf[resid].append(r)

# Compute average RMSF per residue
resids = sorted(res_rmsf.keys())
rmsf_per_residue = [np.mean(res_rmsf[r]) for r in resids]

# -------- SAVE CSV --------
print("Saving CSV:", CSV_FILE)
np.savetxt(
    CSV_FILE,
    np.column_stack((resids, rmsf_per_residue)),
    header="Residue,RMSF_Angstrom",
    delimiter=",",
    fmt="%d,%.4f"
)

# -------- PLOT --------
print("Saving PNG:", PNG_FILE)
plt.figure()
plt.plot(resids, rmsf_per_residue)
plt.xlabel("Residue")
plt.ylabel("RMSF (Å)")
plt.title("RMSF per residue")
plt.tight_layout()
plt.savefig(PNG_FILE, dpi=300)
plt.show()

print("Done ✔")
