#!/usr/bin/env python3

import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import sys
import csv

# -------- USER INPUT --------
TRAJ_FILE = "production.nc"
TOP_FILE  = "solvated.pdb"
ATOM_SELECTION = "protein and backbone"

CSV_FILE = "rmsd.csv"
PNG_FILE = "rmsd.png"
# ----------------------------

print("Loading trajectory...")
traj = md.load(TRAJ_FILE, top=TOP_FILE)

# Time in ps
time_ps = traj.time
if time_ps is None:
    print("WARNING: No time info found – using frame index")
    time_ps = np.arange(len(traj))

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

print("Computing RMSD...")
rmsd_nm = md.rmsd(traj, ref, atom_indices=atom_indices)
rmsd_angstrom = rmsd_nm * 10

# -------- SAVE CSV --------
print("Saving CSV:", CSV_FILE)
with open(CSV_FILE, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Time_ps", "RMSD_Angstrom"])
    for t, r in zip(time_ps, rmsd_angstrom):
        writer.writerow([t, r])

# -------- PLOT --------
print("Saving PNG:", PNG_FILE)
plt.figure()
plt.plot(time_ps, rmsd_angstrom)
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (Å)")
plt.title("RMSD vs Time")
plt.tight_layout()
plt.savefig(PNG_FILE, dpi=300)
plt.show()

print("Done ✔")
