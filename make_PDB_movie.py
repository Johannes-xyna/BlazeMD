#!/usr/bin/env python3

import mdtraj as md
import sys

# -------- USER INPUT --------
TRAJ_FILE = "production.nc"
TOP_FILE  = "solvated.pdb"
OUTPUT_PDB = "movie_protein.pdb"

# Selection: usually "protein" is best for movies to keep file size down.
# Use "all" if you really want to see the water (warning: huge file).
ATOM_SELECTION = "protein" 

# Stride: Save every Nth frame. 
# 1 = save every frame (huge file). 10 = save every 10th frame.
STRIDE = 10 
# ----------------------------

print(f"Loading trajectory: {TRAJ_FILE}")
print(f"Using topology: {TOP_FILE}")

# Load the trajectory
traj = md.load(TRAJ_FILE, top=TOP_FILE)
print(f"Original frames: {traj.n_frames}")

# 1. Apply Stride (Reduce number of frames)
if STRIDE > 1:
    print(f"Applying stride (taking every {STRIDE}th frame)...")
    traj = traj[::STRIDE]
    print(f"Frames to process: {traj.n_frames}")

# 2. Image Molecules (Fix Periodic Boundary Conditions)
# This ensures the protein doesn't "jump" across the edges of the box
print("Centering molecules and fixing periodic boundaries...")
traj.image_molecules()

# 3. Select Atoms (Strip water/ions if desired)
print(f"Selecting atoms: '{ATOM_SELECTION}'...")
atom_indices = traj.topology.select(ATOM_SELECTION)

if len(atom_indices) == 0:
    print("ERROR: atom selection returned 0 atoms!")
    sys.exit(1)

# Create a new trajectory object with only the selected atoms
traj_subset = traj.atom_slice(atom_indices)

# 4. Superpose (Align)
# This aligns all frames to the first frame so the protein doesn't rotate/drift
print("Aligning trajectory to frame 0...")
traj_subset.superpose(traj_subset[0])

# 5. Save
print(f"Saving PDB movie to: {OUTPUT_PDB}")
traj_subset.save_pdb(OUTPUT_PDB)

print("Done ✔")