# Blaze MD 🔥

**A lightweight, highly optimized OpenMM Molecular Dynamics pipeline for speed and precision.**

Blaze MD is a high-performance simulation kit designed to maximize sampling speed (ns/day) on portable hardware and consumer GPUs. It serves as an all-in-one wrapper for **OpenMM**, automating topology fixing, solvation, equilibration, and production runs with state-of-the-art optimizations like Hydrogen Mass Repartitioning (HMR). Installation time is less than 5 minutes.

> **Forked from:** [RoyalMD](https://github.com/TheVisualHub/RoyalMD) | **Core Engine:** OpenMM

## Key Features

*   **Automated Topology Repair:** Automatically converts MSE to MET, fills missing atoms, and adds hydrogens appropriate for pH 7.0.
*   **Smart Solvation:** Supports 3-point and 4-point water models (e.g., TIP4PEW) with automatic virtual site handling.
*   **High-Speed Optimization:** Uses **Hydrogen Mass Repartitioning (HMR)** to allow for a **4 fs timestep**, significantly increasing simulation speed without sampling loss. Expect 2x simulations speed.
*   **Robust Equilibration:** Performs Energy Minimization $\rightarrow$ NVT $\rightarrow$ NPT $\rightarrow$ Gradual Restraint Release $\rightarrow$ Unrestrained Pre-production.
*   **GPU Accelerated:** Auto-detects CUDA, OpenCL, or Metal to utilize maximum hardware performance.

## Performance on different GPUs
The list below shows the avg. ns/day speed achieved with different GPU's. The test was done using the parameters shown below. Higher speeds can be achieved with less detailed force fields and water models. A four point water model (i.e. opc or tip4pew) reduces the speed by around 23% compared to three point water models (i.e. tip3p).

    GPU Name            |   HMR Enabled  |  HMR Disabled
    RTX 4060 Laptop GPU | **540 ns/day** | **270 ns/day**

Contact me if you attempt this on a GPU that isn't listed here.

---

## Installation

I recommend using **Micromamba** or **Conda** to manage dependencies. Blaze MD relies on `OpenMM`, `PDBFixer`, and `MDTraj`.

### 1. Prerequisite: Install Micromamba or Conda
If you don't have a package manager, install [Miniforge](https://github.com/conda-forge/miniforge) or [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

### 2. Create the Environment
Run the following commands in your Linux/macOS terminal:

```bash
# Create a new environment named 'blazemd'
micromamba create -n blazemd -c conda-forge openmm pdbfixer mdtraj matplotlib numpy netcdf4 python=3.12

# Activate the environment
micromamba activate blazemd
```

*(Note: You can replace `micromamba` with `conda` if you prefer standard Anaconda).*

---

## Main Pipeline: `BlazeMD.py`

This is the core execution engine. It takes a raw PDB file (crystal structure, structure prediction etc.) and performs the entire MD simulation process.

### Usage
```bash
python BlazeMD.py your_structure.pdb
```

### Quick test
```bash
python BlazeMD.py ./test/cyclotide_miniprotein.pdb
```

### How does it work?
1.  **Cleaning:** Removes old simulation artifacts. [PDBFixer]
2.  **Fixing:** Repairs the PDB (standardizes residues, removes blacklisted ligands like SO4/EDO). [PDBFixer]
3.  **Solvation:** Builds a Dodecahedron box with padding and neutralizes ionic strength (0.15M). [OpenMM]
4.  **Equilibration:**
    *   *Minimization* (removes steric clashes). [OpenMM]
    *   *NVT Heating* (brings system to 300K). [OpenMM]
    *   *NPT Pressurizing* (stabilizes density at 1 bar). [OpenMM]
    *   *Restraint Release* (gradually releases protein backbone). [OpenMM]
5.  **Production:** Runs the simulation (default 1000ns) and saves trajectory to `production.nc`. [OpenMM]

### Configuration
To change simulation parameters (length, temperature, forcefields), open `BlazeMD.py` in a text editor (or VS Code) and modify the **MASTER SETTINGS** block near the bottom. Make sure you learn the importance of reliable equilibrium parameters. The ones below are not ideal for publication material. Publication grade parameters are shown in the parentheses to the side.

```python
# Simulation config
TARGET_MD_LENGTH =      1000            # Simulation length in ns
USE_HMR =               True            # Use Hydrogen Mass Repartitioning for longer timesteps

# File config
pre_sim_pdb =           'solvated.pdb'  # Pre-simulation prepared structure
output_nc =             'production.nc' # Output trajectory file

# Equilibration config (a bit coarse)
em_tolerance =          2.5             # kJ/mol/nm
nvt_time =              0.5             # ns 
npt_time =              0.5             # ns                                (2.5 ns)
pre_production_time =   1.5             # ns                                (3   ns)
report_interval =       10              # ps - Time interval per frame
restrain_release_time = 200             # ps per restraint release stage    (250 ns)

# HMR config
hydrogen_mass =         3.0             # amu scaling factor, for HMR
time_step_HMR =         0.004           # ps, 4 fs
time_step_noHMR =       0.002           # ps, 2 fs

# Debug / Logging config
log_freq =              1               # Log every % progress in terminal

# Environment config
temperature =           300             # K
pressure =              1.0             # bar
box_padding =           1.5             # nm
box_type =              'dodecahedron'  # Box shape
ionic_strength =        0.15            # M
env_pH =                7.0             # pH
ligands_to_remove =     ['SO4','EDO', 
                            'LIG', 'LIH', 
                            'lig', 'lih',
                            ]           # Blacklist of ligands to remove
# Forcefield config
ff_protein =            'amber19'       # Protein forcefield
ff_water =              'tip4pew'       # Water model
```

---

## Post-Processing & Analysis

Blaze MD includes three scripts to analyze your data immediately after the simulation finishes. Just make sure all the python scripts are in the working directory and when the MD simulation is done then execute the script as shown below. No parameters are passed to the scripts. File inputs can be changed in the top of the scripts.

### 1. RMSD Analysis (`plot_RMSD.py`)
Calculates the Root Mean Square Deviation of the protein backbone over time to assess stability.

*   **Reads:** `production.nc`, `solvated.pdb`
*   **Output:** `rmsd.csv`, `rmsd.png`

**Usage:**
```bash
python plot_RMSD.py
```

### 2. RMSF Analysis (`plot_RMSF.py`)
Calculates the Root Mean Square Fluctuation per residue. This identifies flexible regions (loops) versus rigid regions (helices/sheets).

*   **Reads:** `production.nc`, `solvated.pdb`
*   **Output:** `rmsf.csv`, `rmsf.png`

**Usage:**
```bash
python plot_RMSF.py
```

### 3. Movie Generation (`make_PDB_movie.py`)
Raw trajectories often have visual artifacts due to Periodic Boundary Conditions (atoms jumping across box edges). This script centers the protein, removes water (for smaller file size), and aligns the frames.

*   **Reads:** `production.nc`, `solvated.pdb`
*   **Output:** `movie_protein.pdb` (Load this into PyMOL or ChimeraX)

**Usage:**
```bash
python make_PDB_movie.py
```

*Note: You can adjust the `STRIDE` variable inside the script to reduce file size (e.g., save every 10th frame).*

---

## Output Files

| File | Description |
| :--- | :--- |
| `solvated.pdb` | The prepared system (protein + water + ions) used as topology. |
| `equilibrated.cif` | Snapshot of the system after equilibration, before production. |
| `production.nc` | The main binary trajectory file (NetCDF format). |
| `production.dcd` | Alternative trajectory format (if NetCDF is unavailable). |
| `rmsd.png` / `rmsf.png` | Visual plots of your simulation dynamics. |
| `movie_protein.pdb` | View the trajectory as a movie. Molstar Viewer is advised. |

---

## Credits & License

**Blaze MD** is a fork of [RoyalMD](https://github.com/TheVisualHub/RoyalMD).
*   **Original Author:** Gleb Novikov
*   **Fork Developed by:** PyroTronix 
*   **Contact:** johannes.friis@xyna.bio

This software is provided "as is" under the MIT License. High-speed dynamics are achieved via OpenMM.