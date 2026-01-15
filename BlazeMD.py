# Blaze MD (ver 0.1 alpha): run dynamics simulations on portable hardware with high speed
# Originally created and developed by Gleb Novikov, forked by PyroTronix
# To use on Mac or Linux terminals:
# python BlazeMD.py your_structure.pdb
# (c) PyroTronix

import os
import sys
import glob 
import time
from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import PDBFixer

# --- STEP 0: CLEAN OLD DATA ---
def clean_old_files(target_pdb, activate: bool = False):
    if not activate:
        return
    print(f"🧹Cleaning Old Data:")
    patterns = ["*.dcd", "*.nc", "*cif", target_pdb]
    for pattern in patterns:
        for file in glob.glob(pattern):
            try:
                os.remove(file)
                print(f"Removed: {file}")
            except Exception as e:
                print(f"Error removing {file}: {e}")

# --- STEP 1: FIX & CLEAN TOPOLOGY ---
def fix_topology(input_pdb, env_pH, ligands_to_remove):
    print(f"--- Step 1: Fix & Clean Topology ({input_pdb}) ---")
    fixer = PDBFixer(filename=input_pdb)

    # === NEW BLOCK: CONVERT MSE TO MET ===
    # This must run BEFORE we process atoms or missing residues
    mse_count = 0
    for residue in fixer.topology.residues():
        if residue.name == 'MSE':
            residue.name = 'MET'
            mse_count += 1
            for atom in residue.atoms():
                if atom.element.symbol == 'Se':
                    # Convert Selenium (Se) to Sulfur (S)
                    atom.element = Element.getBySymbol('S')
                    atom.name = 'SD' # Standard name for Sulfur in MET
    
    if mse_count > 0:
        print(f"🧬 Converted {mse_count} 'MSE' residues to 'MET'.")
    # =====================================

    # 1. Handle Gaps
    fixer.findMissingResidues()
    
    # 2. Blacklist Logic (Remove unwanted ligands)
    modeller = Modeller(fixer.topology, fixer.positions)
    blacklist = set(ligands_to_remove)
    
    residues_to_strip = [r for r in modeller.topology.residues() if r.name in blacklist]

    # Optional: Log kept residues
    unique_residues_kept = {r.name for r in modeller.topology.residues() if r.name not in blacklist}
    print(f"ℹ️  Residues being preserved: {sorted(list(unique_residues_kept))}")

    if residues_to_strip:
        print(f"🧹 Cleaning: Removing {len(residues_to_strip)} unwanted residues/ligands.")
        removed_names = {r.name for r in residues_to_strip}
        print(f"   (Removed types: {removed_names})")
        modeller.delete(residues_to_strip)
    else:
        print("✨ No blacklisted residues found. Keeping all original atoms.")

    # Update fixer
    fixer.topology = modeller.topology
    fixer.positions = modeller.positions

    # 3. Add Hydrogens / Missing Atoms
    print(f"⚛️ Filling missing atoms and hydrogens at pH '{env_pH}'...")
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(env_pH)
    return fixer

# --- STEP 2: SOLVATION ---
def solvate_system(fixer, ff_protein, ff_water, ff_map, box_type, box_padding, ionic_strength, pre_sim_pdb):
    print(f"-- Step 2: Solvation 🫧 {box_type} box with {box_padding} nm --")
    protein_xml = ff_map[ff_protein]

    # Load Forcefields
    try:
        water_xml = f"{ff_protein}/{ff_water}.xml"
        forcefield = ForceField(protein_xml, water_xml)
    except ValueError:
        try:
            water_xml = f"{ff_water}.xml"
            forcefield = ForceField(protein_xml, water_xml)
        except ValueError:
            print(f"⚠️ Could not find {water_xml}. Please ensure the water model XML is available.")

    print(f"🥂 Successfully loaded Forcefield for: {ff_protein} & {ff_water}")

    modeller = Modeller(fixer.topology, fixer.positions)
    
    # 1. Remove water 
    modeller.deleteWater()

    # Check for virtual sites (Required for 4 point water models)
    print(f"✨ Checking for virtual sites...")
    modeller.addExtraParticles(forcefield)
    
    # check if water model is 4-point
    solvent_model = ff_water
    if 'opc' in ff_water.lower() or 'tip4p' in ff_water.lower():
        print(f"💧 Detected 4-point water ({ff_water}). Forcing 'tip4pew' geometry generator.")
        solvent_model = 'tip4pew'

    print(f"🌊 Adding solvent ({ff_water})...")
    modeller.addSolvent(forcefield, 
                        model=solvent_model,
                        padding=box_padding * nanometer, 
                        boxShape=box_type, 
                        ionicStrength=ionic_strength * molar)

    print(f"Final System Composition: {modeller.topology.getNumAtoms()} atoms.")
    
    with open(pre_sim_pdb, 'w') as f:
        PDBFile.writeFile(modeller.topology, modeller.positions, f)
    
    return modeller, forcefield

# --- STEP 3: SETUP & MINIMIZATION & EQUILIBRATION ---
def setup_and_minimize(modeller, forcefield, temperature, pressure, timestep, em_tolerance, nvt_steps, npt_steps, steps_per_stage, pre_production_steps, use_hmr, H_mass):
    print(f"--- Step 3: Setup System, Restraints & Equilibration ---")
    
    # --- HMR LOGIC ---
    system_args = {
        'nonbondedMethod': PME,
        'nonbondedCutoff': 1.2 * nanometer,
        'constraints': HBonds,
    }
    
    if use_hmr:
        print(f"🚀 HMR ENABLED: Repartitioning hydrogen mass to {H_mass} amu.")
        print(f"   -> Timestep is {timestep} ps.")
        system_args['hydrogenMass'] = H_mass * amu 
        
        # Safety Check
        if timestep < 0.003:
            print("⚠️  WARNING: HMR is on, but timestep is small (< 3fs). You can increase it to 4fs for speed!")
    else:
        print(f"🐢 HMR DISABLED: Using standard hydrogen mass.")
        print(f"   -> Timestep is {timestep} ps.")
        
        # Safety Check
        if timestep > 0.0025:
            print("⛔ WARNING: HMR is OFF, but timestep is large (> 2.5fs). The simulation will likely CRASH.")

    # Create System using unpacked arguments
    system = forcefield.createSystem(modeller.topology, **system_args, useDispersionCorrection=True)
    
    # --- Add Barostat ---
    barostat = MonteCarloBarostat(pressure*bar, temperature*kelvin, 25)
    system.addForce(barostat)

    # --- Add Restraints ---
    restraint = CustomExternalForce("k_restraint * periodicdistance(x, y, z, x0, y0, z0)^2")
    restraint.addGlobalParameter("k_restraint", 1000.0 * kilojoules_per_mole / nanometer**2) 
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    
    # Restrain CA, C, N
    for atom in modeller.topology.atoms():
        if atom.name in ('CA', 'C', 'N'):
            restraint.addParticle(atom.index, modeller.positions[atom.index])
    system.addForce(restraint)

    # --- Setup Integrator ---
    integrator = LangevinMiddleIntegrator(temperature*kelvin, 1/picosecond, timestep*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # 1. Minimize
    print("✨ Minimizing energy...")
    simulation.minimizeEnergy(tolerance=em_tolerance*kilojoules_per_mole/nanometer)

    # 2. NVT (Heating)
    barostat.setFrequency(0) # Off for NVT
    simulation.context.reinitialize(preserveState=True)
    simulation.context.setVelocitiesToTemperature(temperature*kelvin)
    simulation.context.setParameter("k_restraint", 1000.0)
    
    print(f"🔥 Phase 1: NVT Heating ({nvt_steps * timestep / 1000:.2f} ns)...")
    simulation.step(nvt_steps)

    # 3. NPT (Pressurizing)
    barostat.setFrequency(25) # On for NPT
    simulation.context.reinitialize(preserveState=True)
    
    print(f"⚖️ Phase 2: NPT Pressurization ({npt_steps * timestep / 1000:.2f} ns)...")
    simulation.step(npt_steps)

    # 4. Gradual Restraint Release
    print(f"🔓 Phase 3: Gradual Restraint Release...")
    k_values = [1000.0, 500.0, 250.0, 100.0, 50.0, 10.0, 0.0]

    for k in k_values:
        simulation.context.setParameter("k_restraint", k)
        simulation.step(steps_per_stage)

    # 5. Unrestrained Equilibration (The "Phase 4" fix)
    print(f"🌊 Phase 4: Unrestrained Equilibration ({pre_production_steps * timestep / 1000:.2f} ns)...")
    simulation.step(pre_production_steps)
    
    # Reset step count for production
    simulation.context.setTime(0)

    # Save equilibrated structure
    with open('equilibrated.cif', 'w') as f:
        PDBxFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)

    return simulation, system

# --- STEP 4: PRODUCTION RUN ---
def run_production(simulation, system, total_steps, timestep, report_interval, output_nc, log_freq, temperature, pressure, use_hmr):
    
    # Determine HMR status
    hmr_status = "ACTIVE" if use_hmr else "DISABLED"
    
    print(f"--- Step 4: Production Run (NPT) 👑 ---")
    print(f"   -> HMR: {hmr_status}")
    print(f"   -> Timestep: {timestep} ps")
    print(f"   -> Temperature: {temperature} K")
    print(f"   -> Pressure: {pressure} ATM")

    # Detect GPU
    platform = simulation.context.getPlatform()
    print(f"✨ Running on Platform: {platform.getName()}")
    if platform.getName() in ['CUDA', 'OpenCL', 'Metal']:
        device_name = platform.getPropertyValue(simulation.context, 'DeviceName')
        print(f"💠 Active GPU: {device_name}")
    else:
        print("💻 Running on CPU (no GPU detected or selected)")

    # Trajectory Reporter (NetCDF or DCD)
    try:
        from mdtraj.reporters import NetCDFReporter
        simulation.reporters.append(NetCDFReporter(output_nc, report_interval))
        print(f"📝 Logging to NetCDF: {output_nc}")
    except ImportError:
        simulation.reporters.append(DCDReporter('production.dcd', report_interval))
        print(f"📝 Logging to DCD: production.dcd")

    # Custom Output Logic
    total_ns = total_steps * timestep / 1000
    steps_per_percent = int(total_steps / 100)
    
    print(f"\n🔮 Simulation Duration: {total_ns:.2f} ns ({total_steps} steps)")
    print(f"\n{'%':>4} {'Step':>10} {'PE (kJ/mol)':>15} {'Speed (ns/day)':>15}")
    print("-" * 55)

    production_start = time.time()
    last_time = time.time()

    for i in range(1, 101):
        simulation.step(steps_per_percent)
        
        if i % log_freq == 0:
            state = simulation.context.getState(getEnergy=True)
            pot_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
            
            current_time = time.time()
            elapsed = current_time - last_time
            # Calculate ns/day based on the ACTUAL timestep
            speed = (steps_per_percent * log_freq * timestep * 0.001) / (elapsed / 86400)
            
            print(f"{i:>3}% {simulation.currentStep:>10} {pot_energy:>15.1f} {speed:>15.1f}")
            last_time = current_time

    # Final Summary
    production_end = time.time()
    total_seconds = production_end - production_start
    avg_speed = (total_steps * timestep * 0.001) / (total_seconds / 86400)

    print("-" * 55)
    print(f"WORK COMPLETED!")
    print(f"⏳ Wall-clock time: {int(total_seconds // 60)}m {int(total_seconds % 60)}s")
    print(f"⚜️ Average performance: {avg_speed:.2f} ns/day")
    print(f"🌀 Trajectory saved.")

# --- MAIN CONTROLLER ---
def main():
    if len(sys.argv) < 2:
        print("Usage: python RoyalMD_beta.py your_structure.pdb")
        sys.exit(1)

    input_pdb = sys.argv[1] 
    
    #                   --- MASTER SETTINGS ---

    # Simulation config
    TARGET_MD_LENGTH =      1               # Simulation length in ns
    USE_HMR =               True            # Use Hydrogen Mass Repartitioning for longer timesteps
    
    # File config
    pre_sim_pdb =           'solvated.pdb'  # Pre-simulation prepared structure
    output_nc =             'production.nc' # Output trajectory file

    # Equilibration config (ultra fast and light)
    em_tolerance =          2.5             # kJ/mol/nm
    nvt_time =              0.1             # ns
    npt_time =              0.1             # ns
    pre_production_time =   0.1             # ns
    report_interval =       1               # ps - Time interval per frame
    restrain_release_time = 20              # ps per restraint release stage
    
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
                             ]              # Blacklist of ligands to remove
    # Forcefield config
    ff_protein =            'amber19'       # Protein forcefield
    ff_water =              'tip3p'         # Water model, use tip4pew for 4-point water models (best)
    ff_map = {
        'amber19': 'amber19-all.xml',
        'amber14': 'amber14-all.xml',
        'amber99sb': 'amber99sb.xml',
        'amber99sbildn': 'amber99sbildn.xml',
        'amber03': 'amber03.xml',
        'charmm36': 'charmm36.xml'
    }

    # --- BASE LOGIC ---

    # Calculate Timestep & Steps based on HMR
    timestep = 0
    if USE_HMR:
        timestep = time_step_HMR
        print(f"🚀 HMR MODE: Active. Target Time: {TARGET_MD_LENGTH} ns")
    else:
        timestep = time_step_noHMR
        print(f"🐢 NORMAL MODE: Active. Target Time: {TARGET_MD_LENGTH} ns")

    # Time to steps conversion
    production_steps = int((TARGET_MD_LENGTH * 1000) / timestep)
    report_interval = int(report_interval / timestep)
    nvt_steps = int((nvt_time * 1000) / timestep)
    npt_steps = int((npt_time * 1000) / timestep)
    pre_production_steps = int((pre_production_time * 1000) / timestep)
    restrain_release_steps = int(restrain_release_time / timestep)

    # --- PIPELINE ---

    # Step 0
    clean_old_files(pre_sim_pdb, activate=True) 
    
    # Step 1
    fixer = fix_topology(input_pdb, env_pH, ligands_to_remove) 

    # Step 2
    modeller, forcefield = solvate_system(
        fixer, 
        ff_protein, 
        ff_water, 
        ff_map, 
        box_type, 
        box_padding, 
        ionic_strength, 
        pre_sim_pdb,
        )
    
    # Step 3
    simulation, system = setup_and_minimize(
        modeller, 
        forcefield, 
        temperature, 
        pressure, 
        timestep, 
        em_tolerance, 
        nvt_steps, 
        npt_steps,
        restrain_release_steps,
        pre_production_steps,
        USE_HMR,
        hydrogen_mass,
        )
    
    # Step 4
    run_production(
        simulation, 
        system, 
        production_steps, 
        timestep, 
        report_interval, 
        output_nc, 
        log_freq, 
        temperature, 
        pressure,
        USE_HMR,
        )
    
if __name__ == "__main__":
    main()