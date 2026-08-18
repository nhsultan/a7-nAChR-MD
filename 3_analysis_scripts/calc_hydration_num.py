#!/usr/bin/env python3
"""
calc_hydration_num.py

Description:
This script analyzes the coordination environment (hydration number) of ions 
(Na+, K+, Cl-) as they pass through the protein ion channel.

Key Features:
1. Geometric Filtering: Identifies ions within a cylindrical lumen:
   - Radius: 25 Angstroms
   - Z-range: -40 to 100 Angstroms
2. Coordination Shell Analysis: For every ion in the lumen, it counts:
   - Water Oxygens within 3.5 Angstroms
   - Protein Oxygens within 3.5 Angstroms
3. Input/Output: Seamlessly integrates into the bash pipeline via argparse.

Outputs:
A CSV file containing the frame, ion identity, Z-position, and coordination counts.
"""

import argparse
import numpy as np
import pandas as pd
import mdtraj as md
from timeit import default_timer as timer

def main():
    start_time = timer()

    # 1. Setup Argument Parsing
    parser = argparse.ArgumentParser(description="Analyze ion coordination in a channel lumen")
    parser.add_argument("--traj", required=True, help="Path to the trajectory file (.nc, .xtc)")
    parser.add_argument("--top", required=True, help="Path to the topology file (.psf, .pdb)")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files (e.g., results/hydration)")
    args = parser.parse_args()

    out_file = f"{args.out_prefix}_coordination.csv"

    # 2. Configuration & Unit Conversion
    LUMEN_RADIUS_NM = 25.0 / 10.0       # 25 A -> 2.5 nm
    Z_MIN_NM = -40.0 / 10.0             # -40 A -> -4.0 nm
    Z_MAX_NM = 100.0 / 10.0             # 100 A -> 10.0 nm
    COORDINATION_CUTOFF_NM = 3.5 / 10.0 # 3.5 A -> 0.35 nm

    print(f"Loading trajectory: {args.traj}")
    print(f"Loading topology: {args.top}")
    
    try:
        traj = md.load(args.traj, top=args.top)
    except Exception as e:
        print(f"Error loading trajectory: {e}")
        return

    print(f"Loaded {traj.n_frames} frames.")
    topology = traj.topology

    # 3. Atom Selections
    target_resnames = ['SOD', 'NA', 'POT', 'K', 'CLA', 'CL']
    ion_indices = []

    print("Identifying ions and oxygen atoms...")
    for residue in topology.residues:
        if residue.name in target_resnames:
            # Standard single-site ions: take the first (and usually only) atom
            ion_indices.append(next(residue.atoms).index)

    all_ion_indices = np.array(ion_indices)

    # Select Oxygens
    # Water oxygens (HOH, TIP3, WAT) and Protein oxygens
    water_oxygen_indices = topology.select("water and element O")
    protein_oxygen_indices = topology.select("protein and element O")

    if len(all_ion_indices) == 0:
        print("Error: No ions found. Check your selection syntax or trajectory.")
        return

    print(f"Tracking {len(all_ion_indices)} ions, {len(water_oxygen_indices)} water oxygens, and {len(protein_oxygen_indices)} protein oxygens.")

    results = []

    # 4. Frame Iteration and Coordination Calculation
    print("Analyzing frames...")
    for frame_idx, frame in enumerate(traj):
        
        # Extract coordinates for this frame [n_atoms, 3]
        xyz = frame.xyz[0] 
        ion_coords = xyz[all_ion_indices]
        
        # Calculate radial distance and Z-values
        r_sq = ion_coords[:, 0]**2 + ion_coords[:, 1]**2
        z_vals = ion_coords[:, 2]
        
        # Create boolean mask for ions inside the cylinder
        in_lumen_mask = (r_sq < LUMEN_RADIUS_NM**2) & (z_vals > Z_MIN_NM) & (z_vals < Z_MAX_NM)
        
        # Get specific indices of ions inside the lumen for this frame
        lumen_ion_indices = all_ion_indices[in_lumen_mask]

        if len(lumen_ion_indices) == 0:
            continue

        # Iterate over each ion inside the lumen to calculate its specific coordination
        for ion_idx in lumen_ion_indices:
            atom = topology.atom(ion_idx)
            
            # A. Water Coordination (Oxygens within 3.5A)
            water_neighbors = md.compute_neighbors(
                frame, 
                COORDINATION_CUTOFF_NM, 
                query_indices=[ion_idx], 
                haystack_indices=water_oxygen_indices
            )
            num_water_o = len(water_neighbors[0])

            # B. Protein Coordination (Oxygens within 3.5A)
            protein_neighbors = md.compute_neighbors(
                frame, 
                COORDINATION_CUTOFF_NM, 
                query_indices=[ion_idx], 
                haystack_indices=protein_oxygen_indices
            )
            num_protein_o = len(protein_neighbors[0])

            # Z position in Angstroms
            z_pos_angstrom = xyz[ion_idx, 2] * 10.0

            results.append({
                'Frame': frame_idx,
                'Ion_Name': atom.name,
                'Residue_Name': atom.residue.name,
                'Residue_ID': atom.residue.resSeq, 
                'Z_Position_A': z_pos_angstrom,
                'Water_O_Count': num_water_o,
                'Protein_O_Count': num_protein_o,
                'Total_Coordination': num_water_o + num_protein_o
            })

    # 5. Export Data
    df = pd.DataFrame(results)
    
    if not df.empty:
        df.to_csv(out_file, index=False)
        print(f"Analysis complete. Found {len(df)} ion occurrences in the lumen.")
        print(f"Saved results to {out_file}")
    else:
        print("Analysis complete, but no ions were found inside the defined lumen geometry.")

    end_time = timer()
    print(f"Hydration analysis finished in {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()