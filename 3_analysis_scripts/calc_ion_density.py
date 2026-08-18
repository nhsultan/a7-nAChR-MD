#!/usr/bin/env python3
"""
calc_ion_density.py

Description:
This script calculates the density of specified ions (Na+, K+, Cl-) along the Z-axis, 
within a cylindrical volume centered on the protein channel.

Key Features:
1. Centers the trajectory frame-by-frame on the protein backbone to correct for X-Y drift.
2. Filters ions based on a cylindrical geometry defined by:
   - Radius: 25 Angstroms
   - Z-axis range: Dynamically determined to cover the entire simulation box.
3. Calculates the number density (average count per frame) of ions along the Z-axis 
   in 0.5 Angstrom bins.
   
Outputs:
A single, well-formatted text file containing the Z-coordinates and the calculated 
average density for each tracked ion type.
"""

import argparse
import numpy as np
import mdtraj as md
from timeit import default_timer as timer

def main():
    start_time = timer()

    # 1. Setup Argument Parsing
    parser = argparse.ArgumentParser(description="Calculate ion density along the Z-axis in a cylinder")
    parser.add_argument("--traj", required=True, help="Path to the trajectory file (.nc, .xtc)")
    parser.add_argument("--top", required=True, help="Path to the topology file (.psf, .pdb)")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files (e.g., results/density)")
    args = parser.parse_args()

    out_file = f"{args.out_prefix}_IonCount.txt"

    # 2. Load the trajectory
    print(f"Loading trajectory: {args.traj}")
    print(f"Loading topology: {args.top}")
    try:
        traj = md.load(args.traj, top=args.top)
    except OSError:
        print("Error: Files not found. Please check your file paths.")
        return

    print(f"Loaded {traj.n_frames} frames.")

    # 3. Define Parameters
    # Convert fixed radius from Angstroms to nanometers (mdtraj native unit)
    radius_nm = 25.0 / 10.0

    # Determine Z-axis range dynamically from the trajectory to cover the entire box
    print("Calculating Z-axis limits from trajectory...")
    global_z_min_nm = np.min(traj.xyz[:, :, 2])
    global_z_max_nm = np.max(traj.xyz[:, :, 2])

    # Convert to Angstroms and round to capture the edges perfectly
    z_min_angstrom = np.floor(global_z_min_nm * 10.0)
    z_max_angstrom = np.ceil(global_z_max_nm * 10.0)

    # Update nm limits for filtering
    z_min_nm = z_min_angstrom / 10.0
    z_max_nm = z_max_angstrom / 10.0

    print(f"Dynamic Z-range determined: {z_min_angstrom} A to {z_max_angstrom} A")

    # 4. Select Atoms
    sod_indices = traj.topology.select('resname SOD or resname NA')
    pot_indices = traj.topology.select('resname POT or resname K')
    cal_indices = traj.topology.select('resname CLA or resname CL')
    protein_indices = traj.topology.select('protein and backbone')
    
    print(f"Found {len(sod_indices)} Na ions, {len(pot_indices)} K ions, {len(cal_indices)} Cl ions.")

    # 5. Data Collection
    sod_z_all = []
    pot_z_all = []
    cal_z_all = []

    print("Processing frames and extracting cylindrical densities...")
    
    # Helper function to process each ion type per frame
    def collect_ions(indices, xyz_frame, center_x, center_y, storage_list):
        if len(indices) == 0:
            return
            
        ion_xyz = xyz_frame[indices]
        
        # Shift coordinates relative to dynamic protein center
        ion_x = ion_xyz[:, 0] - center_x
        ion_y = ion_xyz[:, 1] - center_y
        ion_z = ion_xyz[:, 2] 
        
        # Cylindrical Filter: Check if within radius and within Z bounds
        r_squared = ion_x**2 + ion_y**2
        mask = (r_squared < radius_nm**2) & (ion_z >= z_min_nm) & (ion_z <= z_max_nm)
        
        # Extract Z-values of ions inside the cylinder and convert to Angstroms
        ions_in_cylinder_z = ion_z[mask] * 10.0
        storage_list.extend(ions_in_cylinder_z)

    for i in range(traj.n_frames):
        xyz_frame = traj.xyz[i]
        
        # Calculate dynamic protein center for this frame to correct X-Y drift
        if len(protein_indices) > 0:
            prot_xyz = xyz_frame[protein_indices]
            cx = np.mean(prot_xyz[:, 0])
            cy = np.mean(prot_xyz[:, 1])
        else:
            cx = 0.0
            cy = 0.0

        # Process each ion type
        collect_ions(sod_indices, xyz_frame, cx, cy, sod_z_all)
        collect_ions(pot_indices, xyz_frame, cx, cy, pot_z_all)
        collect_ions(cal_indices, xyz_frame, cx, cy, cal_z_all)

    # 6. Calculation (Counts per bin averaged over frames)
    # Create bins along Z-axis (0.5 Angstrom width) covering the detected range
    bin_width = 0.5 
    bins = np.arange(z_min_angstrom, z_max_angstrom + bin_width, bin_width)
    
    # Calculate bin centers
    bin_centers = (bins[:-1] + bins[1:]) / 2

    def get_avg_counts(data_list):
        if not data_list:
            return np.zeros(len(bin_centers))
        hist, _ = np.histogram(data_list, bins=bins)
        return hist / traj.n_frames

    avg_sod = get_avg_counts(sod_z_all)
    avg_pot = get_avg_counts(pot_z_all)
    avg_cal = get_avg_counts(cal_z_all)

    # 7. Save Data
    # Stack data columns for output
    header = "Z(A)    Avg_Na    Avg_K    Avg_Cl"
    data_block = np.column_stack((bin_centers, avg_sod, avg_pot, avg_cal))
    
    np.savetxt(out_file, data_block, header=header, fmt='%.3f', comments='')
    
    end_time = timer()
    print(f"Analysis complete in {end_time - start_time:.2f} seconds. Data saved to {out_file}")

if __name__ == "__main__":
    main()