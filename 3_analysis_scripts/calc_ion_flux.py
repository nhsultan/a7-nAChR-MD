#!/usr/bin/env python3
"""
calc_ion_flux.py

Description:
This script calculates the ion flux (permeation events) through the pore of the 
α7-nAChR. It processes a molecular dynamics trajectory and tracks the movement of 
specific ions (POT, SOD, CLA). 

It defines the channel pore geometry using specific residues:
  - Z-axis reference (z=0): VAL 250 (Equatorial gate)
  - Bottom boundary: GLU 236 (Intracellular end)
  - Top boundary: GLU 257 (Extracellular end)

The script categorizes permeation events as:
  - TB: Top to Bottom permeation
  - LB: Lateral to Bottom permeation (entering through fenestrations)
  - BT: Bottom to Top permeation
  - BL: Bottom to Lateral permeation

Outputs:
Saves three text files containing the tracked ion events, pathway types, and ion identities.
"""

import argparse
import numpy as np
import mdtraj as md
from timeit import default_timer as timer

def main():
    start_time = timer()

    # 1. Setup Argument Parsing to match the master bash script
    parser = argparse.ArgumentParser(description="Calculate ion flux from MD trajectory")
    parser.add_argument("--traj", required=True, help="Path to the trajectory file (e.g., .nc, .xtc)")
    parser.add_argument("--top", required=True, help="Path to the topology file (e.g., .psf, .pdb)")
    parser.add_argument("--out_prefix", required=True, help="Prefix/path for output files (e.g., results/flux)")
    args = parser.parse_args()

    # Define output filenames based on prefix
    event_data_file = f"{args.out_prefix}_IonTracking.txt"
    pathway_file = f"{args.out_prefix}_IonPathway.txt"
    which_ion_file = f"{args.out_prefix}_WhichIon.txt"

    print(f"Loading trajectory: {args.traj}")
    print(f"Loading topology: {args.top}")
    traj = md.load(args.traj, top=args.top)
    top = traj.topology

    # 2. System Selection and Geometric Boundaries
    print("Calculating pore boundaries and center axis...")
    
    # Calculate global X and Y average of the protein to define the pore center axis
    protein_idx = top.select('protein')
    avg_x = np.mean(traj.xyz[:, protein_idx, 0])
    avg_y = np.mean(traj.xyz[:, protein_idx, 1])

    # Define Z-axis boundaries based on specific residues
    val250_idx = top.select('resname VAL and resid 250')
    val_z_avg = np.mean(traj.xyz[:, val250_idx, 2])  # This acts as z = 0

    glu236_bot_idx = top.select('resname GLU and resid 236')
    glu_bot_z_avg = np.mean(traj.xyz[:, glu236_bot_idx, 2])
    glu_bot_norm = glu_bot_z_avg - val_z_avg         # Negative Z boundary

    # Note: Original script named this ARGtop but selected GLU 257. Updated name for clarity.
    glu257_top_idx = top.select('resname GLU and resid 257')
    glu_top_z_avg = np.mean(traj.xyz[:, glu257_top_idx, 2])
    glu_top_norm = glu_top_z_avg - val_z_avg         # Positive Z boundary

    # 3. Track Ions
    target_ions = ['POT', 'SOD', 'CLA']
    
    total_events = 0
    pathway_list = []
    ion_number_list = []
    event_start_list = []
    event_end_list = []
    events_count_list = []
    ion_identity_list = []

    print(f"Tracking ions: {target_ions}...")

    for ion_name in target_ions:
        md_ion_indices = top.select(f'resname {ion_name}')

        for ion_idx in md_ion_indices:
            coords = traj.xyz[:, ion_idx, :]
            
            # Vectorized calculation for radius and normalized Z
            x_norm = coords[:, 0] - avg_x
            y_norm = coords[:, 1] - avg_y
            r_array = np.sqrt(x_norm**2 + y_norm**2)
            z_norm_array = coords[:, 2] - val_z_avg

            # Create a mask for when the ion is strictly inside the pore cylinder
            in_pore_mask = (z_norm_array > glu_bot_norm) & (z_norm_array < glu_top_norm) & (r_array < 2.0)
            
            radius_tracked = np.where(in_pore_mask, r_array, None)
            z_tracked = np.where(in_pore_mask, z_norm_array, None)

            # 4. Identify Continuous Permeation Events
            event_starts = []
            event_ends = []
            
            for i in range(len(z_tracked)):
                if radius_tracked[i] is not None:
                    # Check for event start/end at the very boundaries of the trajectory
                    if i == 0 and radius_tracked[i+1] is not None:
                        event_starts.append(i)
                    if i == len(z_tracked) - 1 and radius_tracked[i-1] is not None:
                        event_ends.append(i)

                    # Check for event start/end during the trajectory
                    if 0 < i < len(z_tracked) - 1:
                        if radius_tracked[i-1] is None and radius_tracked[i+1] is not None:
                            event_starts.append(i)
                        if radius_tracked[i+1] is None and radius_tracked[i-1] is not None:
                            event_ends.append(i)

            # 5. Classify the Type of Event (Positive vs Negative flux)
            event_type = []
            event_start_time = []
            event_end_time = []
            events_count = []
            ion_list = []

            for start, end in zip(event_starts, event_ends):
                # Ensure we are not just looking at a boundary artifact
                if start != 0 and end != len(z_tracked) - 1:
                    
                    z_start = z_tracked[start]
                    z_end = z_tracked[end]
                    r_start = radius_tracked[start]
                    r_end = radius_tracked[end]

                    # Check for Top-to-Bottom (Inward flux)
                    if z_end < 0 and z_start > 0:
                        if r_start <= 1.8:
                            event_type.append('TB')
                        else:
                            event_type.append('LB') # Lateral entry
                            
                        events_count.append(1)
                        ion_list.append(ion_name)
                        event_start_time.append(start)
                        event_end_time.append(end)
                        total_events += 1

                    # Check for Bottom-to-Top (Outward flux)
                    elif z_end > 0 and z_start < 0:
                        if r_end < 2.0:
                            event_type.append('BT')
                        else:
                            event_type.append('BL') # Lateral exit
                            
                        events_count.append(-1)
                        ion_list.append(ion_name)
                        event_start_time.append(start)
                        event_end_time.append(end)
                        total_events += 1

            if event_type:
                pathway_list.append(event_type)
                ion_number_list.append(ion_idx)
                event_start_list.append(event_start_time)
                event_end_list.append(event_end_time)
                events_count_list.append(events_count)
                ion_identity_list.append(ion_list)

    # 6. Format and Save Output Data
    print(f"Total valid permeation events detected: {total_events}")
    
    if total_events > 0:
        event_data_combined = np.zeros((total_events, 4))
        pathway_type_out = []
        which_ion_out = []
        row_idx = 0

        for i in range(len(pathway_list)):
            for j in range(len(pathway_list[i])):
                event_data_combined[row_idx, 0] = ion_number_list[i]
                event_data_combined[row_idx, 1] = events_count_list[i][j]
                event_data_combined[row_idx, 2] = event_start_list[i][j]
                event_data_combined[row_idx, 3] = event_end_list[i][j]
                pathway_type_out.append(pathway_list[i][j])
                which_ion_out.append(ion_identity_list[i][j])
                row_idx += 1

        # Save to text files
        np.savetxt(event_data_file, event_data_combined, fmt='%.0f')

        with open(pathway_file, 'w') as fp:
            fp.write('\n'.join(pathway_type_out))

        with open(which_ion_file, 'w') as fp:
            fp.write('\n'.join(which_ion_out))
            
        print(f"Data successfully written to {args.out_prefix}_* files.")
    else:
        print("No events detected. Output files were not generated.")

    end_time = timer()
    print(f"Analysis completed in {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()