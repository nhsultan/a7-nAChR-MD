#!/usr/bin/env python3
"""
calc_pca_tsne.py

Description:
This script performs dimensionality reduction (PCA and t-SNE) on the Transmembrane 
Domain (TMD) of the α7-nAChR to classify structural states. 

It calculates 10 pairwise distances between the 5 subunits for 7 specific 
lumen-facing residues, generating 70 distance features per frame. These features 
are then standardized (mean=0, variance=1) before being reduced to 2 dimensions 
via PCA and t-SNE.

Outputs:
A single CSV file containing the Frame index, PCA components (1 & 2), and 
t-SNE components (1 & 2).
"""

import argparse
import numpy as np
import pandas as pd
import mdtraj as md
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from timeit import default_timer as timer
import warnings

# Suppress minor warnings
warnings.filterwarnings('ignore')

def main():
    start_time = timer()

    # 1. Setup Argument Parsing
    parser = argparse.ArgumentParser(description="Calculate PCA and t-SNE on TMD distances")
    parser.add_argument("--traj", required=True, help="Path to the trajectory file (.nc, .xtc)")
    parser.add_argument("--top", required=True, help="Path to the topology file (.psf, .pdb)")
    parser.add_argument("--out_prefix", required=True, help="Prefix for output files (e.g., results/pca)")
    args = parser.parse_args()

    out_file = f"{args.out_prefix}_dimensionality.csv"

    print(f"Loading trajectory: {args.traj}")
    print(f"Loading topology: {args.top}")
    traj = md.load(args.traj, top=args.top)
    top = traj.topology
    print(f"Loaded {traj.n_frames} frames.")

    # 2. Define TMD Residues and Subunit Offsets
    # Homopentamer offset assumes subunits are exactly 386 residues apart
    base_res_ids = [237, 240, 244, 247, 251, 254, 258]
    atom_names = ['CD', 'CB', 'CB', 'CG', 'CB', 'CG', 'CD']
    subunit_offsets = [0, 386, 772, 1158, 1544]

    # Pairwise definitions for the 5 subunits (A=0, B=1, C=2, D=3, E=4)
    pair_defs = [
        (0, 1, "ab"), (1, 2, "bc"), (2, 3, "cd"), (3, 4, "de"), (4, 0, "ea"),
        (1, 4, "be"), (1, 3, "bd"), (0, 3, "ad"), (0, 2, "ac"), (2, 4, "ce")
    ]

    atom_pairs = []
    feature_labels = []

    print("Identifying atoms for distance calculations...")
    
    # Extract specific atom indices for all 5 subunits for each target residue
    for base_id, atom_name in zip(base_res_ids, atom_names):
        subunit_atom_indices = []
        
        for offset in subunit_offsets:
            # MDTraj uses 1-based indexing for resSeq (matches Amber resid)
            target_resSeq = base_id + offset
            query = f"resSeq {target_resSeq} and name {atom_name}"
            
            # Find the atom index
            idx = top.select(query)
            if len(idx) != 1:
                print(f"Warning: Expected 1 atom for {query}, but found {len(idx)}. Check topology.")
            subunit_atom_indices.append(idx[0])

        # Generate the 10 pairwise combinations for this residue
        for p1, p2, label in pair_defs:
            atom_pairs.append([subunit_atom_indices[p1], subunit_atom_indices[p2]])
            feature_labels.append(f"Res{base_id}_{atom_name}_{label}")

    # 3. Compute Distances
    print(f"Calculating {len(atom_pairs)} pairwise distances across {traj.n_frames} frames...")
    # md.compute_distances automatically handles periodic boundary conditions (PBC)
    # It returns an array of shape (n_frames, n_pairs) in nanometers
    distances_nm = md.compute_distances(traj, atom_pairs)
    
    # Convert to Angstroms to match original cpptraj output
    distances_angstrom = distances_nm * 10.0

    # 4. Standardization
    print("Standardizing features (mean=0, variance=1)...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(distances_angstrom)

    # 5. Dimensionality Reduction (PCA & t-SNE)
    print("Running Principal Component Analysis (PCA)...")
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    print("Running t-Distributed Stochastic Neighbor Embedding (t-SNE)...")
    # Note: t-SNE is computationally heavy. random_state ensures reproducibility.
    tsne = TSNE(n_components=2, random_state=42)
    X_tsne = tsne.fit_transform(X_scaled)

    # 6. Save Data for MATLAB Plotting
    print(f"Saving coordinates to {out_file}...")
    df = pd.DataFrame({
        'Frame': np.arange(traj.n_frames),
        'PCA1': X_pca[:, 0],
        'PCA2': X_pca[:, 1],
        'tSNE1': X_tsne[:, 0],
        'tSNE2': X_tsne[:, 1]
    })
    
    # Save as a clean CSV
    df.to_csv(out_file, index=False, float_format='%.4f')

    end_time = timer()
    print(f"Dimensionality reduction complete in {end_time - start_time:.2f} seconds!")

if __name__ == "__main__":
    main()