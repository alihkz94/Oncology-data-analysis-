#!/usr/bin/env python

import argparse
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
from scipy.sparse import issparse

def process_spatial(input_dir, output_file, min_genes=100, min_cells=3, 
                    mt_pattern="^MT-", mt_cutoff=20, target_sum=1e4, 
                    n_hvg=2000, n_pcs=30, resolution=0.8, 
                    generate_plots=True, output_dir=None):
    """
    Process and analyze spatial transcriptomics data (Visium) for tumor microenvironment mapping.
    
    Steps:
    1. Load spatial data and calculate QC metrics.
    2. Filter spots and genes based on QC parameters.
    3. Normalize, log-transform and identify highly variable genes.
    4. Perform dimensionality reduction and clustering.
    5. Generate spatial visualization plots.
    6. Save processed data with spatial annotations.
    
    Parameters:
    - input_dir: Path to the spatial data directory (Visium format).
    - output_file: Path to save the processed AnnData object (.h5ad).
    - min_genes: Minimum number of genes per spot (default: 100).
    - min_cells: Minimum spots per gene (default: 3).
    - mt_pattern: Pattern to identify mitochondrial genes (default: "^MT-").
    - mt_cutoff: Maximum percentage of mitochondrial genes (default: 20).
    - target_sum: Target counts per spot for normalization (default: 1e4).
    - n_hvg: Number of highly variable genes to select (default: 2000).
    - n_pcs: Number of principal components to use (default: 30).
    - resolution: Resolution parameter for Leiden clustering (default: 0.8).
    - generate_plots: Whether to generate QC and spatial analysis plots (default: True).
    - output_dir: Directory to save visualization plots (default: same as output_file).
    """
    # Set up the output directory for plots
    if generate_plots:
        if output_dir is None:
            output_dir = os.path.dirname(output_file)
            if not output_dir:
                output_dir = '.'
        os.makedirs(output_dir, exist_ok=True)
        sc.settings.figdir = output_dir
    
    # Load spatial data (Visium format)
    try:
        print(f"Loading spatial data from {input_dir} ...")
        adata = sc.read_visium(path=input_dir)
        print(f"Data loaded: {adata.shape[0]} spots x {adata.shape[1]} genes")
    except Exception as e:
        print("Error reading spatial data:", e)
        sys.exit(1)
    
    print("Calculating quality control metrics ...")
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith(mt_pattern)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Generate QC plots if enabled
    if generate_plots:
        print("Generating QC plots for spatial data...")
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                     jitter=0.4, multi_panel=True, save='_spatial_qc_metrics_before_filtering.pdf')
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', 
                      save='_spatial_mt_vs_counts_before_filtering.pdf')
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', 
                      save='_spatial_genes_vs_counts_before_filtering.pdf')
    
    print("Performing quality control filtering ...")
    # Filter spots based on QC metrics
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs.pct_counts_mt < mt_cutoff]
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After filtering: {adata.shape[0]} spots x {adata.shape[1]} genes")
    
    if generate_plots:
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                     jitter=0.4, multi_panel=True, save='_spatial_qc_metrics_after_filtering.pdf')
    
    # Normalization and log transformation
    print(f"Normalizing data to a target sum of {target_sum} and log-transforming...")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    
    # Identify highly variable genes
    print(f"Identifying {n_hvg} highly variable genes...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_hvg)
    
    if generate_plots:
        sc.pl.highly_variable_genes(adata, save='_spatial_highly_variable_genes.pdf')
    
    # Keep only HVGs for further analysis
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    
    # Scale data
    print("Scaling data...")
    sc.pp.scale(adata, max_value=10)
    
    # Run PCA
    print(f"Running PCA with {n_pcs} components...")
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)
    
    if generate_plots:
        sc.pl.pca_variance_ratio(adata, n_pcs=n_pcs, log=True, save='_spatial_pca_variance_ratio.pdf')
        sc.pl.pca(adata, color='total_counts', save='_spatial_pca.pdf')
    
    # Compute neighborhood graph
    print("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs)
    
    # Run UMAP for visualization
    print("Running UMAP...")
    sc.tl.umap(adata)
    
    # Run Leiden clustering
    print(f"Running Leiden clustering at resolution {resolution}...")
    sc.tl.leiden(adata, resolution=resolution)
    
    # Generate spatial visualization plots
    if generate_plots:
        print("Generating spatial plots...")
        sc.pl.spatial(adata, color=['leiden', 'total_counts'], ncols=2, save='_spatial_clusters.pdf')
        # Optional: if high-resolution images are available in the dataset
        sc.pl.spatial(adata, img_key="hires", color='leiden', save='_spatial_image_clusters.pdf')
    
    # Save processed data
    print(f"Saving processed spatial data to {output_file} ...")
    try:
        adata.write(output_file)
    except Exception as e:
        print("Error saving processed spatial data:", e)
        sys.exit(1)
    
    print("Spatial transcriptomics analysis completed successfully.")
    return adata

def main():
    parser = argparse.ArgumentParser(description="Process and analyze spatial transcriptomics data (Visium).")
    parser.add_argument("input_dir", help="Path to the spatial data directory (Visium format).")
    parser.add_argument("output_file", help="Path to save the processed AnnData object (.h5ad).")
    parser.add_argument("--min_genes", type=int, default=100, help="Minimum genes per spot (default: 100)")
    parser.add_argument("--min_cells", type=int, default=3, help="Minimum spots per gene (default: 3)")
    parser.add_argument("--mt_pattern", type=str, default="^MT-", help="Pattern to identify mitochondrial genes (default: '^MT-')")
    parser.add_argument("--mt_cutoff", type=float, default=20, help="Maximum percentage of mitochondrial genes (default: 20)")
    parser.add_argument("--target_sum", type=float, default=1e4, help="Normalization target sum per spot (default: 1e4)")
    parser.add_argument("--n_hvg", type=int, default=2000, help="Number of highly variable genes to select (default: 2000)")
    parser.add_argument("--n_pcs", type=int, default=30, help="Number of principal components to use (default: 30)")
    parser.add_argument("--resolution", type=float, default=0.8, help="Resolution parameter for Leiden clustering (default: 0.8)")
    parser.add_argument("--no_plots", action="store_false", dest="generate_plots", help="Disable generation of plots")
    parser.add_argument("--output_dir", type=str, default=None, help="Directory to save plots (default: same as output_file)")
    args = parser.parse_args()

    process_spatial(
        args.input_dir, args.output_file, 
        min_genes=args.min_genes, 
        min_cells=args.min_cells,
        mt_pattern=args.mt_pattern,
        mt_cutoff=args.mt_cutoff,
        target_sum=args.target_sum,
        n_hvg=args.n_hvg,
        n_pcs=args.n_pcs,
        resolution=args.resolution,
        generate_plots=args.generate_plots,
        output_dir=args.output_dir
    )

if __name__ == "__main__":
    main()
