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

def process_scRNAseq(input_file, output_file, min_genes=200, min_cells=3, 
                     mt_pattern="^MT-", mt_cutoff=20, target_sum=1e4, 
                     n_hvg=2000, n_pcs=50, resolution=0.5, 
                     generate_plots=True, output_dir=None):
    """
    Process and analyze cancer scRNA-seq data from 10X Chromium.
    
    Steps:
    1. Load raw data and calculate QC metrics
    2. Filter cells and genes based on QC parameters
    3. Normalize, log-transform and identify highly variable genes
    4. Perform dimensionality reduction and clustering
    5. Identify marker genes and produce visualizations
    6. Save processed data with annotations
    
    Parameters:
    - input_file: Path to the raw 10X h5 file.
    - output_file: Path to save the processed AnnData object (.h5ad).
    - min_genes: Minimum number of genes per cell (default: 200).
    - min_cells: Minimum cells per gene (default: 3).
    - mt_pattern: Pattern to identify mitochondrial genes (default: "^MT-").
    - mt_cutoff: Maximum percentage of mitochondrial genes (default: 20).
    - target_sum: Target counts per cell for normalization (default: 1e4).
    - n_hvg: Number of highly variable genes to select (default: 2000).
    - n_pcs: Number of principal components to use (default: 50).
    - resolution: Resolution parameter for Leiden clustering (default: 0.5).
    - generate_plots: Whether to generate QC and analysis plots (default: True).
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
    
    # Load data
    try:
        print(f"Loading raw data from {input_file} ...")
        adata = sc.read_10x_h5(input_file)
        print(f"Data loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")
    except Exception as e:
        print("Error reading input file:", e)
        sys.exit(1)
    
    print("Calculating quality control metrics ...")
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith(mt_pattern)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    # Generate QC plots
    if generate_plots:
        print("Generating QC plots...")
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
                    jitter=0.4, multi_panel=True, save='_qc_metrics_before_filtering.pdf')
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_mt_vs_counts_before_filtering.pdf')
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_genes_vs_counts_before_filtering.pdf')
    
    print("Performing quality control filtering ...")
    # Filter cells based on QC metrics
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs.pct_counts_mt < mt_cutoff]
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    print(f"After filtering: {adata.shape[0]} cells x {adata.shape[1]} genes")
    
    if generate_plots:
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                    jitter=0.4, multi_panel=True, save='_qc_metrics_after_filtering.pdf')
    
    # Normalization and log transformation
    print(f"Normalizing data to a target sum of {target_sum} and log-transforming...")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    
    # Feature selection - identify highly variable genes
    print(f"Identifying {n_hvg} highly variable genes...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=n_hvg)
    
    if generate_plots:
        sc.pl.highly_variable_genes(adata, save='_highly_variable_genes.pdf')
    
    # Keep only HVGs for dimensionality reduction
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    
    # Scale data
    print("Scaling data...")
    sc.pp.scale(adata, max_value=10)
    
    # Run PCA
    print(f"Running PCA with {n_pcs} components...")
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs)
    
    if generate_plots:
        sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save='_pca_variance_ratio.pdf')
        sc.pl.pca(adata, color='total_counts', save='_pca.pdf')
    
    # Compute neighborhood graph
    print("Computing neighborhood graph...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_pcs)
    
    # Run UMAP for visualization
    print("Running UMAP...")
    sc.tl.umap(adata)
    
    # Run clustering
    print(f"Running Leiden clustering at resolution {resolution}...")
    sc.tl.leiden(adata, resolution=resolution)
    
    if generate_plots:
        sc.pl.umap(adata, color=['leiden', 'total_counts', 'n_genes_by_counts', 'pct_counts_mt'], 
                  ncols=2, save='_umap_clustering.pdf')
    
    # Find marker genes for each cluster
    print("Identifying marker genes for each cluster...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
    
    if generate_plots:
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes.pdf')
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=10, groupby='leiden', 
                                       show_gene_labels=True, save='_marker_genes_heatmap.pdf')
    
    # Calculate cell cycle scores (important for cancer analysis)
    print("Calculating cell cycle scores...")
    # Cell cycle gene sets from Tirosh et al.
    s_genes = ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 
               'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'HELLS', 'RFC2', 'RPA2', 'NASP', 
               'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 
               'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 
               'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
    g2m_genes = ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 
                 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 
                 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 
                 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 
                 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 
                 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 
                 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
    
    # Gene sets might need to be filtered to those present in the dataset
    s_genes = [gene for gene in s_genes if gene in adata.var_names]
    g2m_genes = [gene for gene in g2m_genes if gene in adata.var_names]
    
    if len(s_genes) > 0 and len(g2m_genes) > 0:
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
        
        if generate_plots:
            sc.pl.umap(adata, color=['S_score', 'G2M_score', 'phase'], 
                      ncols=3, save='_cell_cycle.pdf')
    else:
        print("Warning: Not enough cell cycle genes found in dataset. Skipping cell cycle scoring.")
    
    # Save processed data
    print(f"Saving processed data to {output_file} ...")
    try:
        adata.write(output_file)
    except Exception as e:
        print("Error saving processed data:", e)
        sys.exit(1)
    
    print("Analysis completed successfully.")
    return adata

def main():
    parser = argparse.ArgumentParser(description="Process and analyze cancer scRNA-seq data.")
    parser.add_argument("input_file", help="Path to the raw 10X h5 file.")
    parser.add_argument("output_file", help="Path to save the processed AnnData object (.h5ad).")
    parser.add_argument("--min_genes", type=int, default=200, help="Minimum genes per cell (default: 200)")
    parser.add_argument("--min_cells", type=int, default=3, help="Minimum cells per gene (default: 3)")
    parser.add_argument("--mt_pattern", type=str, default="^MT-", help="Pattern to identify mitochondrial genes (default: '^MT-')")
    parser.add_argument("--mt_cutoff", type=float, default=20, help="Maximum percentage of mitochondrial genes (default: 20)")
    parser.add_argument("--target_sum", type=float, default=1e4, help="Normalization target sum per cell (default: 1e4)")
    parser.add_argument("--n_hvg", type=int, default=2000, help="Number of highly variable genes to select (default: 2000)")
    parser.add_argument("--n_pcs", type=int, default=50, help="Number of principal components to use (default: 50)")
    parser.add_argument("--resolution", type=float, default=0.5, help="Resolution parameter for Leiden clustering (default: 0.5)")
    parser.add_argument("--no_plots", action="store_false", dest="generate_plots", help="Disable generation of plots")
    parser.add_argument("--output_dir", type=str, default=None, help="Directory to save plots (default: same as output_file)")
    args = parser.parse_args()

    process_scRNAseq(
        args.input_file, args.output_file, 
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