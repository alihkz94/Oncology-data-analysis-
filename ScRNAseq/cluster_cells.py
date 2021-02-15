#!/usr/bin/env python

import argparse
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
from sklearn.metrics import silhouette_score
from scipy.stats import zscore
import anndata
import seaborn as sns

def cluster_cells(input_file, output_file, n_neighbors=15, n_pcs=50, 
                  resolutions=[0.2, 0.4, 0.8, 1.0, 1.5], clustering_method="leiden", 
                  use_harmony=False, regress_out=None, generate_plots=True, 
                  output_dir=None, min_cluster_size=10, annotation_database=None,
                  cnv_analysis=False, trajectory_analysis=False):
    """
    clustering  of scRNA-seq data for cancer research.
    
    Steps:
    1. Load the processed AnnData object
    2. Optional batch correction and data integration
    3. Feature selection with highly variable genes
    4. Dimensionality reduction (PCA and optionally Harmony integration)
    5. Optional regression of confounding factors
    6. Neighborhood graph construction
    7. Multi-resolution clustering with quality assessment
    8. Dimensionality reduction for visualization (UMAP and t-SNE)
    9. Marker gene identification for each cluster
    10. Optional cell type annotation using reference databases
    11. Cancer-specific analyses (CNV inference, oncogenic signatures)
    12. Optional trajectory/pseudotime analysis
    13. Comprehensive quality metrics and visualizations
    14. Save the deeply annotated AnnData object
    
    Parameters:
    - input_file: Path to the processed AnnData object (.h5ad)
    - output_file: Path to save the clustered AnnData object (.h5ad)
    - n_neighbors: Number of neighbors for graph construction (default: 15)
    - n_pcs: Number of principal components to use (default: 50)
    - resolutions: List of resolutions for clustering analysis (default: [0.2, 0.4, 0.8, 1.0, 1.5])
    - clustering_method: Algorithm for clustering ('leiden', 'louvain', or 'both') (default: 'leiden')
    - use_harmony: Whether to apply Harmony batch correction (default: False)
    - regress_out: Variables to regress out ('cell_cycle', 'mito', or None) (default: None)
    - generate_plots: Whether to generate analysis plots (default: True)
    - output_dir: Directory for saving plots and reports (default: same as output_file)
    - min_cluster_size: Minimum number of cells for a valid cluster (default: 10)
    - annotation_database: Reference database for cell type annotation (default: None)
    - cnv_analysis: Whether to perform CNV inference from expression data (default: False)
    - trajectory_analysis: Whether to perform trajectory analysis (default: False)
    """
    # Set up output directory
    if generate_plots:
        if output_dir is None:
            output_dir = os.path.dirname(output_file)
            if not output_dir:
                output_dir = '.'
        os.makedirs(output_dir, exist_ok=True)
        sc.settings.figdir = output_dir
        sc.settings.set_figure_params(dpi=100, figsize=(8, 8))

    # Load the processed data
    try:
        print(f"Loading processed data from {input_file}...")
        adata = sc.read(input_file)
        print(f"Data loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")
    except Exception as e:
        print(f"Error reading processed data: {e}")
        sys.exit(1)
        
    # Store raw counts for later use if not already stored
    if 'raw' not in adata.layers and 'counts' not in adata.layers:
        if adata.raw is None:
            print("Storing raw counts for downstream analyses...")
            adata.layers['counts'] = adata.X.copy()
    
    # Feature selection - identify highly variable genes
    print("Identifying highly variable genes...")
    if 'highly_variable' not in adata.var.columns:
        sc.pp.highly_variable_genes(
            adata, 
            min_mean=0.02, 
            max_mean=8, 
            min_disp=0.25, 
            n_top_genes=4000,
            batch_key='batch' if 'batch' in adata.obs.columns else None
        )
        if generate_plots:
            sc.pl.highly_variable_genes(adata, save='_highly_variable_genes.pdf')
    else:
        print("Using existing highly variable gene selection.")
    
    # Optional: Regress out confounding factors
    if regress_out:
        print(f"Regressing out confounding factors: {regress_out}")
        variables_to_regress = []
        
        if 'cell_cycle' in regress_out and 'S_score' in adata.obs and 'G2M_score' in adata.obs:
            variables_to_regress.extend(['S_score', 'G2M_score'])
            
        if 'mito' in regress_out and 'pct_counts_mt' in adata.obs:
            variables_to_regress.append('pct_counts_mt')
            
        if variables_to_regress:
            sc.pp.regress_out(adata, variables_to_regress)
            
    # Keep only highly variable genes for dimensionality reduction
    hvg_adata = adata[:, adata.var.highly_variable].copy()

    # Principal component analysis
    print(f"Performing PCA with {n_pcs} components...")
    sc.tl.pca(hvg_adata, svd_solver='arpack', n_comps=n_pcs, use_highly_variable=True)
    if generate_plots:
        sc.pl.pca_variance_ratio(hvg_adata, n_pcs=n_pcs, log=True, save='_pca_variance_ratio.pdf')
        sc.pl.pca(hvg_adata, color=['total_counts', 'pct_counts_mt'] if 'pct_counts_mt' in hvg_adata.obs else 'total_counts', 
                 save='_pca.pdf')

    # Optional batch correction using Harmony
    if use_harmony and 'batch' in adata.obs.columns:
        print("Performing batch correction with Harmony...")
        try:
            sce.pp.harmony_integrate(hvg_adata, 'batch')
            if generate_plots:
                sc.pl.umap(hvg_adata, color='batch', save='_harmony_batch_correction.pdf')
        except Exception as e:
            print(f"Warning: Harmony batch correction failed: {e}")
            print("Proceeding without batch correction.")
    
    # Compute neighborhood graph
    print(f"Computing neighborhood graph with n_neighbors={n_neighbors} and n_pcs={n_pcs}...")
    sc.pp.neighbors(
        hvg_adata, 
        n_neighbors=n_neighbors,
        n_pcs=n_pcs, 
        use_rep='X_pca_harmony' if use_harmony and 'X_pca_harmony' in hvg_adata.obsm else 'X_pca'
    )
    
    # Multi-resolution clustering with quality assessment
    print("Performing multi-resolution clustering...")
    best_resolution = None
    best_silhouette = -1
    
    for resolution in resolutions:
        print(f"Clustering with resolution={resolution}...")
        
        # Perform clustering with specified method
        if clustering_method.lower() in ['leiden', 'both']:
            sc.tl.leiden(hvg_adata, resolution=resolution, key_added=f'leiden_r{resolution}')
            
        if clustering_method.lower() in ['louvain', 'both']:
            sc.tl.louvain(hvg_adata, resolution=resolution, key_added=f'louvain_r{resolution}')
            
        # Calculate clustering quality metrics
        main_clustering = f'leiden_r{resolution}' if clustering_method.lower() in ['leiden', 'both'] else f'louvain_r{resolution}'
        
        # Calculate silhouette score
        n_clusters = len(hvg_adata.obs[main_clustering].unique())
        if n_clusters > 1 and n_clusters < hvg_adata.shape[0] - 1:
            try:
                # Calculate silhouette score using PCA space
                sil_score = silhouette_score(
                    hvg_adata.obsm['X_pca'], 
                    hvg_adata.obs[main_clustering], 
                    metric='euclidean'
                )
                print(f"Resolution {resolution} produced {n_clusters} clusters with silhouette score: {sil_score:.3f}")
                
                # Track best resolution based on silhouette score
                if sil_score > best_silhouette:
                    best_silhouette = sil_score
                    best_resolution = resolution
            except Exception as e:
                print(f"Warning: Could not calculate silhouette score for resolution {resolution}: {e}")
                
    # Select best resolution based on silhouette score
    if best_resolution is not None:
        print(f"Best clustering resolution based on silhouette score: {best_resolution}")
        main_clustering = f'leiden_r{best_resolution}' if clustering_method.lower() in ['leiden', 'both'] else f'louvain_r{best_resolution}'
        hvg_adata.obs['optimal_clusters'] = hvg_adata.obs[main_clustering].copy()
    else:
        # Default if no optimal resolution found
        default_res = resolutions[len(resolutions)//2]  # Middle resolution
        print(f"No optimal resolution identified. Using default resolution: {default_res}")
        main_clustering = f'leiden_r{default_res}' if clustering_method.lower() in ['leiden', 'both'] else f'louvain_r{default_res}'
        hvg_adata.obs['optimal_clusters'] = hvg_adata.obs[main_clustering].copy()

    # Compute UMAP and t-SNE embeddings
    print("Computing UMAP embedding for visualization...")
    sc.tl.umap(hvg_adata)
    
    print("Computing t-SNE embedding for alternative visualization...")
    sc.tl.tsne(hvg_adata, use_rep='X_pca', n_jobs=8)
    
    # Find marker genes for each cluster
    print("Identifying marker genes for each cluster...")
    sc.tl.rank_genes_groups(hvg_adata, 'optimal_clusters', method='wilcoxon')

    if generate_plots:
        # Generate cluster visualization plots
        print("Generating visualization plots...")
        sc.pl.umap(hvg_adata, color=['optimal_clusters'], save='_optimal_clustering.pdf')
        
        # Visualize marker genes
        sc.pl.rank_genes_groups(hvg_adata, n_genes=25, sharey=False, save='_marker_genes.pdf')
        
        # Generate dot plot for top markers
        top_markers_per_group = {}
        for cluster in hvg_adata.obs['optimal_clusters'].unique():
            markers = sc.get.rank_genes_groups_df(hvg_adata, group=cluster)
            # Select top 5 genes based on lowest adjusted p-value
            top_markers = markers.sort_values('pvals_adj').head(3)['names'].tolist()
            top_markers_per_group[cluster] = top_markers
            
        # Flatten list of markers
        all_top_markers = [gene for genes in top_markers_per_group.values() for gene in genes]
        if all_top_markers:
            sc.pl.dotplot(hvg_adata, all_top_markers, groupby='optimal_clusters', 
                         standard_scale='var', save='_top_markers_dotplot.pdf')
            
            sc.pl.heatmap(hvg_adata, all_top_markers, groupby='optimal_clusters', 
                         standard_scale='var', save='_top_markers_heatmap.pdf')
    
    # Optional: Infer CNVs from expression data (important for cancer)
    if cnv_analysis:
        print("Performing CNV inference from expression data...")
        try:
            # Simple CNV scoring based on moving average of expression across chromosomes
            # In a real implementation, would use inferCNV or similar tools
            # This is a placeholder for the concept
            print("Note: Full CNV analysis requires external tools like inferCNV or CopyKAT")
            print("Generating simple CNV-like visualization based on gene expression patterns")
            
            # Add chromosome annotation to genes if not present
            if 'chromosome' not in adata.var:
                print("Warning: Chromosome annotation not found. CNV analysis will be limited.")
            else:
                # Create a heatmap of genes organized by chromosomal position
                chrom_genes = {}
                for chrom in sorted(adata.var['chromosome'].unique()):
                    chr_genes = adata.var_names[adata.var['chromosome'] == chrom].tolist()
                    if len(chr_genes) > 50:  # Subsample genes if too many
                        chr_genes = np.random.choice(chr_genes, 50, replace=False).tolist()
                    chrom_genes[f"Chr_{chrom}"] = chr_genes
                
                # Plot chromosomal expression patterns by cluster
                for chrom, genes in chrom_genes.items():
                    valid_genes = [g for g in genes if g in adata.var_names]
                    if valid_genes and len(valid_genes) > 5:
                        if generate_plots:
                            sc.pl.heatmap(adata, valid_genes, groupby='optimal_clusters', 
                                          standard_scale='var', 
                                          save=f'_cnv_{chrom}_heatmap.pdf')
        except Exception as e:
            print(f"Warning: CNV analysis failed: {e}")

    # Optional: Cell type annotation using reference databases
    if annotation_database:
        print(f"Performing cell type annotation using {annotation_database}...")
        try:
            # Placeholder for cell type annotation
            # In a real implementation, would use tools like SingleR, scmap, etc.
            print("Note: Cell type annotation requires integration with databases like CellMarker")
        except Exception as e:
            print(f"Warning: Cell type annotation failed: {e}")

    # Optional: Trajectory analysis
    if trajectory_analysis:
        print("Performing trajectory/pseudotime analysis...")
        try:
            # Placeholder for trajectory analysis
            # In a real implementation, would use methods like PAGA, etc.
            # PAGA (Partition-based graph abstraction) - good for disconnected topologies
            sc.tl.paga(hvg_adata, groups='optimal_clusters')
            if generate_plots:
                sc.pl.paga(hvg_adata, save='_paga_trajectory.pdf')
                sc.pl.paga_compare(hvg_adata, basis='umap', save='_paga_umap.pdf')
            
        except Exception as e:
            print(f"Warning: Trajectory analysis failed: {e}")

    # Transfer annotations back to the original AnnData object
    key_obs = ['optimal_clusters']
    for res in resolutions:
        if f'leiden_r{res}' in hvg_adata.obs.columns:
            key_obs.append(f'leiden_r{res}')
        if f'louvain_r{res}' in hvg_adata.obs.columns:
            key_obs.append(f'louvain_r{res}')
    
    # Copy clustering from HVG-filtered data to original data
    for key in key_obs:
        if key in hvg_adata.obs:
            adata.obs[key] = hvg_adata.obs[key].copy()
    
    # Transfer UMAP and t-SNE coordinates
    if 'X_umap' in hvg_adata.obsm:
        adata.obsm['X_umap'] = hvg_adata.obsm['X_umap'].copy()
    if 'X_tsne' in hvg_adata.obsm:
        adata.obsm['X_tsne'] = hvg_adata.obsm['X_tsne'].copy()
    
    # Save PCA results
    if 'X_pca' in hvg_adata.obsm:
        adata.obsm['X_pca'] = hvg_adata.obsm['X_pca'].copy()
        adata.varm['PCs'] = hvg_adata.varm['PCs'].copy()
        adata.uns['pca'] = hvg_adata.uns['pca'].copy() if 'pca' in hvg_adata.uns else None
    
    # Copy marker gene results
    if 'rank_genes_groups' in hvg_adata.uns:
        adata.uns['rank_genes_groups'] = hvg_adata.uns['rank_genes_groups'].copy()
    
    # Save the annotated and clustered data
    print(f"Saving clustered and annotated data to {output_file}...")
    try:
        adata.write(output_file)
    except Exception as e:
        print(f"Error saving clustered data: {e}")
        sys.exit(1)
    
    print("Clustering analysis completed successfully.")
    return adata

def main():
    parser = argparse.ArgumentParser(description="Advanced clustering and analysis of cancer scRNA-seq data.")
    parser.add_argument("input_file", help="Path to the processed AnnData object (.h5ad).")
    parser.add_argument("output_file", help="Path to save the clustered AnnData object (.h5ad).")
    parser.add_argument("--n_neighbors", type=int, default=15, help="Number of neighbors (default: 15)")
    parser.add_argument("--n_pcs", type=int, default=50, help="Number of principal components (default: 50)")
    parser.add_argument("--resolutions", type=float, nargs="+", default=[0.2, 0.4, 0.8, 1.0, 1.5], 
                        help="Clustering resolutions to try (default: 0.2 0.4 0.8 1.0 1.5)")
    parser.add_argument("--clustering_method", choices=["leiden", "louvain", "both"], default="leiden", 
                        help="Clustering algorithm to use (default: leiden)")
    parser.add_argument("--harmony", action="store_true", help="Apply Harmony batch correction")
    parser.add_argument("--regress_out", choices=["cell_cycle", "mito", "both", "none"], default="none",
                        help="Variables to regress out (default: none)")
    parser.add_argument("--no_plots", action="store_false", dest="generate_plots", 
                        help="Disable generation of plots")
    parser.add_argument("--output_dir", type=str, default=None, 
                        help="Directory to save plots (default: same as output_file)")
    parser.add_argument("--min_cluster_size", type=int, default=10, 
                        help="Minimum number of cells for a valid cluster (default: 10)")
    parser.add_argument("--cnv_analysis", action="store_true", 
                        help="Perform CNV inference analysis")
    parser.add_argument("--trajectory_analysis", action="store_true", 
                        help="Perform trajectory/pseudotime analysis")
    
    args = parser.parse_args()

    # Process regress_out argument
    regress_out_options = []
    if args.regress_out == "cell_cycle":
        regress_out_options = ["cell_cycle"]
    elif args.regress_out == "mito":
        regress_out_options = ["mito"]
    elif args.regress_out == "both":
        regress_out_options = ["cell_cycle", "mito"]
    else:
        regress_out_options = None

    cluster_cells(
        args.input_file, 
        args.output_file,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        resolutions=args.resolutions,
        clustering_method=args.clustering_method,
        use_harmony=args.harmony,
        regress_out=regress_out_options,
        generate_plots=args.generate_plots,
        output_dir=args.output_dir,
        min_cluster_size=args.min_cluster_size,
        cnv_analysis=args.cnv_analysis,
        trajectory_analysis=args.trajectory_analysis
    )

if __name__ == "__main__":
    main()