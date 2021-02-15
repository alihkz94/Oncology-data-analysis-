#!/usr/bin/env python

"""
CNV Analysis Pipeline using infercnvpy

This script performs somatic CNV detection using infercnvpy on a processed scRNA-seq dataset.
It assumes the AnnData object (in .h5ad format) contains a cell annotation column (e.g., 'cell_type')
with labels 'normal' for reference cells and 'tumor' for cells of interest. The script uses hardcoded
parameters for demonstration, including an expression cutoff and a smoothing window length for CNV inference.
Results include a CNV heatmap and an updated AnnData file with inferred CNV scores, saved to an output directory.
"""

import os
import sys
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt

def main():
    # Hardcoded parameters for CNV analysis
    input_file = "output/scanpy/scRNAseq_clustered.h5ad"   # Pre-processed AnnData with clustering and cell annotations
    output_dir = "output/cnv_clonal_evolution"              # Directory to save CNV results and plots
    cell_annotation = "cell_type"                           # Column in adata.obs indicating cell type labels
    reference_label = "normal"                              # Reference label in cell_annotation for non-malignant cells
    tumor_label = "tumor"                                   # Label for tumor cells (used implicitly by selecting reference)
    cutoff = 1.0                                          # Expression cutoff to filter genes for CNV inference
    window_length = 101                                   # Smoothing window length for CNV signal segmentation
    cluster_by = "leiden"                                 # Optional: use clustering labels to group similar cells
    
    # Create output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)
    
    print("Loading AnnData from:", input_file)
    try:
        adata = sc.read_h5ad(input_file)
    except Exception as e:
        print("Error loading AnnData file:", e)
        sys.exit(1)
    
    # Verify that the required cell annotation column exists
    if cell_annotation not in adata.obs.columns:
        print(f"Error: '{cell_annotation}' column not found in adata.obs. Please add cell annotations.")
        sys.exit(1)
    
    print("Performing CNV inference using infercnvpy...")
    # Run CNV inference: normal cells are used as a reference for CNV estimation in tumor cells.
    # The parameter 'scale_data' normalizes gene expression prior to inference.
    cnv.tl.infercnv(
        adata,
        reference_key=cell_annotation,
        reference_cat=reference_label,
        cutoff=cutoff,
        window_length=window_length,
        cluster_by=cluster_by,
        scale_data=True,
        verbose=True
    )
    
    print("Generating CNV heatmap plot...")
    # Generate and save a heatmap of inferred CNVs; the heatmap is saved as a PDF.
    heatmap_file = os.path.join(output_dir, "infercnv_heatmap.pdf")
    cnv.pl.infercnv(adata, output_format="pdf", save=True, show=False, output_file=heatmap_file)
    
    # Save the updated AnnData object with CNV annotations
    output_h5ad = os.path.join(output_dir, "scRNAseq_cnv_inferred.h5ad")
    adata.write(output_h5ad)
    
    print("CNV analysis completed successfully.")
    print("Results saved to:", output_dir)

if __name__ == "__main__":
    main()
