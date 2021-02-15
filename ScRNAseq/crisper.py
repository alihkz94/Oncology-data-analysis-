# Description: This script processes CRISPR UMI data to normalize UMI counts per cell, filter cells based on total UMI count, and identify enriched guide RNAs.
# Usage: python crispr.py
#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def main():
    # Example parameters (hard-coded)
    input_file = "input/crispr_umi/raw_umi_data.csv"  # Input CSV with columns: cell_barcode, guide_id, umi_count
    output_file = "output/crispr_umi/processed_umi_data.csv"  # Processed normalized UMI matrix
    plot_file = "output/crispr_umi/umi_distribution.pdf"  # UMI distribution plot
    min_total_umi = 500         # Minimum total UMI count per cell for filtering
    enrichment_threshold = 1.5  # Example threshold for guide enrichment (mean normalized count)

    # Create output directory if it does not exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    print("Loading CRISPR UMI data from:", input_file)
    try:
        # Assume CSV has columns: cell_barcode, guide_id, umi_count
        df = pd.read_csv(input_file)
    except Exception as e:
        print("Error loading input file:", e)
        return

    print("Aggregating UMI counts per cell and guide...")
    # Pivot table: rows = cell_barcode, columns = guide_id, aggregated sum of umi_count
    umi_matrix = df.pivot_table(index='cell_barcode', columns='guide_id', values='umi_count', aggfunc='sum', fill_value=0)
    
    print("Calculating total UMI counts per cell...")
    umi_matrix['total_umi'] = umi_matrix.sum(axis=1)
    
    print(f"Filtering cells with total UMI count less than {min_total_umi}...")
    filtered_matrix = umi_matrix[umi_matrix['total_umi'] >= min_total_umi].copy()
    total_umi = filtered_matrix.pop('total_umi')  # Remove total_umi column for normalization
    
    print("Normalizing UMI counts by total counts per cell (scaled to 10,000)...")
    normalized_matrix = filtered_matrix.div(total_umi, axis=0) * 1e4

    print("Calculating guide RNA enrichment (mean normalized count)...")
    guide_means = normalized_matrix.mean(axis=0)
    enriched_guides = guide_means[guide_means > enrichment_threshold].index.tolist()
    print("Enriched guides identified:", enriched_guides)

    print("Saving processed UMI data to:", output_file)
    normalized_matrix.to_csv(output_file)

    print("Generating UMI distribution plot...")
    plt.figure(figsize=(10, 6))
    plt.hist(total_umi, bins=50, color='skyblue', edgecolor='black')
    plt.xlabel("Total UMI Count per Cell")
    plt.ylabel("Number of Cells")
    plt.title("Distribution of Total UMI Counts per Cell")
    plt.savefig(plot_file)
    plt.close()

    print("CRISPR UMI analysis completed successfully.")

if __name__ == "__main__":
    main()