# Oncology Data Analysis

This repository contains bioinformatics scripts and pipelines for analyzing oncology-related genomic data, with a focus on single-cell RNA sequencing (scRNA-seq) and whole-exome sequencing (WES) analyses.

## Repository Structure

- **`ScRNAseq/`**: Single-cell RNA sequencing analysis scripts
  - `cluster_cells.py`: Advanced cell clustering with quality assessment using scanpy
  - `CNV.py`: Copy number variation analysis using infercnvpy
  - `crisper.py`: CRISPR UMI data processing, normalization and guide RNA enrichment analysis
  - `mast.R`: Seurat-based scRNA-seq analysis using the MAST differential expression framework
  - `ml_models.py`: Machine learning models for treatment response prediction using scVI for latent representation
  - `scanpy.py`: Standard scanpy-based scRNA-seq analysis pipeline with quality control
  - `seurat.R`: Comprehensive Seurat-based workflow for cancer scRNA-seq with cell cycle analysis
  - `spatial_analysis.py`: Spatial transcriptomics processing and visualization using Scanpy/Visium

- **`WES/`**: Whole-exome sequencing analysis scripts
  - `jwes.sh`: Shell script for WES data processing

## Key Features

- **Cancer genomics**: Tools optimized for tumor heterogeneity and copy number analysis 
- **Single-cell analytics**: Implementation of best practices for scRNA-seq quality control and clustering
- **Treatment response**: ML integration for predictive modeling from transcriptomic signatures
- **Spatial analysis**: Support for spatial transcriptomics in tumor microenvironment studies
- **Reproducible workflows**: Command-line interface with standardized parameters

## Technologies

- **Languages**: Python, R, Shell
- **Core libraries**: 
  - Python: scanpy, scVI, infercnvpy, scikit-learn, PyTorch
  - R: Seurat, MAST, dplyr, ggplot2