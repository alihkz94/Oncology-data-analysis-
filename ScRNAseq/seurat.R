#!/usr/bin/env Rscript
#==============================================================================
# Seurat Pipeline for Cancer Single-Cell RNA-seq Analysis
# 
# This script performs end-to-end analysis of cancer single-cell RNA-seq data
# using the Seurat framework. 
#
# Author: Ali Hakimzadeh
# Date: Februrary 18, 2021
#==============================================================================

#------------------------------------------------------------------------------
# 1. Setup and Library Loading
#------------------------------------------------------------------------------
start_time <- Sys.time()
print(paste("Analysis started at:", format(start_time)))

# Load required libraries with proper error handling
required_packages <- c("Seurat", "dplyr", "ggplot2", "patchwork", "optparse", "RColorBrewer")
invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Package", pkg, "is required but not installed. Please install it using:\n",
               "install.packages('", pkg, "')", sep=""))
  }
  suppressMessages(library(pkg, character.only = TRUE))
}))

# Set default plotting theme
theme_set(theme_bw() + 
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                legend.position = "right"))

#------------------------------------------------------------------------------
# 2. Command-line Arguments
#------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the raw 10X data folder", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="seurat_analysis",
              help="Output directory for analysis results (default: 'seurat_analysis')", metavar="character"),
  make_option(c("--project"), type="character", default="CancerSC",
              help="Project name for labeling (default: 'CancerSC')", metavar="character"),
  make_option(c("--min_features"), type="integer", default=200,
              help="Minimum number of features per cell (default: 200)", metavar="integer"),
  make_option(c("--min_cells"), type="integer", default=3,
              help="Minimum number of cells expressing a gene (default: 3)", metavar="integer"),
  make_option(c("--mt_cutoff"), type="double", default=20,
              help="Maximum percentage of mitochondrial reads (default: 20%)", metavar="double"),
  make_option(c("--nfeatures"), type="integer", default=2000,
              help="Number of variable features to select (default: 2000)", metavar="integer"),
  make_option(c("--n_dims"), type="integer", default=30,
              help="Number of PCA dimensions to use (default: 30)", metavar="integer"),
  make_option(c("--resolution"), type="double", default=0.5,
              help="Resolution for clustering (default: 0.5)", metavar="double"),
  make_option(c("--skip_plots"), type="logical", default=FALSE, action="store_true",
              help="Skip generation of plots (default: FALSE)", metavar="logical")
)

opt_parser <- OptionParser(option_list=option_list, 
                          description="Comprehensive cancer scRNA-seq analysis pipeline using Seurat")
opt <- parse_args(opt_parser)

# Input validation
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Error: Input data folder must be provided.", call.=FALSE)
}

# Create output directory structure
output_dir <- opt$output
fig_dir <- file.path(output_dir, "figures")
data_dir <- file.path(output_dir, "data")
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(fig_dir, showWarnings=FALSE)
dir.create(data_dir, showWarnings=FALSE)

# Initialize log file
log_file <- file.path(output_dir, "analysis_log.txt")
log_con <- file(log_file, "w")
writeLines(paste("Seurat Analysis Log -", format(Sys.time())), log_con)
writeLines(paste("Input data:", opt$input), log_con)
writeLines(paste("Output directory:", output_dir), log_con)
writeLines(paste("Project name:", opt$project), log_con)
writeLines("", log_con)

# Function to log messages to console and log file
log_message <- function(message) {
  cat(message, "\n")
  writeLines(message, log_con)
}

log_message("==== ANALYSIS PARAMETERS ====")
log_message(paste("Minimum features per cell:", opt$min_features))
log_message(paste("Minimum cells per gene:", opt$min_cells))
log_message(paste("Mitochondrial cutoff:", opt$mt_cutoff, "%"))
log_message(paste("Number of variable features:", opt$nfeatures))
log_message(paste("Number of PCA dimensions:", opt$n_dims))
log_message(paste("Clustering resolution:", opt$resolution))
log_message("=============================\n")

#------------------------------------------------------------------------------
# 3. Data Loading and Initial QC
#------------------------------------------------------------------------------
log_message("\n==== DATA LOADING AND PROCESSING ====")
log_message("Step 1: Loading raw 10X data...")

# Read the 10X data with error handling
tryCatch({
  seurat_data <- Read10X(data.dir = opt$input)
  log_message(paste("  Data loaded successfully with", dim(seurat_data)[2], "cells and", 
                    dim(seurat_data)[1], "genes"))
}, error = function(e) {
  log_message(paste("Error loading data:", e$message))
  stop("Failed to load data. Check the input path and file format.")
})

# Create the Seurat object
log_message("Step 2: Creating Seurat object...")
seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                project = opt$project,
                                min.cells = opt$min_cells,
                                min.features = opt$min_features)

# Calculate QC metrics
log_message("Step 3: Calculating QC metrics...")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Summarize data before filtering
log_message("\n==== PRE-FILTERING SUMMARY ====")
log_message(paste("Total cells:", ncol(seurat_obj)))
log_message(paste("Total genes:", nrow(seurat_obj)))
log_message(paste("Median genes per cell:", median(seurat_obj$nFeature_RNA)))
log_message(paste("Median UMIs per cell:", median(seurat_obj$nCount_RNA)))
log_message(paste("Median % mitochondrial:", median(seurat_obj$percent.mt)))

# Save QC plots if requested
if(!opt$skip_plots) {
  log_message("Creating QC plots...")
  
  # QC violin plots
  p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 3, pt.size = 0.1) + 
        ggtitle("Quality Control Metrics Before Filtering")
  
  # Feature scatter plots
  p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
        ggtitle("Gene Count vs UMI Count")
  
  p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
        ggtitle("Mitochondrial % vs UMI Count")
  
  # Save plots
  ggsave(file.path(fig_dir, "01_qc_violins_pre_filtering.png"), plot = p1, width = 12, height = 5)
  ggsave(file.path(fig_dir, "02_feature_scatter1_pre_filtering.png"), plot = p2, width = 8, height = 6)
  ggsave(file.path(fig_dir, "03_feature_scatter2_pre_filtering.png"), plot = p3, width = 8, height = 6)
}

#------------------------------------------------------------------------------
# 4. Cell Filtering 
#------------------------------------------------------------------------------
log_message("\n==== FILTERING CELLS ====")
log_message(paste("Filtering cells: nFeature_RNA >", opt$min_features, 
                 "and percent.mt <", opt$mt_cutoff))

# Store metrics before filtering
pre_filter_metrics <- data.frame(
  cells = ncol(seurat_obj),
  genes = nrow(seurat_obj),
  median_genes = median(seurat_obj$nFeature_RNA),
  median_umi = median(seurat_obj$nCount_RNA),
  median_mt = median(seurat_obj$percent.mt)
)

# Apply filters
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > opt$min_features & 
                                          percent.mt < opt$mt_cutoff)

# Store metrics after filtering
post_filter_metrics <- data.frame(
  cells = ncol(seurat_obj),
  genes = nrow(seurat_obj),
  median_genes = median(seurat_obj$nFeature_RNA),
  median_umi = median(seurat_obj$nCount_RNA),
  median_mt = median(seurat_obj$percent.mt)
)

# Summarize filtering effects
log_message("\n==== POST-FILTERING SUMMARY ====")
log_message(paste("Cells retained:", ncol(seurat_obj), 
                 sprintf("(%.1f%%)", 100 * ncol(seurat_obj) / pre_filter_metrics$cells)))
log_message(paste("Genes retained:", nrow(seurat_obj)))
log_message(paste("Median genes per cell:", median(seurat_obj$nFeature_RNA)))
log_message(paste("Median UMIs per cell:", median(seurat_obj$nCount_RNA)))
log_message(paste("Median % mitochondrial:", median(seurat_obj$percent.mt)))

# Plot post-filtering QC if requested
if(!opt$skip_plots) {
  log_message("Creating post-filtering QC plots...")
  
  p4 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
               ncol = 3, pt.size = 0.1) + 
        ggtitle("Quality Control Metrics After Filtering")
  
  ggsave(file.path(fig_dir, "04_qc_violins_post_filtering.png"), plot = p4, width = 12, height = 5)
}

#------------------------------------------------------------------------------
# 5. Normalization and Feature Selection
#------------------------------------------------------------------------------
log_message("\n==== NORMALIZATION AND FEATURE SELECTION ====")
log_message("Step 1: Normalizing data (LogNormalize method)...")
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

log_message(paste("Step 2: Identifying", opt$nfeatures, "highly variable features..."))
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", 
                                  nfeatures = opt$nfeatures)

# Display top variable features
top_features <- head(VariableFeatures(seurat_obj), 20)
log_message(paste("Top variable features:", paste(top_features, collapse=", ")))

if(!opt$skip_plots) {
  log_message("Creating variable feature plot...")
  
  # Variable features plot
  p5 <- VariableFeaturePlot(seurat_obj) +
        ggtitle("Highly Variable Genes")
  
  # Top variable features with labels
  p6 <- LabelPoints(plot = p5, points = top_features, repel = TRUE, xnudge = 0, ynudge = 0)
  
  ggsave(file.path(fig_dir, "05_variable_features.png"), plot = p5, width = 8, height = 6)
  ggsave(file.path(fig_dir, "06_variable_features_labeled.png"), plot = p6, width = 8, height = 6)
}

#------------------------------------------------------------------------------
# 6. Scaling and PCA
#------------------------------------------------------------------------------
log_message("\n==== DATA SCALING AND DIMENSIONALITY REDUCTION ====")
log_message("Step 1: Scaling data...")

# Scale using variable features plus specific genes of interest for cancer
cancer_genes <- c("MKI67", "PCNA", "TOP2A", "EPCAM", "KRT19", "EGFR", "MYC", "CCND1")
features_to_scale <- union(VariableFeatures(seurat_obj), 
                          intersect(cancer_genes, rownames(seurat_obj)))

log_message(paste("Scaling", length(features_to_scale), "features..."))
seurat_obj <- ScaleData(seurat_obj, features = features_to_scale)

log_message(paste("Step 2: Running PCA with", opt$n_dims, "dimensions..."))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), 
                     npcs = opt$n_dims, verbose = FALSE)

# Print top genes associated with PCs
log_message("\n==== TOP GENES FOR PRINCIPAL COMPONENTS ====")
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 10)
pc_loadings <- as.data.frame(seurat_obj[["pca"]]@feature.loadings)
write.csv(pc_loadings, file.path(data_dir, "pca_feature_loadings.csv"))

if(!opt$skip_plots) {
  log_message("Creating PCA plots...")
  
  # PCA plots
  p7 <- DimPlot(seurat_obj, reduction = "pca") + ggtitle("PCA Visualization")
  p8 <- ElbowPlot(seurat_obj, ndims = opt$n_dims) + ggtitle("PCA Elbow Plot")
  
  ggsave(file.path(fig_dir, "07_pca_plot.png"), plot = p7, width = 8, height = 7)
  ggsave(file.path(fig_dir, "08_pca_elbow.png"), plot = p8, width = 8, height = 6)
}

#------------------------------------------------------------------------------
# 7. Clustering
#------------------------------------------------------------------------------
log_message("\n==== CLUSTERING ====")
log_message(paste("Step 1: Finding nearest neighbors using first", opt$n_dims, "PCs..."))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:opt$n_dims)

log_message(paste("Step 2: Finding clusters with resolution =", opt$resolution, "..."))
seurat_obj <- FindClusters(seurat_obj, resolution = opt$resolution)

# Summarize clusters
cluster_counts <- table(Idents(seurat_obj))
log_message("\n==== CLUSTER SUMMARY ====")
log_message(paste("Total clusters identified:", length(cluster_counts)))
for (i in names(cluster_counts)) {
  log_message(paste("Cluster", i, ":", cluster_counts[i], "cells", 
                   sprintf("(%.1f%%)", 100 * cluster_counts[i] / sum(cluster_counts))))
}

#------------------------------------------------------------------------------
# 8. Dimensional Reduction for Visualization
#------------------------------------------------------------------------------
log_message("\n==== DIMENSIONAL REDUCTION FOR VISUALIZATION ====")
log_message("Step 1: Running UMAP...")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:opt$n_dims)

log_message("Step 2: Running t-SNE...")
seurat_obj <- RunTSNE(seurat_obj, dims = 1:opt$n_dims)

if(!opt$skip_plots) {
  log_message("Creating dimensional reduction visualizations...")
  
  # UMAP plot
  p9 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE) + 
       ggtitle("UMAP Clustering Visualization") + 
       theme(legend.position = "right")
  
  # t-SNE plot
  p10 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE) + 
        ggtitle("t-SNE Clustering Visualization") + 
        theme(legend.position = "right")
  
  ggsave(file.path(fig_dir, "09_umap_clusters.png"), plot = p9, width = 10, height = 8)
  ggsave(file.path(fig_dir, "10_tsne_clusters.png"), plot = p10, width = 10, height = 8)
}

#------------------------------------------------------------------------------
# 9. Finding Marker Genes
#------------------------------------------------------------------------------
log_message("\n==== IDENTIFYING MARKER GENES ====")
log_message("Finding markers for each cluster...")

# Find markers for all clusters
all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, 
                             logfc.threshold = 0.25, test.use = "wilcox")

log_message(paste("Total differentially expressed markers found:", nrow(all_markers)))
write.csv(all_markers, file.path(data_dir, "all_cluster_markers.csv"), row.names = FALSE)

# Extract top markers per cluster
top_n_markers <- 5
top_markers_by_cluster <- all_markers %>%
    group_by(cluster) %>%
    top_n(n = top_n_markers, wt = avg_log2FC)

log_message("\n==== TOP MARKER GENES BY CLUSTER ====")
for (i in unique(top_markers_by_cluster$cluster)) {
  cluster_markers <- top_markers_by_cluster %>% filter(cluster == i)
  log_message(paste("Cluster", i, "top markers:", 
                   paste(cluster_markers$gene, collapse=", ")))
}

if(!opt$skip_plots) {
  log_message("Creating marker gene visualizations...")
  
  # Visualize top markers for each cluster
  for (i in unique(all_markers$cluster)) {
    markers <- all_markers %>% 
      filter(cluster == i) %>% 
      top_n(n = top_n_markers, wt = avg_log2FC)
    
    if(nrow(markers) > 0) {
      # Feature plot for top genes
      genes_to_plot <- head(markers$gene, 4)  # Plot top 4 genes
      if(length(genes_to_plot) > 0) {
        p_feature <- FeaturePlot(seurat_obj, features = genes_to_plot, ncol = 2)
        ggsave(file.path(fig_dir, paste0("11_cluster", i, "_top_markers_featureplot.png")), 
               plot = p_feature, width = 10, height = 8)
      }
    }
  }
  
  # Create overall heatmap of top markers
  top_n_overall <- 10
  markers_to_plot <- all_markers %>%
      group_by(cluster) %>%
      top_n(n = top_n_overall, wt = avg_log2FC) %>%
      pull(gene) %>%
      unique()
  
  if(length(markers_to_plot) > 0) {
    p_heatmap <- DoHeatmap(seurat_obj, features = markers_to_plot, size = 3, angle = 0) + 
                 ggtitle("Top Marker Genes by Cluster")
    ggsave(file.path(fig_dir, "12_marker_genes_heatmap.png"), plot = p_heatmap, 
           width = 14, height = 10)
  }
  
  # Create dot plot
  p_dotplot <- DotPlot(seurat_obj, features = head(markers_to_plot, 25), dot.scale = 8) + 
               ggtitle("Expression of Top Marker Genes") +
               theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(fig_dir, "13_marker_genes_dotplot.png"), plot = p_dotplot, 
         width = 12, height = 8)
}

#------------------------------------------------------------------------------
# 10. Cell Cycle Scoring (Important for Cancer Analysis)
#------------------------------------------------------------------------------
log_message("\n==== CELL CYCLE ANALYSIS ====")
log_message("Performing cell cycle scoring...")

# Cell cycle gene sets
s_genes <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", 
             "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", 
             "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", 
             "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", 
             "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")
g2m_genes <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", 
               "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", 
               "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", 
               "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "CDC20", "TTK", 
               "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", 
               "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", 
               "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")

# Filter gene sets to genes in the dataset
s_genes <- intersect(s_genes, rownames(seurat_obj))
g2m_genes <- intersect(g2m_genes, rownames(seurat_obj))

if(length(s_genes) > 0 && length(g2m_genes) > 0) {
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s_genes, g2m.features = g2m_genes)
  
  # Log cell cycle distribution
  cc_distribution <- table(seurat_obj$Phase)
  log_message(paste("G1 phase cells:", cc_distribution["G1"], 
                   sprintf("(%.1f%%)", 100 * cc_distribution["G1"] / sum(cc_distribution))))
  log_message(paste("S phase cells:", cc_distribution["S"], 
                   sprintf("(%.1f%%)", 100 * cc_distribution["S"] / sum(cc_distribution))))
  log_message(paste("G2M phase cells:", cc_distribution["G2M"], 
                   sprintf("(%.1f%%)", 100 * cc_distribution["G2M"] / sum(cc_distribution))))
  
  if(!opt$skip_plots) {
    log_message("Creating cell cycle plots...")
    
    p_cc1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "Phase") +
            ggtitle("Cell Cycle Phases")
    
    p_cc2 <- FeaturePlot(seurat_obj, features = c("S.Score", "G2M.Score"), ncol = 2) +
            ggtitle("Cell Cycle Scores")
    
    ggsave(file.path(fig_dir, "14_cell_cycle_phases.png"), plot = p_cc1, width = 10, height = 8)
    ggsave(file.path(fig_dir, "15_cell_cycle_scores.png"), plot = p_cc2, width = 12, height = 6)
  }
} else {
  log_message("Warning: Not enough cell cycle genes found in dataset. Skipping cell cycle scoring.")
}

#------------------------------------------------------------------------------
# 11. Cancer-Specific Feature Visualization
#------------------------------------------------------------------------------
log_message("\n==== CANCER-SPECIFIC GENE EXPRESSION ====")

# Define cancer-related gene sets
proliferation_genes <- c("MKI67", "PCNA", "TOP2A", "AURKA", "AURKB", "BIRC5")
epithelial_genes <- c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1")
mesenchymal_genes <- c("VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1")
stemness_genes <- c("CD44", "PROM1", "ALDH1A1", "ALDH1A3", "ITGA6", "SOX2", "NANOG", "POU5F1")
oncogenes <- c("MYC", "KRAS", "EGFR", "PIK3CA", "AKT1", "CCND1")
tumor_suppressors <- c("TP53", "PTEN", "RB1", "CDKN2A", "CDKN1A", "CDKN1B")

# Filter for genes present in the dataset
proliferation_genes <- intersect(proliferation_genes, rownames(seurat_obj))
epithelial_genes <- intersect(epithelial_genes, rownames(seurat_obj))
mesenchymal_genes <- intersect(mesenchymal_genes, rownames(seurat_obj))
stemness_genes <- intersect(stemness_genes, rownames(seurat_obj))
oncogenes <- intersect(oncogenes, rownames(seurat_obj))
tumor_suppressors <- intersect(tumor_suppressors, rownames(seurat_obj))

log_message("Visualizing cancer-relevant gene sets...")
cancer_feature_sets <- list(
  proliferation = proliferation_genes,
  epithelial = epithelial_genes,
  mesenchymal = mesenchymal_genes,
  stemness = stemness_genes,
  oncogenes = oncogenes,
  tumor_suppressors = tumor_suppressors
)

if(!opt$skip_plots) {
  # Create feature plots for each gene set
  for (set_name in names(cancer_feature_sets)) {
    genes <- cancer_feature_sets[[set_name]]
    if(length(genes) > 0) {
      log_message(paste("Plotting", set_name, "genes:", paste(genes, collapse=", ")))
      
      # Calculate module score for gene set
      seurat_obj <- AddModuleScore(seurat_obj, features = list(genes), 
                                  name = paste0(set_name, "_score"))
      
      # Feature plots - limited to 6 genes per plot
      if(length(genes) > 0) {
        for(i in seq(1, length(genes), 6)) {
          end_idx <- min(i+5, length(genes))
          genes_subset <- genes[i:end_idx]
          
          p_genes <- FeaturePlot(seurat_obj, features = genes_subset, 
                                ncol = min(3, length(genes_subset)))
          ggsave(file.path(fig_dir, paste0("16_", set_name, "_genes_", i, ".png")), 
                 plot = p_genes, width = 12, height = 8)
        }
      }
      
      # Module score on UMAP
      p_score <- FeaturePlot(seurat_obj, 
                           features = paste0(set_name, "_score1"), 
                           cols = c("lightgrey", "red")) + 
                ggtitle(paste(set_name, "score"))
      ggsave(file.path(fig_dir, paste0("17_", set_name, "_score.png")), 
             plot = p_score, width = 8, height = 7)
    }
  }
}

#------------------------------------------------------------------------------
# 12. Save Results
#------------------------------------------------------------------------------
log_message("\n==== SAVING RESULTS ====")

# Save Seurat object
log_message(paste("Saving Seurat object to", file.path(data_dir, "seurat_object.rds")))
saveRDS(seurat_obj, file = file.path(data_dir, "seurat_object.rds"))

# Save final cell metadata
log_message("Saving cell metadata...")
cell_metadata <- seurat_obj@meta.data
write.csv(cell_metadata, file = file.path(data_dir, "cell_metadata.csv"))

# Save cluster marker summary
log_message("Saving cluster marker summary...")
marker_summary <- all_markers %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
write.csv(marker_summary, file = file.path(data_dir, "top10_markers_by_cluster.csv"), row.names = FALSE)

#------------------------------------------------------------------------------
# 13. Final Summary
#------------------------------------------------------------------------------
end_time <- Sys.time()
time_taken <- round(difftime(end_time, start_time, units="mins"), 2)

log_message("\n==== ANALYSIS SUMMARY ====")
log_message(paste("Analysis completed in", time_taken, "minutes"))
log_message(paste("Total cells analyzed:", ncol(seurat_obj)))
log_message(paste("Total genes analyzed:", nrow(seurat_obj)))
log_message(paste("Number of clusters identified:", length(unique(Idents(seurat_obj)))))
log_message(paste("Results saved to:", output_dir))

if(!opt$skip_plots) {
  log_message(paste("Plots saved to:", fig_dir))
}

log_message("\nThank you for using the Seurat scRNA-seq analysis pipeline!")

# Close log file
close(log_con)