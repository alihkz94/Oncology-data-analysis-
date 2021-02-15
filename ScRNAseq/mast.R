#!/usr/bin/env Rscript

###############################################################################
# Seurat Pipeline for Cancer Single-Cell RNA-seq Analysis                     #
# Author: Ali Hakimzadeh                                                      #
# Date: March 2021                                                            #
###############################################################################

# ============================================================================ #
# SECTION 1: SETUP AND INITIALIZATION                                          #
# ============================================================================ #

# Load required packages for analysis
# Each package serves a specific purpose in our analysis workflow
message("Loading required packages...")
required_packages <- c(
  "Seurat",      # Core single-cell analysis framework
  "dplyr",       # Data manipulation
  "ggplot2",     # Data visualization
  "patchwork",   # Combining plots
  "optparse",    # Command-line argument parsing
  "RColorBrewer" # Color palettes for visualization
)

# Load each package, installing if necessary
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing package: ", pkg))
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Set default plotting theme for publication quality figures
theme_set(theme_bw(base_size = 12) + 
          theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                legend.position = "right",
                panel.grid.minor = element_blank()))

# Create a color palette for visualization
cluster_colors <- colorRampPalette(brewer.pal(9, "Set1"))(20)

# ============================================================================ #
# SECTION 2: PARSE COMMAND-LINE ARGUMENTS                                      #
# ============================================================================ #

# Define command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Path to the raw 10X data folder", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="seurat_results",
              help="Output directory for results (default: seurat_results)", metavar="character"),
  make_option(c("-p", "--project"), type="character", default="CancerSC",
              help="Project name (default: CancerSC)", metavar="character"),
  make_option(c("--min_features"), type="integer", default=200,
              help="Minimum number of genes per cell (default: 200)", metavar="integer"),
  make_option(c("--max_features"), type="integer", default=6000,
              help="Maximum number of genes per cell (default: 6000)", metavar="integer"),
  make_option(c("--min_cells"), type="integer", default=3,
              help="Minimum number of cells per gene (default: 3)", metavar="integer"),
  make_option(c("--mt_pattern"), type="character", default="^MT-",
              help="Pattern to identify mitochondrial genes (default: ^MT-)", metavar="character"),
  make_option(c("--mt_cutoff"), type="double", default=20,
              help="Maximum percentage of mitochondrial reads (default: 20%)", metavar="double"),
  make_option(c("--n_variable"), type="integer", default=2000,
              help="Number of variable features to select (default: 2000)", metavar="integer"),
  make_option(c("--n_dims"), type="integer", default=30,
              help="Number of PCA dimensions to compute and use (default: 30)", metavar="integer"),
  make_option(c("--resolution"), type="double", default=0.5,
              help="Clustering resolution parameter (default: 0.5)", metavar="double"),
  make_option(c("--seed"), type="integer", default=42,
              help="Random seed for reproducibility (default: 42)", metavar="integer")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list, 
                          description="Seurat pipeline for cancer scRNA-seq analysis")
opt <- parse_args(opt_parser)

# Check for required input
if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Error: Input data path is required.", call.=FALSE)
}

# ============================================================================ #
# SECTION 3: CREATE OUTPUT DIRECTORY AND SETUP LOGGING                         #
# ============================================================================ #

# Create output directory structure
message(paste0("Creating output directory structure in: ", opt$output))
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(opt$output, "figures"), showWarnings = FALSE)
dir.create(file.path(opt$output, "data"), showWarnings = FALSE)

# Set up logging
log_file <- file.path(opt$output, "analysis_log.txt")
log_conn <- file(log_file, open = "wt")
writeLines(paste("Seurat Analysis Log -", format(Sys.time())), log_conn)

# Function to log messages to both console and log file
log_message <- function(message) {
  cat(message, "\n")
  writeLines(message, log_conn)
}

# Set random seed for reproducibility
set.seed(opt$seed)

# Print analysis parameters
log_message("==== CANCER SINGLE-CELL RNA-SEQ ANALYSIS ====")
log_message(paste("Input data:", opt$input))
log_message(paste("Output directory:", opt$output))
log_message(paste("Project name:", opt$project))
log_message(paste("Minimum features per cell:", opt$min_features))
log_message(paste("Maximum features per cell:", opt$max_features))
log_message(paste("Minimum cells per gene:", opt$min_cells))
log_message(paste("Mitochondrial pattern:", opt$mt_pattern))
log_message(paste("Mitochondrial cutoff:", opt$mt_cutoff, "%"))
log_message(paste("Number of variable features:", opt$n_variable))
log_message(paste("Number of PCA dimensions:", opt$n_dims))
log_message(paste("Clustering resolution:", opt$resolution))
log_message(paste("Random seed:", opt$seed))
log_message("==========================================")

# ============================================================================ #
# SECTION 4: DATA LOADING AND QUALITY CONTROL                                  #
# ============================================================================ #

# Step 1: Load the 10X data
log_message("\n==== DATA LOADING ====")
log_message(paste("Loading 10X data from:", opt$input))

# Read the 10X data
tryCatch({
  # Determine if input is a directory or h5 file
  if (dir.exists(opt$input)) {
    seurat_data <- Read10X(data.dir = opt$input)
    log_message("Data read from 10X directory")
  } else if (grepl("\\.h5$", opt$input)) {
    seurat_data <- Read10X_h5(filename = opt$input)
    log_message("Data read from 10X H5 file")
  } else {
    stop("Input must be either a 10X directory or an H5 file")
  }
  
  log_message(paste("Data matrix dimensions:", dim(seurat_data)[1], "genes and", 
                   dim(seurat_data)[2], "cells"))
}, error = function(e) {
  log_message(paste("Error loading data:", e$message))
  stop("Failed to load data. Check the input path and file format.")
})

# Step 2: Create a Seurat object
log_message("\n==== CREATING SEURAT OBJECT ====")
log_message("Creating Seurat object with initial filtering...")

# Create Seurat object with basic filtering criteria
seurat_obj <- CreateSeuratObject(
  counts = seurat_data,
  project = opt$project,
  min.cells = opt$min_cells,     # Only keep genes detected in at least this many cells
  min.features = opt$min_features # Only keep cells with at least this many genes
)

# Step 3: Calculate Quality Control (QC) metrics
log_message("\n==== QUALITY CONTROL METRICS ====")
log_message("Calculating quality control metrics...")

# Calculate percentage of mitochondrial genes
# High mitochondrial content often indicates cell stress or dying cells
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = opt$mt_pattern)

# Calculate percentage of ribosomal genes
# This can help identify cells with high translational activity
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Print initial QC statistics
log_message("Initial dataset statistics:")
log_message(paste("  Total cells:", ncol(seurat_obj)))
log_message(paste("  Total genes:", nrow(seurat_obj)))
log_message(paste("  Median genes per cell:", median(seurat_obj$nFeature_RNA)))
log_message(paste("  Median UMIs per cell:", median(seurat_obj$nCount_RNA)))
log_message(paste("  Median % mitochondrial:", round(median(seurat_obj$percent.mt), 2)))

# Create QC plots for visualization
# These plots help identify appropriate thresholds for filtering
log_message("Creating QC visualization plots...")

# Violin plots of key QC metrics
qc_violin <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  ncol = 4,
  pt.size = 0.1
) + ggtitle("Quality Control Metrics")

ggsave(
  file.path(opt$output, "figures", "01_qc_violins.png"),
  plot = qc_violin,
  width = 12,
  height = 5,
  dpi = 300
)

# Create gene count vs UMI count scatter plot
scatter1 <- FeatureScatter(
  seurat_obj, 
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  pt.size = 0.5
) + ggtitle("Gene Count vs UMI Count")

# Create mitochondrial percentage vs UMI count scatter plot
scatter2 <- FeatureScatter(
  seurat_obj, 
  feature1 = "nCount_RNA",
  feature2 = "percent.mt",
  pt.size = 0.5
) + ggtitle("Mitochondrial % vs UMI Count")

# Combine scatter plots using patchwork
combined_scatter <- scatter1 + scatter2

ggsave(
  file.path(opt$output, "figures", "02_qc_scatter_plots.png"),
  plot = combined_scatter,
  width = 12,
  height = 5,
  dpi = 300
)

# Step 4: Filter cells based on QC metrics
log_message("\n==== FILTERING CELLS ====")
log_message(paste("Filtering cells based on QC thresholds:"))
log_message(paste("  Minimum genes per cell:", opt$min_features))
log_message(paste("  Maximum genes per cell:", opt$max_features))
log_message(paste("  Maximum mitochondrial %:", opt$mt_cutoff))

# Store pre-filtering metrics for reporting
pre_filter_cells <- ncol(seurat_obj)

# Apply QC thresholds to filter low-quality cells
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= opt$min_features &
           nFeature_RNA <= opt$max_features &
           percent.mt <= opt$mt_cutoff
)

# Report filtering results
log_message(paste("Cells after filtering:", ncol(seurat_obj), 
                 sprintf("(%.1f%% retained)", 100 * ncol(seurat_obj) / pre_filter_cells)))
log_message(paste("Genes after filtering:", nrow(seurat_obj)))

# ============================================================================ #
# SECTION 5: NORMALIZATION AND FEATURE SELECTION                               #
# ============================================================================ #

# Step 5: Normalize data
log_message("\n==== NORMALIZATION ====")
log_message("Normalizing data with LogNormalize method...")
log_message("This transforms feature counts to account for library size differences")

# LogNormalize: log(1 + x) transformation after scaling by library size
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Step 6: Identify highly variable features
log_message("\n==== FEATURE SELECTION ====")
log_message(paste("Identifying top", opt$n_variable, "variable features..."))
log_message("Variable features show high cell-to-cell variation and are used for dimensionality reduction")

# Find variable features using variance stabilizing transformation (VST)
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = opt$n_variable
)

# Report number of variable features
variable_features <- VariableFeatures(seurat_obj)
log_message(paste("Selected", length(variable_features), "variable features"))

# Plot variable features
var_plot <- VariableFeaturePlot(seurat_obj)

# Label top 10 most variable genes
top10_var_genes <- head(variable_features, 10)
var_plot_labeled <- LabelPoints(
  plot = var_plot, 
  points = top10_var_genes, 
  repel = TRUE
)

ggsave(
  file.path(opt$output, "figures", "03_variable_features.png"),
  plot = var_plot_labeled,
  width = 10,
  height = 8,
  dpi = 300
)

# Step 7: Scale data
log_message("\n==== DATA SCALING ====")
log_message("Scaling data (zero mean and unit variance)...")
log_message("This is necessary before dimensional reduction to give equal weight to all genes")

# Scale data to zero mean and unit variance
seurat_obj <- ScaleData(
  seurat_obj,
  features = rownames(seurat_obj),
  vars.to.regress = c("nCount_RNA", "percent.mt") # Regress out technical factors
)

# ============================================================================ #
# SECTION 6: DIMENSIONAL REDUCTION AND CLUSTERING                              #
# ============================================================================ #

# Step 8: Run PCA
log_message("\n==== PRINCIPAL COMPONENT ANALYSIS ====")
log_message(paste("Running PCA with", opt$n_dims, "dimensions..."))
log_message("PCA captures the primary sources of heterogeneity in the dataset")

# Perform PCA on scaled data using variable features
seurat_obj <- RunPCA(
  seurat_obj,
  features = variable_features,
  npcs = opt$n_dims,
  verbose = FALSE
)

# Create an elbow plot to determine significant PCs
pca_elbow <- ElbowPlot(
  seurat_obj, 
  ndims = opt$n_dims
) + ggtitle("PCA Elbow Plot")

ggsave(
  file.path(opt$output, "figures", "04_pca_elbow.png"),
  plot = pca_elbow,
  width = 8,
  height = 6,
  dpi = 300
)

# Step 9: Determine number of PCs to use for clustering
# For educational purposes, we'll use a heuristic method 
# that looks at where the elbow in the PC variance plot occurs
pct_variance <- seurat_obj[["pca"]]@stdev^2 / sum(seurat_obj[["pca"]]@stdev^2) * 100
cumulative_variance <- cumsum(pct_variance)

dims_use <- min(opt$n_dims, which(cumulative_variance > 80)[1])
if (is.na(dims_use) || dims_use < 10) dims_use <- min(20, opt$n_dims)

log_message(paste("Using", dims_use, "PCs for downstream analysis based on variance explained"))

# Step 10: Construct KNN graph
log_message("\n==== CLUSTERING ANALYSIS ====")
log_message("Building nearest neighbor graph...")
log_message("This graph captures relationships between cells in PC space")

# Find nearest neighbors
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:dims_use,
  k.param = 20
)

# Step 11: Perform clustering
log_message(paste("Clustering cells using Louvain algorithm at resolution", opt$resolution))
log_message("The resolution parameter controls clustering granularity (higher = more clusters)")

# Find clusters using the Louvain algorithm
seurat_obj <- FindClusters(
  seurat_obj,
  resolution = opt$resolution,
  algorithm = 1,
  random.seed = opt$seed
)

# Report cluster information
cluster_counts <- table(Idents(seurat_obj))
log_message(paste("Identified", length(cluster_counts), "clusters"))
log_message("Cells per cluster:")
for (i in 1:length(cluster_counts)) {
  log_message(paste("  Cluster", names(cluster_counts)[i], ":", cluster_counts[i], "cells"))
}

# Step 12: Run UMAP for visualization
log_message("\n==== VISUALIZATION WITH UMAP ====")
log_message("Running UMAP dimensional reduction...")
log_message("UMAP provides a 2D representation that preserves both local and global structure")

# Run UMAP algorithm
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:dims_use,
  seed.use = opt$seed
)

# Create UMAP plot colored by cluster
umap_clusters <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.5,
  label.size = 4
) + 
ggtitle("UMAP: Cell Clusters") +
scale_color_manual(values = cluster_colors)

ggsave(
  file.path(opt$output, "figures", "05_umap_clusters.png"),
  plot = umap_clusters,
  width = 10,
  height = 8,
  dpi = 300
)

# ============================================================================ #
# SECTION 7: BIOMARKER IDENTIFICATION AND CHARACTERIZATION                     #
# ============================================================================ #

# Step 13: Identify marker genes for each cluster
log_message("\n==== MARKER GENE IDENTIFICATION ====")
log_message("Finding differentially expressed genes for each cluster...")
log_message("These marker genes help define the biological identity of each cell type")

# Find all cluster marker genes using Wilcoxon Rank Sum test
all_markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,           # Only positive markers
  min.pct = 0.25,            # Expressed in at least 25% of cells in the cluster
  logfc.threshold = 0.25,    # At least 0.25 log2FC difference
  test.use = "wilcox",       # Statistical test
  return.thresh = 0.05       # Only return significant markers (p<0.05)
)

# Save marker genes to file
write.csv(
  all_markers,
  file = file.path(opt$output, "data", "cluster_marker_genes.csv"),
  row.names = FALSE
)

# Report number of markers found
log_message(paste("Found", nrow(all_markers), "marker genes across all clusters"))
log_message("Top markers per cluster saved to cluster_marker_genes.csv")

# Get top markers for each cluster
top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

log_message("Top markers by cluster:")
for (cluster_id in unique(top_markers$cluster)) {
  cluster_markers <- top_markers %>% filter(cluster == cluster_id)
  log_message(paste("  Cluster", cluster_id, ":", 
                   paste(cluster_markers$gene[1:min(3,nrow(cluster_markers))], collapse=", "), "..."))
}

# Step 14: Create a heatmap of top markers
log_message("\n==== MARKER VISUALIZATION ====")
log_message("Generating marker gene heatmap...")
log_message("The heatmap shows relative expression of top markers across clusters")

# Create a heatmap of top markers
marker_heatmap <- DoHeatmap(
  seurat_obj,
  features = top_markers$gene,
  group.by = "seurat_clusters",
  angle = 0
) + ggtitle("Top Marker Genes by Cluster")

ggsave(
  file.path(opt$output, "figures", "06_marker_heatmap.png"),
  plot = marker_heatmap,
  width = 12,
  height = 10,
  dpi = 300
)

# Step 15: Create feature plots for top markers
log_message("Generating feature plots for top markers...")
log_message("These plots show the expression distribution of marker genes on the UMAP")

# Plot top 2 markers for each of the first 4 clusters (for demonstration)
top_markers_to_plot <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) %>%
  filter(as.numeric(cluster) <= 4) %>%
  pull(gene)

# Generate feature plots
feature_plots <- FeaturePlot(
  seurat_obj,
  features = top_markers_to_plot,
  min.cutoff = "q10",
  max.cutoff = "q90",
  ncol = 4,
  order = TRUE
) & theme(plot.title = element_text(size = 10))

ggsave(
  file.path(opt$output, "figures", "07_marker_feature_plots.png"),
  plot = feature_plots,
  width = 16,
  height = 10,
  dpi = 300
)

# ============================================================================ #
# SECTION 8: CANCER-SPECIFIC ANALYSIS                                          #
# ============================================================================ #

# Step 16: Calculate cell cycle scores
# Cell cycle is often dysregulated in cancer cells
log_message("\n==== CANCER-SPECIFIC ANALYSES ====")
log_message("Calculating cell cycle scores...")
log_message("Cell cycle dysregulation is a hallmark of cancer")

# Calculate cell cycle scores
seurat_obj <- CellCycleScoring(
  seurat_obj,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

# Plot cell cycle phase
cell_cycle_plot <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "Phase",
  pt.size = 0.5
) + ggtitle("Cell Cycle Phase")

ggsave(
  file.path(opt$output, "figures", "08_cell_cycle.png"),
  plot = cell_cycle_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Step 17: Check for cancer-related gene signatures
# Here we examine common oncogenes and tumor suppressors
log_message("Checking expression of known cancer-related genes...")

# Define common cancer-related genes
oncogenes <- c("MYC", "EGFR", "KRAS", "ERBB2", "CCND1", "MDM2")
tumor_suppressors <- c("TP53", "PTEN", "RB1", "CDKN2A", "APC")

# Filter for genes present in the dataset
oncogenes_in_data <- oncogenes[oncogenes %in% rownames(seurat_obj)]
tumor_suppressors_in_data <- tumor_suppressors[tumor_suppressors %in% rownames(seurat_obj)]
cancer_genes <- c(oncogenes_in_data, tumor_suppressors_in_data)

log_message(paste("Found", length(cancer_genes), "cancer-related genes in the dataset"))

# If cancer genes are found, create visualization
if (length(cancer_genes) > 0) {
  # Visualize expression of cancer genes
  cancer_plot <- FeaturePlot(
    seurat_obj,
    features = cancer_genes[1:min(6, length(cancer_genes))],
    pt.size = 0.5,
    ncol = 3,
    order = TRUE
  )
  
  ggsave(
    file.path(opt$output, "figures", "09_cancer_genes.png"),
    plot = cancer_plot,
    width = 15,
    height = 10,
    dpi = 300
  )
  
  # Create violin plots of cancer genes by cluster
  cancer_violins <- VlnPlot(
    seurat_obj,
    features = cancer_genes[1:min(4, length(cancer_genes))],
    ncol = 2,
    pt.size = 0
  ) & theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    file.path(opt$output, "figures", "10_cancer_genes_violins.png"),
    plot = cancer_violins,
    width = 12,
    height = 10,
    dpi = 300
  )
}

# ============================================================================ #
# SECTION 9: SAVE RESULTS AND FINISH                                           #
# ============================================================================ #

# Step 18: Save the Seurat object for future analysis
log_message("\n==== SAVING RESULTS ====")
log_message(paste("Saving Seurat object to:", file.path(opt$output, "data", "seurat_object.rds")))
saveRDS(seurat_obj, file = file.path(opt$output, "data", "seurat_object.rds"))

# Save cluster information as a simple table
cluster_info <- data.frame(
  Cell = colnames(seurat_obj),
  Cluster = Idents(seurat_obj),
  UMAP_1 = seurat_obj@reductions$umap@cell.embeddings[,1],
  UMAP_2 = seurat_obj@reductions$umap@cell.embeddings[,2]
)

write.csv(
  cluster_info,
  file = file.path(opt$output, "data", "cell_clusters.csv"),
  row.names = FALSE
)

# Calculate and report runtime
end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
log_message("\n==== ANALYSIS COMPLETE ====")
log_message(paste("Analysis completed in", round(runtime, 2), "minutes"))
log_message(paste("Results saved in:", opt$output))
log_message(paste("Date and time:", format(end_time)))

# Close the log connection
close(log_conn)

# ============================================================================ #
# END OF SCRIPT                                                                #
# ============================================================================ #