
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("celldex")
BiocManager::install("SingleR")
BiocManager::install("SingleCellExperiment")

library(Seurat)
library(Matrix)
library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(ggplot2)
library(patchwork)
library(SingleCellExperiment)
library(SingleR)
library(celldex)


options(Seurat.object.assay.version = "v3")

# -------------------------
# Input files
# -------------------------
pbmc_files <- sprintf(
  "./DATA/multi_config_%s_sample_filtered_feature_bc_matrix.h5",
  1:5
)

names(pbmc_files) <- paste0("multi_config_", 1:5)

demux_files <- sprintf(
  "./DEMUX/multi_config_%s_demux/assignments_refined.tsv.gz",
  1:5
)

names(demux_files) <- paste0("multi_config_", 1:5)

print(pbmc_files)
print(demux_files)

dir.create("PBMC_figures", showWarnings = FALSE)

# -------------------------
# Read one PBMC H5
# -------------------------
read_pbmc_h5 <- function(h5_file, sample_name) {
  
  x <- Read10X_h5(h5_file)
  
  print(sample_name)
  print(names(x))
  
  obj <- CreateSeuratObject(
    counts = x[["Gene Expression"]],
    project = sample_name
  )
  
  obj$orig.ident <- sample_name
  obj$Lane <- str_replace(sample_name, "multi_config_", "")
  
  if ("Antibody Capture" %in% names(x)) {
    obj[["ADT"]] <- CreateAssayObject(
      counts = x[["Antibody Capture"]]
    )
  }
  
  return(obj)
}

# -------------------------
# Read and merge 5 lanes
# -------------------------
pbmc_list <- mapply(
  FUN = read_pbmc_h5,
  h5_file = pbmc_files,
  sample_name = names(pbmc_files),
  SIMPLIFY = FALSE
)

pbmc_merged <- merge(
  pbmc_list[[1]],
  y = pbmc_list[2:5],
  add.cell.ids = names(pbmc_list)
)

print(table(pbmc_merged$orig.ident))
print(table(pbmc_merged$Lane))

# -------------------------
# Load demux assignments
# -------------------------
read_demux_one <- function(path, sample_name) {
  
  df <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
  
  barcode_col <- intersect(c("BARCODE", "barcode", "Barcode", "cell", "CELL"), names(df))[1]
  
  if (is.na(barcode_col)) {
    stop(paste("No barcode column found in:", path))
  }
  
  subj_col <- setdiff(names(df), barcode_col)[1]
  
  if (is.na(subj_col)) {
    stop(paste("No subject column found in:", path))
  }
  
  df %>%
    transmute(
      core_bc = sub("-.*$", "", .data[[barcode_col]]),
      subject = .data[[subj_col]],
      cell = paste0(sample_name, "_", core_bc, "-1"),
      demux_lane = sample_name
    ) %>%
    distinct(cell, .keep_all = TRUE)
}

demux_tbl <- purrr::map2_dfr(
  demux_files,
  names(demux_files),
  read_demux_one
)

print(head(demux_tbl))
print(table(demux_tbl$demux_lane))
print(table(demux_tbl$subject, useNA = "ifany"))

# -------------------------
# Add demux metadata to Seurat object
# -------------------------
meta <- pbmc_merged@meta.data %>%
  mutate(cell = rownames(.))

meta2 <- meta %>%
  left_join(demux_tbl, by = "cell") %>%
  mutate(
    SNP_call = case_when(
      is.na(subject) ~ "Unassigned",
      str_detect(as.character(subject), "\\+") ~ "Doublet",
      TRUE ~ "Singleton"
    )
  )

meta2 <- meta2[match(colnames(pbmc_merged), meta2$cell), , drop = FALSE]

pbmc_merged$subject <- meta2$subject
pbmc_merged$demux_lane <- meta2$demux_lane
pbmc_merged$SNP_call <- meta2$SNP_call

print(table(pbmc_merged$orig.ident, pbmc_merged$SNP_call, useNA = "ifany"))
print(table(pbmc_merged$subject, useNA = "ifany"))

# -------------------------
# QC metrics
# -------------------------
pbmc_merged <- PercentageFeatureSet(pbmc_merged, "^MT-", col.name = "percent_mito")
pbmc_merged <- PercentageFeatureSet(pbmc_merged, "^RP[SL]", col.name = "percent_ribo")

hb_genes <- grep("^HB", rownames(pbmc_merged), value = TRUE)
hb_genes <- hb_genes[!grepl("^HBP", hb_genes)]

pbmc_merged <- PercentageFeatureSet(
  pbmc_merged,
  features = hb_genes,
  col.name = "percent_hb"
)

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

p_qc_before <- VlnPlot(
  pbmc_merged,
  features = feats,
  group.by = "orig.ident",
  pt.size = 0,
  ncol = 3
) + NoLegend() + ggtitle("PBMC RNA QC Before Filtering")

# -------------------------
# ADT QC if present
# -------------------------
if ("ADT" %in% names(pbmc_merged@assays)) {
  
  pbmc_merged$nCount_ADT <- Matrix::colSums(
   # GetAssayData(pbmc_merged, assay = "ADT", slot = "counts")
    GetAssayData(pbmc_merged, assay = "ADT", layer = "counts")
    )
  
  pbmc_merged$nFeature_ADT <- Matrix::colSums(
  #  GetAssayData(pbmc_merged, assay = "ADT", slot = "counts") > 0
    GetAssayData(pbmc_merged, assay = "ADT", layer = "counts") > 0
    )
  
  print(table(pbmc_merged$orig.ident, pbmc_merged$nCount_ADT > 0))
  print(table(pbmc_merged$orig.ident, pbmc_merged$nCount_ADT > 10))
}

# -------------------------
# Filter cells
# Keep singlets only + RNA QC
# -------------------------
pbmc_filt <- subset(
  pbmc_merged,
  subset = SNP_call == "Singleton" &
    nFeature_RNA > 200 &
    percent_mito < 10
)

print(dim(pbmc_filt))
print(table(pbmc_filt$orig.ident))
print(table(pbmc_filt$subject, useNA = "ifany"))

p_qc_after <- VlnPlot(
  pbmc_filt,
  features = feats,
  group.by = "orig.ident",
  pt.size = 0,
  ncol = 3
) + NoLegend() + ggtitle("PBMC RNA QC After Filtering")



rm(pbmc_merged)
rm(pbmc_list)

# -------------------------
# Normalize ADT if present
# -------------------------
if ("ADT" %in% names(pbmc_filt@assays)) {
  DefaultAssay(pbmc_filt) <- "ADT"
  
  pbmc_filt <- NormalizeData(
    pbmc_filt,
    normalization.method = "CLR",
    margin = 2
  )
  
  DefaultAssay(pbmc_filt) <- "RNA"
}

# -------------------------
# RNA processing
# -------------------------
DefaultAssay(pbmc_filt) <- "RNA"

pbmc_filt <- NormalizeData(pbmc_filt)
pbmc_filt <- FindVariableFeatures(pbmc_filt)
pbmc_filt <- ScaleData(pbmc_filt)
pbmc_filt <- RunPCA(pbmc_filt, verbose = FALSE)
pbmc_filt <- FindNeighbors(pbmc_filt, dims = 1:30)
pbmc_filt <- FindClusters(pbmc_filt, resolution = 0.8)
pbmc_filt <- RunUMAP(pbmc_filt, dims = 1:30)

p_umap_cluster <- DimPlot(
  pbmc_filt,
  label = TRUE
) + NoLegend() + ggtitle("PBMC RNA UMAP Clusters")

p_umap_lane <- DimPlot(
  pbmc_filt,
  group.by = "orig.ident"
) + ggtitle("PBMC lanes")

p_umap_subject <- DimPlot(
  pbmc_filt,
  group.by = "subject"
) + ggtitle("PBMC subjects")

# -------------------------
# Monaco SingleR annotation
# -------------------------
DefaultAssay(pbmc_filt) <- "RNA"

monaco.ref <- celldex::MonacoImmuneData()
pbmc_sce <- as.SingleCellExperiment(pbmc_filt)

pbmc_monaco_main <- SingleR(
  test = pbmc_sce,
  assay.type.test = 1,
  ref = monaco.ref,
  labels = monaco.ref$label.main
)

pbmc_monaco_fine <- SingleR(
  test = pbmc_sce,
  assay.type.test = 1,
  ref = monaco.ref,
  labels = monaco.ref$label.fine
)

pbmc_filt$monaco_main <- pbmc_monaco_main$pruned.labels
pbmc_filt$monaco_fine <- pbmc_monaco_fine$pruned.labels

print(table(pbmc_filt$monaco_main, useNA = "ifany"))
print(table(pbmc_filt$monaco_fine, useNA = "ifany"))

p_monaco_main <- DimPlot(
  pbmc_filt,
  reduction = "umap",
  group.by = "monaco_main",
  label = TRUE,
  repel = TRUE
) + ggtitle("PBMC Monaco main annotation")

p_monaco_fine <- DimPlot(
  pbmc_filt,
  reduction = "umap",
  group.by = "monaco_fine",
  label = TRUE,
  repel = TRUE
) + ggtitle("PBMC Monaco fine annotation")

# -------------------------
# Save figures
# -------------------------
ggsave(
  "PBMC_figures/PBMC_QC_before_filtering.png",
  p_qc_before,
  width = 14,
  height = 7,
  dpi = 300
)

ggsave(
  "PBMC_figures/PBMC_QC_after_filtering.png",
  p_qc_after,
  width = 14,
  height = 7,
  dpi = 300
)

umap_basic <- (p_umap_cluster | p_umap_lane) / p_umap_subject

ggsave(
  "PBMC_figures/PBMC_UMAP_clusters_lanes_subjects.png",
  umap_basic,
  width = 14,
  height = 10,
  dpi = 600
)

annotation_fig <- p_monaco_main + p_monaco_fine

ggsave(
  "PBMC_figures/PBMC_Monaco_annotation_UMAP.png",
  annotation_fig,
  width = 14,
  height = 6,
  dpi = 600
)

ggsave(
  "PBMC_figures/PBMC_Monaco_annotation_UMAP.pdf",
  annotation_fig,
  width = 14,
  height = 6
)

###############################

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
library(tidyr)

# -------------------------
# Choose object
# -------------------------
obj <- pbmc_filt   # or use data.filt for Ag sorted cells

# Must have:
# obj$subject
# obj$orig.ident or condition
# obj$monaco_fine

# -------------------------
# Rename / simplify B cell labels
# -------------------------
obj$bcell_type <- obj$monaco_fine

bcell_keep <- c(
  "Exhausted B cells",
  "Naive B cells",
  "Non-switched memory B cells",
  "Switched memory B cells"
)

bcell_obj <- subset(obj, subset = bcell_type %in% bcell_keep)

# Optional cleaner names
bcell_obj$bcell_type <- gsub("Exhausted B cells", "Exhausted B cell", bcell_obj$bcell_type)
bcell_obj$bcell_type <- gsub("Naive B cells", "Naive", bcell_obj$bcell_type)
bcell_obj$bcell_type <- gsub("Non-switched memory B cells", "Non-switched Mem", bcell_obj$bcell_type)
bcell_obj$bcell_type <- gsub("Switched memory B cells", "Switched Mem", bcell_obj$bcell_type)

# -------------------------
# Panel C-style stacked barplot
# Percent B cell types per sample/subject
# -------------------------
freq_df <- bcell_obj@meta.data %>%
  dplyr::count(subject, orig.ident, bcell_type) %>%
  dplyr::group_by(subject, orig.ident) %>%
  dplyr::mutate(percent = n / sum(n) * 100) %>%
  ungroup()

freq_df$sample_label <- paste(freq_df$subject, freq_df$orig.ident, sep = "_")

pC <- ggplot(freq_df, aes(x = sample_label, y = percent, fill = bcell_type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  ylab("Percentage of cells") +
  xlab("") +
  ggtitle("B cell population frequencies") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_blank()
  )

pC

ggsave(
  "PBMC_figures/PanelC_Bcell_population_frequencies.png",
  pC,
  width = 8,
  height = 5,
  dpi = 600
)

# -------------------------
# Panel D-style heatmap
# Marker expression across B cell populations
# -------------------------
DefaultAssay(bcell_obj) <- "RNA"

markers <- c(
  "CD19", "MS4A1", "CD79A", "CD79B",
  "IGHD", "IGHM", "IGHG1", "IGHG2", "IGHG3", "IGHA1",
  "CD27", "TNFRSF13B", "BANK1",
  "MKI67", "TBX21", "FCRL5", "ITGAX",
  "PRDM1", "XBP1", "MZB1"
)

markers <- markers[markers %in% rownames(bcell_obj)]

# Set order like paper panel
bcell_obj$bcell_type <- factor(
  bcell_obj$bcell_type,
  levels = c(
    "Naive",
    "Non-switched Mem",
    "Switched Mem",
    "Exhausted B cell"
  )
)

Idents(bcell_obj) <- "bcell_type"

pD <- DoHeatmap(
  bcell_obj,
  features = markers,
  group.by = "bcell_type",
  size = 3
) +
  ggtitle("B cell marker expression")

pD

ggsave(
  "PBMC_figures/PanelD_Bcell_marker_heatmap.png",
  pD,
  width = 10,
  height = 6,
  dpi = 600
)

ggsave(
  "PBMC_figures/PanelD_Bcell_marker_heatmap.pdf",
  pD,
  width = 10,
  height = 6
)

# -------------------------
# Save combined C + D figure
# -------------------------
panel_CD <- pC / pD + plot_layout(heights = c(1, 1.4))

ggsave(
  "PBMC_figures/PanelC_D_Bcell_frequency_heatmap.png",
  panel_CD,
  width = 11,
  height = 11,
  dpi = 600
)

ggsave(
  "PBMC_figures/PanelC_D_Bcell_frequency_heatmap.pdf",
  panel_CD,
  width = 11,
  height = 11
)




# -------------------------
# Save object
# -------------------------
saveRDS(pbmc_merged, "PBMC_5lanes_merged_with_demux.rds")
saveRDS(pbmc_filt, "PBMC_5lanes_singlets_filtered_Monaco_annotated.rds")