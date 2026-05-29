# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("dittoSeq")
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)
library(SingleCellExperiment)
library(SingleR)
library(celldex)
library(dittoSeq)
library(dplyr)
library(scales)
library(readr)
library(stringr)
library(purrr)
library(tidyr)

options(Seurat.object.assay.version = "v3")

dir.create("AG_figures", showWarnings = FALSE)

# -------------------------
# Input files
# config 6 = AG negative
# config 7 = AG positive
# -------------------------
AG_files <- c(
  AG_neg = "./DATA/multi_config_6_sample_filtered_feature_bc_matrix.h5",
  AG_pos = "./DATA/multi_config_7_sample_filtered_feature_bc_matrix.h5"
)

demux_files <- c(
  AG_neg = "./DEMUX/multi_config_6_demux/assignments_refined.tsv.gz",
  AG_pos = "./DEMUX/multi_config_7_demux/assignments_refined.tsv.gz"
)

sample_to_config <- c(
  AG_neg = "multi_config_6",
  AG_pos = "multi_config_7"
)

# -------------------------
# Read one AG H5
# -------------------------
read_ag_h5 <- function(h5_file, sample_name) {
  
  x <- Read10X_h5(h5_file)
  print(sample_name)
  print(names(x))
  
  obj <- CreateSeuratObject(
    counts = x[["Gene Expression"]],
    project = sample_name
  )
  
  obj$orig.ident <- sample_name
  obj$condition <- ifelse(sample_name == "AG_pos", "Antigen positive", "Antigen negative")
  obj$config <- sample_to_config[[sample_name]]
  obj$Lane <- ifelse(sample_name == "AG_pos", "7", "6")
  
  if ("Antibody Capture" %in% names(x)) {
    obj[["ADT"]] <- CreateAssayObject(counts = x[["Antibody Capture"]])
  }
  
  return(obj)
}

AG_list <- mapply(
  FUN = read_ag_h5,
  h5_file = AG_files,
  sample_name = names(AG_files),
  SIMPLIFY = FALSE
)

# -------------------------
# Merge AG+ and AG-
# -------------------------
alldata <- merge(
  AG_list[[1]],
  y = AG_list[[2]],
  add.cell.ids = names(AG_list)
)

print(table(alldata$orig.ident))
print(table(alldata$condition))

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
      demux_sample = sample_name
    ) %>%
    distinct(cell, .keep_all = TRUE)
}

demux_tbl <- purrr::map2_dfr(
  demux_files,
  names(demux_files),
  read_demux_one
)

print(head(demux_tbl))
print(table(demux_tbl$demux_sample))
print(table(demux_tbl$subject, useNA = "ifany"))

# -------------------------
# Add demux metadata
# -------------------------
meta <- alldata@meta.data %>%
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

meta2 <- meta2[match(colnames(alldata), meta2$cell), , drop = FALSE]

alldata$subject <- meta2$subject
alldata$demux_sample <- meta2$demux_sample
alldata$SNP_call <- meta2$SNP_call

print(table(alldata$condition, alldata$SNP_call, useNA = "ifany"))
print(table(alldata$subject, useNA = "ifany"))

# -------------------------
# QC metrics
# -------------------------
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")

hb_genes <- grep("^HB", rownames(alldata), value = TRUE)
hb_genes <- hb_genes[!grepl("^HBP", hb_genes)]

alldata <- PercentageFeatureSet(
  alldata,
  features = hb_genes,
  col.name = "percent_hb"
)

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

p_qc_before <- VlnPlot(
  alldata,
  features = feats,
  group.by = "condition",
  pt.size = 0,
  ncol = 3
) + NoLegend() + ggtitle("AG RNA QC Before Filtering")

# -------------------------
# ADT QC
# -------------------------
if ("ADT" %in% names(alldata@assays)) {
  
  adt_counts <- GetAssayData(alldata, assay = "ADT", layer = "counts")
  
  alldata$nCount_ADT <- Matrix::colSums(adt_counts)
  alldata$nFeature_ADT <- Matrix::colSums(adt_counts > 0)
  
  print(table(alldata$condition, alldata$nCount_ADT > 0))
  print(table(alldata$condition, alldata$nCount_ADT > 10))
}

# -------------------------
# Filter cells
# -------------------------
data.filt <- subset(
  alldata,
  subset = SNP_call == "Singleton" &
    nFeature_RNA > 200 &
    percent_mito < 10
)

print(dim(data.filt))
print(table(data.filt$condition))
print(table(data.filt$subject, useNA = "ifany"))

p_qc_after <- VlnPlot(
  data.filt,
  features = feats,
  group.by = "condition",
  pt.size = 0,
  ncol = 3
) + NoLegend() + ggtitle("AG RNA QC After Filtering")

# -------------------------
# Normalize ADT
# -------------------------
if ("ADT" %in% names(data.filt@assays)) {
  DefaultAssay(data.filt) <- "ADT"
  data.filt <- NormalizeData(data.filt, normalization.method = "CLR", margin = 2)
  DefaultAssay(data.filt) <- "RNA"
}

# -------------------------
# RNA processing
# -------------------------
DefaultAssay(data.filt) <- "RNA"

data.filt <- NormalizeData(data.filt)
data.filt <- FindVariableFeatures(data.filt)
data.filt <- ScaleData(data.filt)
data.filt <- RunPCA(data.filt, verbose = FALSE)
data.filt <- FindNeighbors(data.filt, dims = 1:30)
data.filt <- FindClusters(data.filt, resolution = 0.8)
data.filt <- RunUMAP(data.filt, dims = 1:30)

p_cluster <- DimPlot(
  data.filt,
  label = TRUE
) + NoLegend() + ggtitle("AG RNA UMAP Clusters")

p_condition <- DimPlot(
  data.filt,
  group.by = "condition"
) + ggtitle("AG+ vs AG-")

p_subject <- DimPlot(
  data.filt,
  group.by = "subject"
) + ggtitle("AG demux subjects")

# -------------------------
# Monaco annotation
# -------------------------
DefaultAssay(data.filt) <- "RNA"

monaco.ref <- celldex::MonacoImmuneData()
sce <- as.SingleCellExperiment(data.filt)

monaco.main <- SingleR(
  test = sce,
  assay.type.test = 1,
  ref = monaco.ref,
  labels = monaco.ref$label.main
)

monaco.fine <- SingleR(
  test = sce,
  assay.type.test = 1,
  ref = monaco.ref,
  labels = monaco.ref$label.fine
)

data.filt$monaco_main <- monaco.main$pruned.labels
data.filt$monaco_fine <- monaco.fine$pruned.labels

print(table(data.filt$monaco_main, useNA = "ifany"))
print(table(data.filt$monaco_fine, useNA = "ifany"))

p_monaco_main <- DimPlot(
  data.filt,
  group.by = "monaco_main",
  label = TRUE,
  repel = TRUE
) + ggtitle("AG Monaco main annotation")

p_monaco_fine <- DimPlot(
  data.filt,
  group.by = "monaco_fine",
  label = TRUE,
  repel = TRUE
) + ggtitle("AG Monaco fine annotation")

# -------------------------
# Panel C-style B cell frequency plot
# -------------------------
bcell_types <- c(
  "Exhausted B cells",
  "Naive B cells",
  "Non-switched memory B cells",
  "Switched memory B cells"
)

bcell_obj <- subset(
  data.filt,
  subset = monaco_fine %in% bcell_types
)

bcell_obj$bcell_type <- bcell_obj$monaco_fine
bcell_obj$bcell_type <- gsub("Exhausted B cells", "Exhausted B cell", bcell_obj$bcell_type)
bcell_obj$bcell_type <- gsub("Naive B cells", "Naive", bcell_obj$bcell_type)
bcell_obj$bcell_type <- gsub("Non-switched memory B cells", "Non-switched Mem", bcell_obj$bcell_type)
bcell_obj$bcell_type <- gsub("Switched memory B cells", "Switched Mem", bcell_obj$bcell_type)

freq_df <- bcell_obj@meta.data %>%
  dplyr::count(subject, condition, bcell_type) %>%
  dplyr::group_by(subject, condition) %>%
  dplyr::mutate(percent = n / sum(n) * 100) %>%
  ungroup()

freq_df$sample_label <- paste(freq_df$subject, freq_df$condition, sep = "_")

pC <- ggplot(freq_df, aes(x = sample_label, y = percent, fill = bcell_type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  ylab("Percentage of cells") +
  xlab("") +
  ggtitle("C. B cell population frequencies by demux subject") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.title = element_blank()
  )

# -------------------------
# Panel D-style marker heatmap
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

bcell_obj$bcell_type <- factor(
  bcell_obj$bcell_type,
  levels = c("Naive", "Non-switched Mem", "Switched Mem", "Exhausted B cell")
)

Idents(bcell_obj) <- "bcell_type"

pD <- DoHeatmap(
  bcell_obj,
  features = markers,
  group.by = "bcell_type",
  size = 3
) + ggtitle("D. B cell marker expression")

# -------------------------
# Save figures
# -------------------------
ggsave(
  "AG_figures/AG_QC_before_filtering.png",
  p_qc_before,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  "AG_figures/AG_QC_after_filtering.png",
  p_qc_after,
  width = 12,
  height = 7,
  dpi = 300
)

ggsave(
  "AG_figures/AG_UMAP_condition_subject.png",
  (p_cluster | p_condition) / p_subject,
  width = 14,
  height = 10,
  dpi = 600
)

ggsave(
  "AG_figures/AG_Monaco_annotation.png",
  p_monaco_main + p_monaco_fine,
  width = 14,
  height = 6,
  dpi = 600
)

ggsave(
  "AG_figures/PanelC_AG_Bcell_population_frequencies.png",
  pC,
  width = 9,
  height = 5,
  dpi = 600
)

ggsave(
  "AG_figures/PanelD_AG_Bcell_marker_heatmap.png",
  pD,
  width = 10,
  height = 6,
  dpi = 600
)

ggsave(
  "AG_figures/PanelC_D_AG_Bcell_frequency_heatmap.png",
  pC / pD,
  width = 11,
  height = 11,
  dpi = 600
)
###################



dir.create("AG_figures", showWarnings = FALSE)

DefaultAssay(data.filt) <- "ADT"

# -------------------------
# AG antigen FeaturePlots
# -------------------------
p_ag_feature <- FeaturePlot(
  data.filt,
  features = c("DV1", "DV2", "DV3", "DV4", "HSA"),
  reduction = "umap",
  cols = c("lightgrey", "red"),
  ncol = 3
)

ggsave(
  "AG_figures/AG_DV_HSA_FeaturePlots.png",
  p_ag_feature,
  width = 12,
  height = 8,
  dpi = 600
)

# -------------------------
# AG antigen signal by condition
# -------------------------
p_ag_vln_condition <- VlnPlot(
  data.filt,
  features = c("DV1", "DV2", "DV3", "DV4", "HSA"),
  group.by = "condition",
  pt.size = 0,
  ncol = 3
)

ggsave(
  "AG_figures/AG_DV_HSA_VlnPlot_by_condition.png",
  p_ag_vln_condition,
  width = 12,
  height = 8,
  dpi = 600
)

# -------------------------
# AG antigen signal by subject
# -------------------------
p_ag_vln_subject <- VlnPlot(
  data.filt,
  features = c("DV1", "DV2", "DV3", "DV4", "HSA"),
  group.by = "subject",
  pt.size = 0,
  ncol = 3
)

ggsave(
  "AG_figures/AG_DV_HSA_VlnPlot_by_subject.png",
  p_ag_vln_subject,
  width = 14,
  height = 8,
  dpi = 600
)

# -------------------------
# DV-positive UMAP
# -------------------------
p_dv_umap <- DimPlot(
  data.filt,
  group.by = "DV_positive",
  cols = c("DV_neg" = "grey80", "DV_pos" = "red")
)

ggsave(
  "AG_figures/AG_DV_positive_UMAP.png",
  p_dv_umap,
  width = 6,
  height = 5,
  dpi = 600
)

# -------------------------
# DV-positive split by subject
# -------------------------
p_dv_subject <- DimPlot(
  data.filt,
  group.by = "DV_positive",
  split.by = "subject",
  cols = c("DV_neg" = "grey80", "DV_pos" = "red"),
  ncol = 4
)

ggsave(
  "AG_figures/AG_DV_positive_by_subject.png",
  p_dv_subject,
  width = 14,
  height = 10,
  dpi = 600
)

# -------------------------
# Monaco annotation split by DV status
# -------------------------
DefaultAssay(data.filt) <- "RNA"

p_monaco_dv <- DimPlot(
  data.filt,
  group.by = "monaco_fine",
  split.by = "DV_positive",
  label = TRUE,
  repel = TRUE
)

ggsave(
  "AG_figures/AG_Monaco_split_by_DV_status.png",
  p_monaco_dv,
  width = 14,
  height = 6,
  dpi = 600
)

# -------------------------
# DV-positive frequency per subject
# -------------------------
dv_subject_freq <- data.filt@meta.data |>
  dplyr::count(subject, condition, DV_positive) |>
  dplyr::group_by(subject, condition) |>
  dplyr::mutate(freq = n / sum(n)) |>
  dplyr::ungroup()

p_dv_subject_freq <- ggplot(
  dv_subject_freq,
  aes(x = subject, y = freq, fill = DV_positive)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  facet_wrap(~ condition) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percent of cells") +
  xlab("Subject") +
  ggtitle("DV-positive frequency per subject") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(
  "AG_figures/AG_DV_positive_frequency_by_subject.png",
  p_dv_subject_freq,
  width = 10,
  height = 5,
  dpi = 600
)



table(data.filt$subject, data.filt$DV_positive)

prop.table(
  table(data.filt$subject, data.filt$DV_positive),
  margin = 1
)




ggsave(
  "AG_figures/AG_DV_HSA_VlnPlot_by_subject.png",
  p_ag_vln_subject,
  width = 14,
  height = 8,
  dpi = 600
)

ggsave(
  "AG_figures/AG_DV_HSA_VlnPlot_by_subject.pdf",
  p_ag_vln_subject,
  width = 14,
  height = 8
)

ggsave(
  "AG_figures/current_plot.png",
  width = 14,
  height = 8,
  dpi = 600
)


f<-FeaturePlot(
  data.filt,
  features = c("DV1"),
  split.by = "subject",
  cols = c("lightgrey", "red"),
  ncol = 4
)

ggsave(
  "AG_figures/DV1_feature_plot.png",
  f,
  width = 12,
  height = 3,
  dpi = 600
)


DimPlot(
  data.filt,
  group.by = "DV_positive",
  split.by = "condition",
  cols = c("grey80", "red")
)


###################


# -------------------------
# DV-positive split by condition
# -------------------------
p_dv_condition <- DimPlot(
  data.filt,
  group.by = "DV_positive",
  split.by = "condition",
  cols = c("DV_neg" = "grey80", "DV_pos" = "red"),
  ncol = 2
) + ggtitle("DV-positive cells by condition")

ggsave(
  "AG_figures/AG_DV_positive_by_condition.png",
  p_dv_condition,
  width = 10,
  height = 5,
  dpi = 600
)

ggsave(
  "AG_figures/AG_DV_positive_by_condition.pdf",
  p_dv_condition,
  width = 10,
  height = 5
)

# -------------------------
# DV1 FeaturePlot split by subject
# -------------------------
p_dv1_subject <- FeaturePlot(
  data.filt,
  features = "DV1",
  split.by = "subject",
  cols = c("lightgrey", "red"),
  ncol = 4
)

ggsave(
  "AG_figures/AG_DV1_split_by_subject.png",
  p_dv1_subject,
  width = 14,
  height = 4,
  dpi = 600
)

ggsave(
  "AG_figures/AG_DV1_split_by_subject.pdf",
  p_dv1_subject,
  width = 14,
  height = 4
)

# -------------------------
# DV2 FeaturePlot split by subject
# -------------------------
p_dv2_subject <- FeaturePlot(
  data.filt,
  features = "DV2",
  split.by = "subject",
  cols = c("lightgrey", "red"),
  ncol = 4
)

ggsave(
  "AG_figures/AG_DV2_split_by_subject.png",
  p_dv2_subject,
  width = 14,
  height = 4,
  dpi = 600
)

# -------------------------
# DV3 FeaturePlot split by subject
# -------------------------
p_dv3_subject <- FeaturePlot(
  data.filt,
  features = "DV3",
  split.by = "subject",
  cols = c("lightgrey", "red"),
  ncol = 4
)

ggsave(
  "AG_figures/AG_DV3_split_by_subject.png",
  p_dv3_subject,
  width = 14,
  height = 4,
  dpi = 600
)

# -------------------------
# DV4 FeaturePlot split by subject
# -------------------------
p_dv4_subject <- FeaturePlot(
  data.filt,
  features = "DV4",
  split.by = "subject",
  cols = c("lightgrey", "red"),
  ncol = 4
)

ggsave(
  "AG_figures/AG_DV4_split_by_subject.png",
  p_dv4_subject,
  width = 14,
  height = 4,
  dpi = 600
)

# -------------------------
# HSA control split by subject
# -------------------------
p_hsa_subject <- FeaturePlot(
  data.filt,
  features = "HSA",
  split.by = "subject",
  cols = c("lightgrey", "blue"),
  ncol = 4
)

ggsave(
  "AG_figures/AG_HSA_split_by_subject.png",
  p_hsa_subject,
  width = 14,
  height = 4,
  dpi = 600
)
# -------------------------
# Save updated object with DV_positive metadata
# -------------------------



# -------------------------
# Save objects
# -------------------------
saveRDS(alldata, "AG_merged_with_demux.rds")
saveRDS(data.filt, "AG_singlets_filtered_Monaco_annotated.rds")
saveRDS(bcell_obj, "AG_Bcell_only_Monaco_demux.rds")