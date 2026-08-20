#-----
# 1. Install and load packages
#-----

library(DESeq2)
library(dplyr)
library(readr)
library(tibble)
library(Matrix)

#-----
# 2. Load raw counts and metadata
#-----

raw_path <- "GSM4304272_IL13_chronic_raw_matrix.txt.gz"
meta_path <- "GSM4304272_IL13_chronic_metadata.txt.gz"

raw_counts <- read.delim(gzfile(raw_path),row.names = 1,check.names = FALSE)
metadata <- read.delim(gzfile(meta_path),row.names = 1,check.names = FALSE)

dim(raw_counts)
dim(metadata)
head(metadata)
colnames(metadata)

#-----
# 3. Check that cell names match + align metadata
#-----

head(colnames(raw_counts))
head(rownames(metadata))
sum(colnames(raw_counts) %in% rownames(metadata))
ncol(raw_counts)
colnames(metadata)
rownames(raw_counts)

metadata <- metadata[colnames(raw_counts), ]
stopifnot(all(colnames(raw_counts) == rownames(metadata)))

# rename some columns for my own sake
metadata$condition <- metadata$treatment
metadata$sample_id <- paste(metadata$donor, metadata$condition, sep = "_")
colnames(metadata)

#-----
# 4. Map clusters to cell types, based on the original cluster mapping from Jackson et al. (2020)
#-----

metadata$celltype <- "other"
metadata$celltype[metadata$clusters %in% c("c1", "c2", "c3", "c4")] <- "ciliated"
metadata$celltype[metadata$clusters %in% c("c7", "c8")] <- "goblet"
metadata$celltype[metadata$clusters %in% c("c5", "c6")] <- "club"

metadata$sample_id <- paste(metadata$celltype,metadata$donor,metadata$condition,sep = "_")

table(metadata$celltype, metadata$clusters)
table(metadata$celltype, metadata$condition)
table(metadata$celltype, metadata$donor)
table(metadata$sample_id)

#-----
# 5. Pseudobulk analysis
#-----

run_pseudobulk_deseq2_matrix <- function(raw_counts,metadata,target_celltype,output_dir = "deg_deseq2") {
  
  message("Running ", target_celltype)
  
  cells_keep <- rownames(metadata)[metadata$celltype == target_celltype]
  
  counts_sub <- raw_counts[, cells_keep, drop = FALSE]
  meta_sub <- metadata[cells_keep, , drop = FALSE]
  
  sample_ids <- unique(meta_sub$sample_id)
  
  pb_counts <- sapply(sample_ids, function(s) {
    cells <- rownames(meta_sub)[meta_sub$sample_id == s]
    Matrix::rowSums(as.matrix(counts_sub[, cells, drop = FALSE]))
  })
  
  pb_counts <- as.matrix(pb_counts)
  
  coldata <- meta_sub %>%
    distinct(sample_id, condition, donor) %>%
    as.data.frame()
  
  rownames(coldata) <- coldata$sample_id
  
  pb_counts <- pb_counts[, rownames(coldata)]
  
  coldata$condition <- factor(
    coldata$condition,
    levels = c("BSA", "IL13")
  )
  
  coldata$donor <- factor(coldata$donor)
  
  keep <- rowSums(pb_counts) >= 10
  pb_counts <- pb_counts[keep, ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_counts),
    colData = coldata,
    design = ~ donor + condition
  )
  
  dds <- DESeq(dds)
  
  res <- results(
    dds,
    contrast = c("condition", "IL13", "BSA")
  )
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)
  
  res_df$celltype <- target_celltype
  
  out_path <- file.path(
    output_dir,
    paste0("jackson_pseudobulk_deseq2_", target_celltype, ".csv")
  )
  
  write_csv(res_df, out_path)
  
  message("Saved: ", out_path)
  
  return(res_df)
}

#-----
# 6. Run it
#-----

output_dir <- "deg_deseq2"

goblet_res <- run_pseudobulk_deseq2_matrix(
  raw_counts,
  metadata,
  "goblet",
  output_dir
)

club_res <- run_pseudobulk_deseq2_matrix(
  raw_counts,
  metadata,
  "club",
  output_dir
)

ciliated_res <- run_pseudobulk_deseq2_matrix(
  raw_counts,
  metadata,
  "ciliated",
  output_dir
)

all_res <- bind_rows(
  goblet_res,
  club_res,
  ciliated_res
)

write_csv(
  all_res,
  "jackson_pseudobulk_deseq2_all_celltypes.csv"
)

# check
head(all_res)
table(all_res$celltype)
sum(all_res$padj < 0.05, na.rm = TRUE)

