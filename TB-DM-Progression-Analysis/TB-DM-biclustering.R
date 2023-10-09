# link dataset: https://ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193978

library(GEOquery)
library(SummarizedExperiment)
library(dplyr)


if (!dir.exists("GEO")) {
  dir.create("GEO")
}

# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE193978", "file=GSE193978_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

#View(tbl)

# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]

# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
dat <- log10(tbl + 1)

#View(dat)

# check platform
#gpl <- getGEO(filename = "./GEO/GPL18573.soft.gz")
#attributes(gpl)
#View(gpl)

# read metadata
gse <- getGEO(GEO = "GSE193978", GSEMatrix = TRUE, destdir = "GEO")

metadata <- pData(phenoData(gse[[1]]))
View(head(metadata))
metadata <- metadata %>% 
  select(2, 40, 41, 42, 43) %>%
  rename(disease_category = "disease_category:ch1") %>%
  rename(field_site = "field_site:ch1") %>%
  rename(patient_id = "patient id:ch1") %>%
  rename(timepoint = "timepoint:ch1")

if (!dir.exists("GSE193978")) {
  dir.create("GSE193978")
}

write.csv2(metadata, "./GSE193978/GSE193978_metadata.csv")
write.csv2(dat, "./GSE193978/GSE193978_dat.csv")


# biclustering
library(pheatmap)

# filter by gene counts
dat_sample <- dat[1:100, ]
pheatmap(
  dat_sample,
  cluster_rows = T,
  cluster_cols = T,
)

# add annotation color bar
annot_cols <- data.frame(
  Country = metadata[colnames(dat), "field_site"],
  Disease_Category = metadata[colnames(dat), "disease_category"],
  row.names = colnames(dat)
)

dat_sample <- dat[1:30, 0:100]
pheatmap(
  dat_sample,
  cluster_rows = T,
  cluster_cols = T,
  annotation_col = annot_cols
)

# customize annotation color
# Ref: https://stackoverflow.com/questions/48678910/change-colors-in-multiple-annotations
annot_colors <- list(
  Disease_Category = c("TB-DM"="red", "TB-IH"="pink", "TB-only"="darkred", "TB-preDM"="orange"),
  Country = c(Indonesia="red", Romania="yellow", South_Africa="green")
)

dat_sample <- dat[1:100, 0:200]
# color palette ref: https://stackoverflow.com/questions/71298942/making-a-continuous-color-chart-for-heatmap-using-pheatmap
pheatmap(
  dat_sample,
  #scale = "row",
  cluster_rows = T,
  cluster_cols = T,
  annotation_col = annot_cols,
  annotation_colors = annot_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("darkgreen", "darkred"))(100)
)



# analyzing data from Indonesia

metadata_indo <- metadata[colnames(dat), ]
metadata_indo <- metadata_indo[metadata_indo$field_site=="Indonesia", ]
dat_sample_indo = dat[1:100, rownames(metadata_indo)]

annot_cols <- data.frame(
  Disease_Category = metadata[colnames(dat_sample_indo), "disease_category"],
  row.names = colnames(dat_sample_indo)
)

annot_colors_indo <- list(
  Disease_Category = c("TB-DM"="red", "TB-IH"="pink", "TB-only"="darkred", "TB-preDM"="orange")
)

pheatmap(
  dat_sample_indo,
  main = "TB-DM Progression Analysis (Indonesia)",
  scale = "row",
  cluster_rows = T,
  cluster_cols = T,
  annotation_col = annot_cols,
  annotation_colors = annot_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("darkred", "darkgreen"))(100),
  cutree_rows = 3,
  cutree_cols = 4,
)

# additional biclustering ref: https://www.youtube.com/watch?v=crkXYfc0tf0
