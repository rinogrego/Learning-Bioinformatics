# DATASET LINK: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103512

BiocManager::install("GEOquery")
BiocManager::install("SummarizedExperiment")

library(GEOquery)
library(SummarizedExperiment)

if (!dir.exists("GEO")) {
  dir.create("GEO")
}

gse <- getGEO(GEO = "GSE103512", destdir = "GEO")
#gse <- getGEO(filename="GEO/GSE103512_series_matrix.txt.gz", AnnotGPL = TRUE)

class(gse)
View(gse)

gse103512 <- as(gse$GSE103512_series_matrix.txt.gz, "SummarizedExperiment")
class(assay(gse103512, "exprs"))
assay_gse103512 <- assay(gse103512, "exprs")
dim(assay_gse103512)
View(assay_gse103512)
View(rowData(gse103512))
View(colData(gse103512))

if (!dir.exists("GSE103512")) {
  dir.create("GSE103512")
}
# write.cv2: separator is ';' (titik koma)
write.csv2(assay_gse103512, "GSE103511/assay_gse103512.csv")
write.csv2(colData(gse103512), "GSE103511/colData_gse103512.csv")
write.csv2(rowData(gse103512), "GSE103511/rowData_gse103512.csv")
