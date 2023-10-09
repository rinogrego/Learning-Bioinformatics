#install.packages("tidyverse")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("SummarizedExperiment")

library("GEOquery")
library("SummarizedExperiment")

if (!dir.exists("GEO")) {
  dir.create("GEO")
}

# download GEO records
#geo <- getGEO(GEO = "GSE68310", destdir = "GEO", )
#geo <- getGEO(GEO = "GSE103512", destdir = "GEO")

# download manual
# ref: https://www.biostars.org/p/390035/
gse <- getGEO(filename="GEO/GSE68310_series_matrix.txt.gz", AnnotGPL = TRUE)


# checking class
class(gse)
typeof(gse)
isS4(gse)
is(gse, "ExpressionSet")
# inspect gse
gse
# checking attribute
attributes(gse)
# check slot names
slotNames(gse)

# getting SummarizedExperiment data
## 1 sample
#gse68310 <- as(gse[[1]], "SummarizedExperiment")
## entire records (depends on the data shape though)
gse68310 <- as(gse, "SummarizedExperiment")
gse68310
typeof(gse68310)
is(gse68310, "SummarizedExperiment")

# getting MIAME data
miame <- metadata(gse68310)[["experimentData"]]
miame

# getting assay data
assays(gse68310)
class(assay(gse68310, "exprs"))
assay_gse68310 <- assay(gse68310, "exprs")
dim(assay_gse68310)

# checking column
class(rowData(gse68310))
class(colData(gse68310))

# checking the data
metadata(assay(gse68310))
tibble::as_tibble(assay(gse68310, "exprs"))
tibble::as_tibble(rowData(gse68310))
tibble::as_tibble(colData(gse68310))

# save to csv or watever
write.csv2(assay_gse68310, "assay_gse68310.csv")
write.csv2(colData(gse68310), "colData_gse68310.csv")
write.csv2(rowData(gse68310), "rowData_gse68310.csv")

