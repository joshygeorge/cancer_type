library(Biobase)
library(ggplot2)
library("factoextra")
library(scatterplot3d)
library(plotly)



tcga.eset <- readRDS("./data/TCGA/tcga_log2_transformed.Rds")
tismo.eset <- readRDS("./data/TISMO/tismo_data_tcga_compatible.Rds")

common.cancer.types <- intersect(tcga.eset$type, tismo.eset$Cancer_type_TCGA)
common.cancer.types <- setdiff(common.cancer.types,"MESO")
#common.cancer.types  <- c("BRCA","SKCM")
tismo.eset.sel <- tismo.eset[,tismo.eset$Cancer_type_TCGA %in% common.cancer.types]
tcga.eset.sel <- tcga.eset[,tcga.eset$type %in% common.cancer.types]


tcga.probes <- featureData(tcga.eset)@data
tcga.exp <- exprs(tcga.eset.sel)
tismo.sel.exp <- exprs(tismo.eset.sel)
genes <- intersect(tcga.probes$symbol,rownames(tismo.sel.exp))

tcga.probes.sel <- tcga.probes[tcga.probes$symbol %in% genes, ]

tcga.final <- tcga.exp[rownames(tcga.probes.sel),]
rownames(tcga.final) <- tcga.probes.sel$symbol

tcga.final.eset <- ExpressionSet(assayData = as.matrix(tcga.final), phenoData = new("AnnotatedDataFrame", data= pData(tcga.eset.sel)))
tismo.final.exp <- tismo.sel.exp[tcga.probes.sel$symbol,]
tismo.final.eset <- ExpressionSet(assayData = as.matrix(tismo.final.exp), phenoData = new("AnnotatedDataFrame", data= pData(tismo.eset.sel)))


saveRDS(tismo.final.eset, file = "./data/TISMO/tismo_data_matched_tumor_type_and_genes.Rds")
saveRDS(tcga.final.eset, file = "./data/TCGA/tcga_data_matched_tumor_type_and_genes.Rds")

