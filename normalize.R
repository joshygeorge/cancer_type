library(Biobase)
library(ggplot2)
library(Matrix)
library(data.table)
library("factoextra")
library(scatterplot3d)
library(plotly)
library(RNOmni)

tismo.eset <- readRDS("./data/TISMO/tismo_data_matched_tumor_type_and_genes.Rds")
tcga.eset  <- readRDS("./data/TCGA/tcga_data_matched_tumor_type_and_genes.Rds")

tcga.exp <- exprs(tcga.eset)
tismo.exp <- exprs(tismo.eset)
tcga.anns <- pData(tcga.eset)
tismo.anns <- pData(tismo.eset)


tcga.exp.nor <- tcga.exp
for(i in 1: ncol(tcga.exp))
{
  tcga.exp.nor[,i] <- RankNorm(tcga.exp[,i])
}


tismo.exp.nor <- tismo.exp
for(i in 1: ncol(tismo.exp))
{
  tismo.exp.nor[,i] <- RankNorm(tismo.exp[,i])
}

tcga.exp.t <- data.frame(t(tcga.exp.nor),CANCER_TYPE = tcga.anns$type)
tismo.exp.t <- data.frame(t(tismo.exp.nor),CANCER_TYPE = tismo.anns$Cancer_type_TCGA)


write.csv(tismo.exp.t,"tismo_data_normalized.csv")
write.csv(tcga.exp.t,"tcga_data_normalized.csv")


