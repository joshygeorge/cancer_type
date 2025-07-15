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


write.csv(tismo.exp.nor,"tismo_data_normalized.csv")
write.csv(tcga.exp.nor,"tcga_data_normalized.csv")


