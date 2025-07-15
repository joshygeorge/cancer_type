library(Biobase)

tismo.dat <- readRDS("./data/tismo/TISMO_expressionvivo_profiles.RDS")
s.anns <- read.csv("./data/tismo/TISMO_vivosample_annotations.csv")
hom <- read.table("./data/homology/HOM_ProteinCoding.rpt",sep="\t")
rownames(s.anns) <- s.anns$SampleName
pData <- new("AnnotatedDataFrame", data= s.anns)
tismo.eset <- ExpressionSet(assayData = as.matrix(tismo.dat), phenoData = pData)

tcga.eset <- readRDS("./data/tcga/tcga_log2_transformed.Rds")
hom <- read.table("./data/homology/HOM_ProteinCoding.rpt",sep="\t")
rownames(hom) <- hom$V2

tismo.sanns <- pData(tismo.eset)
tismo.types <- unique(tismo.sanns$Cancer_type)
tismo.types.mod <- c("LUAD","SKCM","PAAD","GBM","BRCA","BLCA","BRCA","SARC","OV","COAD","MASTO","BRCA","LIHC","PRAD","MESO","KIRC","Follicular","HNSC","STAD")

tismo.tcga <- tismo.sanns$Cancer_type
for(i in 1: length(tismo.types))
{
  tismo.tcga[which(tismo.sanns$Cancer_type == tismo.types[i])] <- tismo.types.mod[i]
}


tismo.sanns$Cancer_type_TCGA <- tismo.tcga
tismo.exp <- exprs(tismo.eset)

genes.sel <- intersect(rownames(tismo.exp),rownames(hom))
hom.sel <- hom[genes.sel,]
tismo.exp.sel <- tismo.exp[genes.sel,]
rownames(tismo.exp.sel) <- hom.sel$V5

pData <- new("AnnotatedDataFrame", data= tismo.sanns)
tismo.tcga.eset <- ExpressionSet(assayData = as.matrix(tismo.exp.sel), phenoData = pData)
saveRDS(tismo.tcga.eset, file = "./data/tismo/tismo_data_tcga_compatible.Rds")


