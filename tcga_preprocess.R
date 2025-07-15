
library("readxl")
library(Biobase)

exp.dat <- read.csv("./data/tcga/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv",sep="\t",row.names=1)
clin.dat <- as.data.frame(read_excel("./data/tcga/TCGA-CDR-SupplementalTableS1.xlsx", col_types = "text"))
clin.dat$ID <- gsub("-",".", clin.dat$bcr_patient_barcode, fixed=T)
sample.type <- substr(colnames(exp.dat),14,15)

#tum.exp <- exp.dat[, sample.type %in% c('01','02', '03', '05', '06' )]
tum.exp <- exp.dat[, sample.type %in% c('01')]
tum.id  <- substr(colnames(tum.exp),1,12)
tum.map <- data.frame(ext.name = colnames(tum.exp), ID = tum.id)
tum.map.uniq <- tum.map[!duplicated(tum.map$ID),]
tum.exp.uniq <- tum.exp[, tum.map.uniq$ext.name]

clin.data.final <- merge(tum.map.uniq,clin.dat,all.x =F)
tum.data.final <- tum.exp[,clin.data.final$ext.name]
tum.data.final[is.na(tum.data.final)] <- 0.001
tum.data.final[tum.data.final < 0.001] <- 0.001
tum.data.log2 <- log2(tum.data.final)
rownames(clin.data.final) <- clin.data.final$ext.name

pData <- new("AnnotatedDataFrame", data= clin.data.final)

symbol <- sapply(rownames(tum.data.final), function(dat) { strsplit(as.character(dat),"|",fixed=T)[[1]][1]})
entid <- sapply(rownames(tum.data.final), function(dat) { strsplit(as.character(dat),"|",fixed=T)[[1]][2]})
probe.data <- data.frame(symbol,entid)
rownames(probe.data) <- rownames(tum.data.final)

annData <- new("AnnotatedDataFrame", data= probe.data)
tcga.eset <- ExpressionSet(assayData = as.matrix(tum.data.log2), phenoData = pData , featureData  = annData)

saveRDS(tcga.eset, file = "./data/tcga/tcga_log2_transformed.Rds")




