
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")
BiocManager::install("EDASeq")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(EDASeq)

# load tumor and normal transcription data on glioblastoma from TCGA
query.exp <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor","Solid Tissue Normal")
)

GDCdownload(query=query.exp,files.per.chunk = 100)

gbm.exp <- GDCprepare(
  query = query.exp, 
  save = F
)

# save barcodes of experiments
samplesDown <- getResults(query.exp,cols=c("cases"))

# save barcodes for tumor samples
dataTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                typesample = "TP")

# save barcodes for "healthy" samples
dataNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                typesample = "NT")

# remove duplicated samples from the same patients
i <- which(duplicated(substr(dataTP,1,12)))
dataTP <- dataTP[-i] 

# Array Array Intensity correlation analysis to identify possible outliers among samples
dataPrep <- TCGAanalyze_Preprocessing(gbm.exp)
rownames(dataPrep) <- gsub("\\..*", "", rownames(dataPrep))
dataPrep <- dataPrep[ , colnames(dataPrep) %in% c(dataTP, dataNT)]

# load data on gene identifiers and types
gene.info <- rowRanges(gbm.exp)
gene.info <- gene.info@elementMetadata@listData
gene.info <- data.frame(id = gene.info$gene_id,
                      gene=gene.info$gene_name,
                      type=gene.info$gene_type
                      )
gene.info$id<-gsub("\\..*", "", gene.info$id)


# selection of the protein-coding genes for further analysis
gene.info<-unique(gene.info[gene.info$type=="protein_coding",-3])
gene.ids<-unique(gene.info$id)
dataPrep <- dataPrep[rownames(dataPrep) %in% gene.ids,]


# normalization of transcription profiles
dataNorm <- TCGAanalyze_Normalization(
  tabDF = dataPrep,
  geneInfo = geneInfoHT,
  method = "gcContent"
)

# filter genes with low expression
dataFilt <- TCGAanalyze_Filtering(
  tabDF = dataNorm,
  method = "quantile", 
  qnt.cut =  0.25
)

# create table with additional data on samples
sample_id <- substr(colnames(dataFilt),1,12)
col_data <- data.frame(colData(gbm.exp),stringsAsFactors = F)
col_data <- col_data[col_data$barcode %in% colnames(dataFilt),]
i <- match(sample_id,col_data$patient)
col_data <- col_data[i,]

# create table with additional data on genes
gene.info <- rowRanges(gbm.exp)
gene.info <- gene.info@elementMetadata@listData
gene.info <- as.data.frame(gene.info)
gene.info <- gene.info[,5:7]
gene.info$gene_id <- gsub("\\..*", "", gene.info$gene_id)
gene.info <- unique(gene.info)
gene.info <- gene.info[gene.info$gene_id %in% rownames(dataFilt),]
i <- match(rownames(dataFilt),gene.info$gene_id)
gene.info <- gene.info[i,]

# save all data
save(dataFilt,col_data,gene.info,file="GBM_transcription_data.RData")
