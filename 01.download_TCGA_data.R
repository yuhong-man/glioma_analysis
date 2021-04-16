#!/data2/wangb/anaconda2/bin/Rscript
# example
# Rscript 01.download_TCGA_data.R BRCA 20160128


args=commandArgs(T)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RTCGAToolbox")

# 下载乳腺癌数据
library(RTCGAToolbox)
# brcaData = getFirehoseData (dataset="BRCA", runDate="20160128",RNASeqGene=TRUE)
brcaData = getFirehoseData (dataset=args[1], runDate=args[2],RNASeqGene=TRUE)
