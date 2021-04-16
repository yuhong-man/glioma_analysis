#!/data2/wangb/anaconda2/bin/Rscript
library(DESeq2)
library(ggplot2)
library(edgeR)
args=commandArgs(T)
# args[1] 为项目名称 如LUSC,BRCA
print('args')
print(args)

#############################################################################################
### 函数定义
#############################################################################################
source('/data2/wangb/pipelines/deseq2/code/Utils.R')


#############################################################################################
### 数据预处理
#############################################################################################
print('read.table')
inputData <- as.matrix(read.table("01.inputdata.txt",sep="\t",header = T,row.names=1, check.names=F, stringsAsFactors =F))
print('read inputData done')
# print(inputData)
colData <- read.table("01.coldata.txt", sep='\t', header = T, stringsAsFactors =F)
print('read colData done')
# print(colData)

# print('rownames(colData)')
# print(rownames(colData))
# print('colData$sample_id')
# print(colData$sample_id)
rownames(colData) <- colData$sample_id
### 记得删除01.coldata.txt里面的最后一行

# print(rownames(colData))
# print(colnames(inputData))

print('all(rownames(colData) == colnames(inputData)) ')
all(rownames(colData) == colnames(inputData)) 
print('condition')
print(colData$condition)
colData$condition = factor(colData$condition,c("N","T"))
print('DESeqDataSetFromMatrix')
print(inputData)
print(colData)

if (nrow(inputData) >= 1000){
	nsub = 1000 
}else{
	nsub = as.integer(nrow(inputData) / 2)
}
dds <- DESeqDataSetFromMatrix(countData = inputData, colData, design = ~ condition)
### Plot PCA

pcaData <- my_plotPCA(dds, paste(args[1],'DESeq2_PCA.pdf', sep = "_", collapse = NULL), nsub=nsub)

save(pcaData,inputData,colData, file="04.step1.DESeq2.RData")


# edgeR
inputData <- inputData[rowSums(cpm(inputData) > 1) >= 2,]
inputData <- DGEList(counts = inputData, group = colData$condition)
inputData <- calcNormFactors(inputData)

inputData <- estimateCommonDisp(inputData)
inputData <- estimateTagwiseDisp(inputData)
et <- exactTest(inputData)
tTag <- topTags(et, n=nrow(inputData))
tTag <- as.data.frame(tTag)
write.table(tTag, file = "edgeR.TvsN.list", sep='\t', quote = FALSE)
