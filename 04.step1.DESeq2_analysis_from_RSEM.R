#!/data2/wangb/anaconda2/bin/Rscript
library(DESeq2)
library(tximport)
args=commandArgs(T)
# args[1] 为项目名称 如LUSC,BRCA
print('args')
print(args)

#############################################################################################
### 函数定义
#############################################################################################
my_plotPCA <- function(dds,filename){
	vsd <- vst(dds,blind=FALSE)
	pcaData <- plotPCA(vsd,intgroup=c("condition"),returnData=TRUE)
	percentVar <- round(100*attr(pcaData,"percentVar"))
	PcaPlotFig <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
	  geom_point(size=3) +
	  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	  coord_fixed()
	ggsave(PcaPlotFig, filename = filename)
  	return(pcaData)
}



#############################################################################################
### 数据预处理
#############################################################################################
# print('read.table')
# inputData <- as.matrix(read.table("01.inputdata.txt",sep="\t",header = T,row.names=1, check.names=F, stringsAsFactors =F))
# Changed to: 
# https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#rsem
print('read inputData done')
# print(inputData)
colData <- read.table("01.coldata.txt", sep='\t', header = T, stringsAsFactors =F)
#
files <- file.path('results', paste0(colData$sample_id, ".genes.results"))

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE) 
inputData <- txi.rsem$counts


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
dds <- DESeqDataSetFromTximport(txi.rsem, colData, design = ~ condition)
### Plot PCA
library(ggplot2)
pcaData <- my_plotPCA(dds,paste(args[1],'DESeq2_PCA.pdf', sep = "_", collapse = NULL))

save(pcaData,inputData,colData, file="04.step1.DESeq2.RData")