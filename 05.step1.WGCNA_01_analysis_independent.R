#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(WGCNA)
# project_name 为项目名称 如LUSC,BRCA
# selected_condition 为分析类型 如T,N
path_inputData = args[1]
path_colData = args[2]
project_name <- args[3]
selected_condition <- args[4]



### The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# load_file <- paste(project_name,'DESeq2_inputData_colData.RData', sep = "_", collapse = NULL)
## 导入变量 inputData colData
# load(file=load_file) 
inputData <- as.matrix(read.table(path_inputData,sep="\t",header = T,row.names=1, check.names=F, stringsAsFactors =F))


colData <- read.table(path_colData, sep='\t', header = T, stringsAsFactors =F)
rownames(colData) <- colData$sample_id
# colData$condition = factor(colData$condition,c("N","T"))
### 分类不同样本条件

# If it's all, ignore this. 
print('selected_condition')
print(selected_condition)
if(selected_condition != 'all'){ 
	class <- which(colData$condition %in% c(selected_condition))
	inputData <- inputData[,class]
}
print(head(inputData))
#############################################################################################
### 数据预处理
#############################################################################################
dataPrepare <- function(inputData){
	datExpr0 = as.data.frame(t(inputData))
	names(datExpr0) = rownames(inputData)
	rownames(datExpr0) = colnames(inputData)
	gsg = goodSamplesGenes(datExpr0, verbose = 3);
	if (!gsg$allOK)
	{
	  # Optionally, print the gene and sample names that were removed:
	  if (sum(!gsg$goodGenes)>0) 
	     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
	  if (sum(!gsg$goodSamples)>0) 
	     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
	  # Remove the offending genes and samples from the data:
	  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
	}
  	return(datExpr0)
}
datExpr <- dataPrepare(inputData)

#############################################################################################
### 画树观察离群点
#############################################################################################
sampleTree <- hclust(dist(datExpr), method = "average")
# pdf("WGCNA_SampleTree.pdf
pdf(paste(project_name,selected_condition,'WGCNA','SampleTree.pdf', sep = "_", collapse = NULL),height=8,width=10)
par(cex = 0.6)
par(mar = c(0,8,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()

datExpr[] <- lapply(datExpr, as.numeric)
save(sampleTree,datExpr, file="05.step1.WGCNA.RData")

