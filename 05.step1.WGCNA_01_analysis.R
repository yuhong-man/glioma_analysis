#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(WGCNA)
# project_name 为项目名称 如LUSC,BRCA
# selected_condition 为分析类型 如T,N
project_name <- args[1]
selected_condition <- args[2]

path_use_gene_list <- args[3]


### The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
load_file <- paste(project_name,'DESeq2_inputData_colData.RData', sep = "_", collapse = NULL)
### 导入变量 inputData colData
load(file=load_file) 
### 分类不同样本条件

# If it's all, ignore this. 
if(selected_condition != 'all'){ 
	class <- which(colData$condition %in% c(selected_condition))
	inputData <- inputData[,class]
}
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


if(file.exists(path_use_gene_list)){
	gene_list = read.table(path_use_gene_list)['V1']
	gene_list = unlist(gene_list)
}else{
	gene_list = rownames(inputData)
}

inputData <- inputData[gene_list,]

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
# To avoid integer exception in blockwiseModules
datExpr[] <- lapply(datExpr, as.numeric)
save(sampleTree,datExpr, file="05.step1.WGCNA.RData")

