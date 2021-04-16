#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(DESeq2)

# args[1] 为项目名称 如LUSC,BRCA
# args[2] 是否是对已处理的数据进行处理 T为是 F为否
# args[3] group类型
# args[4] PC1 大于'g' 还是小于'l'
# args[5] PC1 值筛选
# args[6] PC2 大于'g' 还是小于'l'
# args[7] PC2 值筛选

if(args[2] == 'T'){
	load(file = "04.step2.DESeq2.RData");
} else if(args[2] == 'F'){
	load(file = "04.step1.DESeq2.RData");	
} else{
	stop("miss the argument 'F'/'T'")
}

source('/data2/wangb/pipelines/deseq2/code/Utils.R')
filter_data <- function(pcaData,group,PC1_d,PC1_v,PC2_d,PC2_v){
	if(PC1_d == 'g'&PC2_d =='l'){
		outlier_name <- pcaData[pcaData$group == group&pcaData$PC1 > PC1_v&pcaData$PC2 < PC2_v,]
	}
	else if(PC1_d == 'g'&PC2_d =='g'){
		outlier_name <- pcaData[pcaData$group == group&pcaData$PC1 > PC1_v&pcaData$PC2 > PC2_v,]
	}
	else if(PC1_d == 'l'&PC2_d =='g'){
		outlier_name <- pcaData[pcaData$group == group&pcaData$PC1 < PC1_v&pcaData$PC2 > PC2_v,]
	}
	else if(PC1_d == 'l'&PC2_d =='l'){
		outlier_name <- pcaData[pcaData$group == group&pcaData$PC1 < PC1_v&pcaData$PC2 < PC2_v,]
	}
	else{
		outlier_name <- ''
		stop("group/PC1_d/PC1_v//PC2_d/PC2_v ERROR")
	}
	return(outlier_name)
}
### 手动调整
### 剔除一些离群值 "T" "40" 为自选项

outlier_name <- filter_data(pcaData,args[3],args[4],as.numeric(args[5]),args[6],as.numeric(args[7]))
outlier <- which(colnames(inputData) %in% row.names(outlier_name))
print('length(outlier)')
print(length(outlier))
print('inputData before filtering outlier')
print(dim(inputData))
if(length(outlier) > 0){
	inputData <- inputData[,-outlier]
	colData <- colData[-outlier,]
}

print('dim(inputData)')
print(dim(inputData))
dds <- DESeqDataSetFromMatrix(countData = inputData, colData, design = ~ condition)
library(ggplot2)

if (nrow(inputData) >= 1000){
	nsub = 1000 
}else{
	nsub = as.integer(nrow(inputData) / 2)
}

pcaData <- my_plotPCA(dds, paste(args[1],'DESeq2_PCA_clear.pdf', sep = "_", collapse = NULL), nsub=nsub)

save(pcaData,inputData,colData, file="04.step2.DESeq2.RData")
