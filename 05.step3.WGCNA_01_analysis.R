#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(WGCNA)
load(file='05.step2.WGCNA.RData')
# args[1] 为项目名称 如LUSC,BRCA
# args[2] 为条件类型 如T,N
# args[3] 剔除的临床字段，如1#2#3 字段为clinical表里对应字段行数-1
project_name = args[1]
condition_type = args[2] 
clinical_fields_to_remove = args[3]


#############################################################################################
### 追踪临床数据
#############################################################################################
path_clinical_file = "01.clinical.txt"
if(file.exists(path_clinical_file)){
	traitData <- read.table(path_clinical_file,sep="\t",row.names=1)
	# remove columns that hold information we do not need.
	# allTraits <- traitData[, -c(2,7,8,9,10,11,12,15,17,20,21,22)]; #LUSC
	# allTraits <- traitData[, -c(2,7,8,9,10,11,12,15,16,18,19)]; #BRCA_T 
	
	# There might be spaces in columns or index. 
	colnames(traitData) <- gsub(' ', '_', colnames(traitData))
	rownames(traitData) <- gsub(' ', '_', rownames(traitData))
	print('colnames(traitData) replaced')
	print(colnames(traitData))
	traitData <- as.data.frame(t(traitData))

	dim(traitData)
	names(traitData)

	remove_list <- as.numeric(unlist(strsplit(clinical_fields_to_remove,split='#')))
	allTraits <- traitData[, -remove_list]; #BRCA_N
}else{
	# Not exists, just create a empty dataframe with HYBRIDIZATION_REF field. 
	allTraits <- data.frame(rownames(datExpr))
	colnames(allTraits) <- c('HYBRIDIZATION_REF')
	allTraits[,'sex'] <- seq(nrow(datExpr))
	allTraits[,'attribute1'] <- c('HYBRIDIZATION_REF')
	allTraits[,'attribute2'] <- c('HYBRIDIZATION_REF')
}


dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
print('datExpr')
# print(datExpr)
print('rownames(datExpr)')
print(rownames(datExpr))

# print('allTraits$HYBRIDIZATION_REF')
# print(allTraits['HYBRIDIZATION_REF'])


traitRows <- match(rownames(datExpr),allTraits$HYBRIDIZATION_REF);
print('traitRows')
print(traitRows)
# datTraits <- allTraits[traitRows, -1] # Why -1? 
datTraits <- allTraits[traitRows,]  
datTraits <- as.data.frame(datTraits)
print('datTraits')
print(datTraits)
rownames(datTraits) = allTraits[traitRows, 1];
# datTraits <- as.data.frame(lapply(datTraits,as.numeric))
# ??? failed. 

datTraits <- as.data.frame(lapply(datTraits,as.factor))
datTraits <- as.data.frame(lapply(datTraits,as.numeric))
collectGarbage();

# Re-cluster samples
sampleTree2 <- hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors <- numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.
pdf(paste(args[1],args[2],'WGCNA_heatmap.pdf', sep = "_", collapse = NULL),height=8,width=10)
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap",marAll = c(4, 5, 3, 1))
dev.off()

# save(datExpr, datTraits, file = "FemaleLiver-01-dataInput.RData")
save(datExpr, datTraits,file = paste(args[1],args[2],'step3_WGCNA_01.RData', sep = "_", collapse = NULL))