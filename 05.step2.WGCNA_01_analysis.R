#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(WGCNA)
load(file='05.step1.WGCNA.RData')
# args[1] 为项目名称 如LUSC,BRCA
# args[2] 为分析类型 如T,N ??? Seems to be useless but only generate a file for output. 
# args[3] 为样本树的分割阈值

# args[4] 是为样本树的分割阈值minSize, default=10
if (length(args) > 3){
	minSize = as.numeric(args[4])
}else{
	minSize = 10 
}




### The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


print(dim(datExpr))
pdf(paste(args[1],args[2],'WGCNA_SampleTree2.pdf', sep = "_", collapse = NULL),height=8,width=10)
par(cex = 0.6)
par(mar = c(0,8,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = as.numeric(args[3]), col = "red")
dev.off()

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = as.numeric(args[3]), minSize = minSize)
print(table(clust))
print(clust)
# clust 1 contains the samples we want to keep.
keepSamples <- (clust==1)
datExpr <- datExpr[keepSamples, ]
print('keeping samples: ')
print(length(keepSamples))
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

save(sampleTree,datExpr, file="05.step2.WGCNA.RData")
