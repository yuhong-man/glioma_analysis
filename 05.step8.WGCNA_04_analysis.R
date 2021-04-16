#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(WGCNA)
enableWGCNAThreads()
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# args[1] 为项目名称 如LUSC,BRCA
# args[2] 为条件类型 如T,N
# args[3] power阈值


# Load the expression and trait data saved in the first part
load(file = paste(args[1],args[2],'step3_WGCNA_01.RData', sep = "_", collapse = NULL))
load(file = paste(args[1],args[2],'WGCNA_networkConstruction_auto_02.RData', sep = "_", collapse = NULL))

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = as.numeric(args[3]));
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^as.numeric(args[3]);
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

save(plotTOM, file = paste(args[1],args[2],'05.step8.WGCNA.TOM.RData', sep = "_", collapse = NULL))


if(nGenes < 5000){
	# To avoid long time processing. 
	# Call the plot function
	print('running the full TOM')
	pdf(paste(args[1],args[2],'WGCNA_TOM.pdf', sep = "_", collapse = NULL),height=15,width=15)
	# sizeGrWindow(9,9)
	TOMplot(plotTOM, net$dendrograms, moduleColors, main = "Network heatmap plot, all genes")
	dev.off()
}else{
	print('too many genes, not running the full TOM. ')
}
#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^as.numeric(args[3]);
diag(plotDiss) = NA;
save(plotDiss, file = paste(args[1],args[2],'05.step8.WGCNA.Diss.RData', sep = "_", collapse = NULL))


pdf(paste(args[1],args[2],'WGCNA_TOM_random.pdf', sep = "_", collapse = NULL),height=15,width=15)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()
