#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# args[1] 为项目名称 如LUSC,BRCA
# args[2] 为条件类型 如T,N
# args[3] 临床数据
clinical_trait_of_interest = args[3]
enableWGCNAThreads()
# Load the expression and trait data saved in the first part
load(file = paste(args[1],args[2],'step3_WGCNA_01.RData', sep = "_", collapse = NULL))
load(file = paste(args[1],args[2],'WGCNA_networkConstruction_auto_02.RData', sep = "_", collapse = NULL))

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
print('datTraits')
print(datTraits)

# Isolate weight from the clinical traits
traitArgs = as.data.frame(datTraits[clinical_trait_of_interest]);
print('traitArgs') 
print(traitArgs)
names(traitArgs) = clinical_trait_of_interest
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, traitArgs))
# Plot the relationships among the eigengenes and the trait
pdf(paste(args[1],args[2],'WGCNA_EigengeneNetworks.pdf', sep = "_", collapse = NULL),height=8,width=10)
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Plot the dendrogram
pdf(paste(args[1],args[2],'Eigengene_dendrogram.pdf', sep = "_", collapse = NULL),height=8,width=10)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
dev.off()


pdf(paste(args[1],args[2],'Eigengene_adjacency_heatmap.pdf', sep = "_", collapse = NULL),height=8,width=10)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

