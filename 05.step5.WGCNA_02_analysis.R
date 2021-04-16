#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(WGCNA)
# args[1] 为项目名称 如LUSC,BRCA
# args[2] 为条件类型 如T,N
# args[3] 为power设置的阈值
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

prefix = args[1]
condition = args[2]
power = args[3]

load(file = '05.step2.WGCNA.RData')

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


net = blockwiseModules(datExpr, power = as.numeric(power),
TOMType = "unsigned", minModuleSize = 30,
reassignThreshold = 0, mergeCutHeight = 0.25,
numericLabels = TRUE, pamRespectsDendro = FALSE,
saveTOMs = TRUE,
saveTOMFileBase = paste(prefix, condition, 'TOM', sep = "_", collapse = NULL),
verbose = 3)


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# open a graphics window
pdf(paste(prefix, condition, power, 'WGCNA_ModuleColor.pdf', sep = "_", collapse = NULL),height=8,width=10)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

x <- data.frame(name = names(moduleLabels),label = moduleLabels, color = moduleColors)
write.table(x,file=paste0(condition, power, ".module_color_gene.list"),sep="\t",quote= FALSE,row.names = FALSE)

save(MEs, net, moduleLabels, moduleColors, geneTree, file = paste(prefix, condition,'WGCNA_networkConstruction_auto_02.RData', sep = "_", collapse = NULL))
