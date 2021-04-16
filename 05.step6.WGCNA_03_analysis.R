#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# args[1] 为项目名称 如LUSC,BRCA
# args[2] 为条件类型 如T,N


# Load the expression and trait data saved in the first part
load(file = paste(args[1],args[2],'step3_WGCNA_01.RData', sep = "_", collapse = NULL))
load(file = paste(args[1],args[2],'WGCNA_networkConstruction_auto_02.RData', sep = "_", collapse = NULL))



#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                           signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(paste(args[1],args[2],'WGCNA_Module_Trait_Relationships.pdf', sep = "_", collapse = NULL), height=30,width=15)
par(mar = c(15, 12, 3, 3));

print('moduleTraitCor')
print(moduleTraitCor)
print('textMatrix')
print(textMatrix)
print('datTraits')
print(datTraits)
print('MEs')
print(MEs)

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               yColorWidth = 2 * strwidth("M"),
               xLabelsAngle = 15,
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

save(datTraits, MEs, datExpr, nSamples, file = paste(args[1],args[2],'05.step6.WGCNA.RData', sep = "_", collapse = NULL))

