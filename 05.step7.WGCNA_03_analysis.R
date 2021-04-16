#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# project_name 为项目名称 如LUSC,BRCA
# condition 为条件类型 如T,N
# trait_column_of_interest 选择感兴趣的表现数据
# args[4] module???    

project_name <- args[1]
condition <- args[2]
trait_column_of_interest <- args[3] 
module <- args[4] 


# Adjust the loading order of the two load().
load(file = paste(project_name, condition, 'WGCNA_networkConstruction_auto_02.RData', sep = "_", collapse = NULL))

load(file = paste(project_name, condition, '05.step6.WGCNA.RData', sep = "_", collapse = NULL))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits[trait_column_of_interest]);
# print('weight')
# print(weight)
names(weight) = trait_column_of_interest
# names (colors) of the modules


# Added on 20200629. 
# MEs <- MEs[rownames(datExpr),]
print('dim(MEs)')
print(dim(MEs))

print('modNames before substring')
print(names(MEs))

modNames = substring(names(MEs), 3)
print('modNames after substring')
print(modNames)

print('corPvalueStudent')
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");

print('geneModuleMembership')
# print(geneModuleMembership)

names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


print('module') 
print(module)
print('modNames')
print(modNames)

column = match(module, modNames);
moduleGenes = moduleColors==module;
print('moduleGenes')
print(sum(moduleGenes))
print('column')
print(column)

pdf(paste(project_name,condition,'Module_membership_vs_gene_significance.pdf', sep = "_", collapse = NULL), height=8,width=10)
par(mfrow = c(1,1));
print('verboseScatterplot')

print('elements')
print(abs(geneModuleMembership[moduleGenes, column]))
print(abs(geneTraitSignificance[moduleGenes, 1]))

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste('Gene significance for',trait_column_of_interest, sep = " ", collapse = NULL),
                   main = paste("Module membership vs gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# names(datExpr)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


# names(datExpr)[moduleColors==args[4]]



# #todo
# #=====================================================================================
# #
# #  Code chunk 8
# #
# #=====================================================================================


# annot = read.csv(file = "GeneAnnotation.csv");
# dim(annot)
# names(annot)
# probes = names(datExpr)
# probes2annot = match(probes, annot$substanceBXH)
# # The following is the number or probes without annotation:
# sum(is.na(probes2annot))
# # Should return 0.


# #=====================================================================================
# #
# #  Code chunk 9
# #
# #=====================================================================================


# # Create the starting data frame
# # geneInfo0 = data.frame(substanceBXH = probes,
# #                       geneSymbol = annot$gene_symbol[probes2annot],
# #                       LocusLinkID = annot$LocusLinkID[probes2annot],
# #                       moduleColor = moduleColors,
# #                       geneTraitSignificance,
# #                       GSPvalue)
# geneInfo0 = data.frame(substanceBXH = probes,
#                       moduleColor = moduleColors,
#                       geneTraitSignificance,
#                       GSPvalue)
# # Order modules by their significance for weight
# modOrder = order(-abs(cor(MEs, weight, use = "p")));
# # Add module membership information in the chosen order
# for (mod in 1:ncol(geneModuleMembership))
# {
#   oldNames = names(geneInfo0)
#   geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
#                          MMPvalue[, modOrder[mod]]);
#   names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
#                        paste("p.MM.", modNames[modOrder[mod]], sep=""))
# }
# # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.trait_column_of_interest));
# geneInfo = geneInfo0[geneOrder, ]


# #=====================================================================================
# #
# #  Code chunk 10
# #
# #=====================================================================================


# write.csv(geneInfo, file = "geneInfo.csv")

