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
project = args[1]
condition = args[2]

load(file = paste(prefix, condition, 'WGCNA_networkConstruction_auto_02.RData', sep='_',  collapse = NULL))
load(file = paste(prefix, condition, '05.step6.WGCNA.RData', sep = "_", collapse = NULL))
load(file = paste(prefix, condition, '05.step8.WGCNA.TOM.RData', sep = "_", collapse = NULL))


exportNetworkToCytoscape(plotTOM, edgeFile='WGCNA_edges.txt', nodeFile='WGCNA_nodes.txt', weighted=T, threshold=0.02, nodeNames=colnames(datExpr))


nTop = 30;
IMConn = softConnectivity(datExpr[, colnames(datExpr)]);
top = (rank(-IMConn) <= nTop)
plotTOM_filtered <- plotTOM[top, top]

exportNetworkToCytoscape(plotTOM_filtered, edgeFile='WGCNA_edges_first30.txt', nodeFile='WGCNA_nodes_first30.txt', weighted=T, threshold=0.02, nodeNames=colnames(datExpr)[top])
