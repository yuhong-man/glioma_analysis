#!/usr/bin/env Rscript
library(org.Hs.eg.db)
library(org.Mm.eg.db)

source("/data2/wangb/tools/Fun_AD.r")
source("/data2/wangb/tools/FunctionsForConvertingGeneSymbols.R")

args <- commandArgs(T)

path_input_network = args[1]
path_genes = args[2]

if(length(args) > 2){
	translate_gene_to_protein = as.logical(args[3])
}else{
	translate_gene_to_protein = T
}
if (length(args) > 3){
	translate_from = args[4]
}else{
	translate_from = 'SYMBOL'
}


if (length(args) > 4){
	orgdb = eval(parse(text=args[5]))
}else{
	orgdb = org.Hs.eg.db
}


network <- read.table(path_input_network, header=T, sep='\t')
genes <- scan(path_genes, what=' ', sep='\t')

print(head(network))

if(translate_gene_to_protein){
	# proteins = convertGeneAliasToSymbolAccordingToOrgDb(genes, fromType='ENSEMBL', toType='ENSEMBLPROT', use_from_type_if_no_match=F)
	proteins = convertGeneAliasToSymbolAccordingToOrgDb(genes, fromType=translate_from, toType='ENSEMBLPROT', use_from_type_if_no_match=F, OrgDb=orgdb)

}else{
	proteins = genes 
}

output_prefix = paste(path_genes, '.KDA', sep='')

df_proteins = data.frame(as.matrix(unlist(proteins)), stringsAsFactors=F)
df_proteins = df_proteins[!is.na(df_proteins),]
list_proteins = as.vector(df_proteins)

write.table(list_proteins, paste(path_genes, '.KDA.proteins', sep=''), sep='\t')


keyDriverAnalysis(network, list_proteins, outputfile=output_prefix)
