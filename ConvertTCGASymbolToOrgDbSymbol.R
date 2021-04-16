#!/data2/wangb/anaconda2/lib/R/bin/Rscript
  # [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
  # [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
  # [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
  # [25] "UNIGENE"      "UNIPROT"  

# If failed to map, use the original symbol. 

library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)

file_in <- 'D:/projects/primary_tracing/Manuscripts/random_forest_80_gene_and_GO/rf--max_features_auto--n_estimators_2000_first80_intersected_80.genes.txt'
# file_in <- 'D:/projects/primary_tracing/GEO_20190729/downloaded/Extracted_Required/GSE56087_RAW.FPKM_conf_hi.xls'
file_in <- 'Z:/rf--max_features_auto--n_estimators_2000_first77_intersected_77_delfrom80.genes_df'
file_in <- 'Z:/head_merged20'
file_out <- 'z:/out.xls'
fromType='ALIAS'

which_col_to_convert = 'gene'
remove_initial_symbol_col = F
keep_previous_col_name = T
multiVals = 'first'
args <- commandArgs(T)
if(length(args)>1){
  fromType <- args[1]
  toType <- args[2]
  file_in <- args[3]
  file_out <- args[4]
  if (toType == 'ALIAS'){
	multiVals = 'CharacterList'
  }
}
if(length(args)>=5){
  which_col_to_convert <- args[5]
  print('convert')
  print(which_col_to_convert)
}
if(length(args)>=6){
  remove_initial_symbol_col <- as.logical(args[6])
}
if(length(args)>=7){
  keep_previous_col_name <- as.logical(args[7])
}
if(length(args)>=8){
  orgdb = eval(parse(text = args[8]))
}else{
  orgdb = org.Hs.eg.db
}
print('using orgdb')
print(orgdb)

last <- function(x){x[[length(x)]]}

# table_to_convert <- read.table(file_in, header=F, stringsAsFactors =FALSE, row.names = 1)

table_to_convert <- read.table(file_in, header=1, sep='\t', stringsAsFactors=FALSE, encoding='UTF-8', check.names=F)
print('colnames(table_to_convert)')
print(colnames(table_to_convert))

convertGeneAliasToSymbolAccordingToOrgDb <- function(initial_symbols, fromType='ALIAS', toType='SYMBOL',OrgDb=org.Hs.eg.db, multiVals='first', SYMBOL_TO_UPPER=F){
  # Convert from alias to symbol. If failed, remain the previous. 
  print('initial sybmols: ')
  print(initial_symbols)
  if(SYMBOL_TO_UPPER){
    if (fromType == 'SYMBOL'){
      initial_symbols <- toupper(initial_symbols)
      }
  }
  symbols <- mapIds(OrgDb, keys=initial_symbols, keytype = fromType, column=toType, multiVals=multiVals)
  if (multiVals == 'CharacterList'){
	symbols <- apply(as.matrix(symbols), 1, function(x){
	s = paste0(x, collapse=', ')
	s = gsub(', NA', '', s)  
	})
  }
  for(i in seq(length(initial_symbols))){
    if(is.na(symbols[i])){
      symbols[i] <- initial_symbols[i]
    }
  }
  
  print('to: ')
  print(symbols)
  
  
  return(symbols)
  
  
  # return(df_alias_symbol$SYMBOL)
}

# xs <- lapply(table_to_convert$fake_symbol, fromSingleSymbolToFormalSymbol)
# table_to_convert[,'symbol'] <- unlist(xs)
 
table_to_convert[, toType] <- convertGeneAliasToSymbolAccordingToOrgDb(table_to_convert[,which_col_to_convert], fromType=fromType, toType=toType, multiVals = multiVals, OrgDb=orgdb)
# Remove the first column and move the last column (symbol) to first. 



if(remove_initial_symbol_col){
  print('remove_initial_symbol_col')
  table_to_convert[[which_col_to_convert]] <- NULL
}

print('head(table_to_convert)')
print(head(table_to_convert))
original_col_names <- colnames(table_to_convert)
print('original_col_names')
print(original_col_names)
# print(original_col_names[2:length(original_col_names)])
# print(original_col_names[length(original_col_names)])
# print(original_col_names[1])
# print(original_col_names[3:length(original_col_names)-1])
modified_col_names <- c(original_col_names[length(original_col_names)], original_col_names[1:length(original_col_names)-1])
print('modified_col_names')
print(modified_col_names)
# modified_col_names <- original_col_names

df_out <- table_to_convert[,modified_col_names]
if(keep_previous_col_name){
  colnames(df_out)[1] <- which_col_to_convert
}

file_out <- file(file_out, "w", encoding="UTF-8") 
write.table(df_out, file_out, sep='\t', quote=F, row.names = F)

