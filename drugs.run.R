#!/usr/bin/env Rscript
# help
# Rscript drugs.run.R 基因差异分析 药物与基因文件 结果目录
args=commandArgs(T)
source("/data2/wangb/tools/Fun_AD.r")

########################################################
## 1. 读取差异基因的结果文件，获取flodChange和padj信息
########################################################
DEG.fname = args[1] # 输入的DESEQ2差异基因结果
DEG.matrix <-  read.delim(DEG.fname, header=TRUE, quote = "", stringsAsFactors=FALSE, row.names = 1)
FC.vec <- DEG.matrix[, "log2FoldChange"]
FDR.vec <- DEG.matrix[, "padj"]


###############################################################
## 2. 读取药物基因调控信息，重新命名药物名称，获取药物关联的基因
###############################################################
drug.signature.fname = args[2] # 输入药物基因文件
Data.Fold <- args[3] # 输出的文件夹
if (length(args) > 3){
	SIG.THRESH <- args[4] # 显著性筛选阈值 
}else{
	SIG.THRESH <- 0.05
}

if (length(args) > 4){
	LOG2FOLDCHANGE.THRESH <- args[5] # 显著性筛选阈值 
}else{
	LOG2FOLDCHANGE.THRESH <- 0
}


drug.signature.string <- readLines(drug.signature.fname, -1)
drug.signature.list <- strsplit(drug.signature.string,'\t')
drug.name <- unlist(lapply(drug.signature.list, function(y)
	{
		return(paste0(y[2], ":", y[1]))
	}))
length(drug.name)
length(unique(drug.name))
drug.gene <- lapply(drug.signature.list, function(y)
	{
		 y.len <- length(y)
		 y.vec <- toupper(y[3:y.len])
		 y.vec <- unique(setdiff(y.vec, ""))
		 return(y.vec)
	})
names(drug.gene)<- drug.name
drug.all.fname <- paste0(Data.Fold, "drug.signature.txt")
write.geneset(drug.gene, drug.all.fname)


###############################################################
## 3. 获取药物导致的基因上调的信息 获取药物导致的基因下调信息
###############################################################
flag.vec <- as.character(sapply(drug.name, function(current.name)
	{
		current.vec <- unlist(strsplit(current.name, split = "-"))	
		return(current.vec[length(current.vec)])
	}))

up.ind <- which(flag.vec == "up")
down.ind <- which(flag.vec == "dn")

### up
drug.up.fname <- paste0(Data.Fold, "drug.signature.up.txt")
drug.up.gene <- drug.gene[up.ind]
write.geneset(drug.up.gene, drug.up.fname)

### down
drug.down.fname <- paste0(Data.Fold, "drug.signature.down.txt")
drug.down.gene <- drug.gene[down.ind]
write.geneset(drug.down.gene, drug.down.fname)

###############################################################
## 4. 读取差异基因结果文件，获取上下调显著基因。
###############################################################

deg_genes <- read.delim(DEG.fname, header = T, quote = "", row.names = 1,  stringsAsFactors=FALSE)

gene_symbols <- lapply(row.names(deg_genes), function(current.id)
{
	gene_id <- current.id
	y.vec <- toupper(unlist(strsplit(gene_id, split = "\\|"))[1])
	return(y.vec)
})
gene_symbols <- unique(setdiff(gene_symbols, NA))
g <- length(gene_symbols)
top_up = deg_genes[deg_genes$padj< SIG.THRESH & !is.na(deg_genes$padj) & deg_genes$log2FoldChange>LOG2FOLDCHANGE.THRESH,];
top_down = deg_genes[deg_genes$padj< SIG.THRESH & !is.na(deg_genes$padj) & deg_genes$log2FoldChange<LOG2FOLDCHANGE.THRESH,];

top_up_gene_id = rownames(top_up);
top_up_genes <- lapply(top_up_gene_id, function(y)
{
	y.vec <- toupper(unlist(strsplit(y, split = "\\|"))[1])
	return(y.vec)
})
top_up_genes <- setdiff(unique(top_up_genes), NA);

top_down_gene_id = rownames(top_down);
top_down_genes <- lapply(top_down_gene_id, function(y)
{
	y.vec <- toupper(unlist(strsplit(y, split = "\\|"))[1])
	return(y.vec)
})
top_down_genes <- setdiff(unique(top_down_genes), NA);


###############################################################
## 5. 基因富集分析，了解上下调基因富集的靶药物。
###############################################################
# 所有药物
overlap.infor <- t(sapply(names(drug.gene), function(current.drug)
{
	drug_genes = unique(toupper(setdiff(drug.gene[[current.drug]],"")));
	a = length(intersect(drug_genes,gene_symbols));
	b = length(intersect(drug_genes, top_up_genes));
	d = length(intersect(drug_genes, top_down_genes));
	
	up_sig_pvalue =  1- phyper(b-1,a,g-a, length(top_up_genes),lower.tail=T)
	dn_sig_pvalue =  1- phyper(d-1,a,g-a, length(top_down_genes),lower.tail=T)
	
	top.string <- paste(intersect(drug_genes, top_up_genes), collapse = " ")
	down.string <- paste(intersect(drug_genes, top_down_genes), collapse = " ")
	
	return(c(g, length(drug_genes), a, b, top.string, up_sig_pvalue, d, down.string, dn_sig_pvalue))
}))

colnames(overlap.infor) <- c("background_genes", "drug_genes", "overlap", "overlap_up", "overlap_up_genes", "up_pvalue", "overlap_dn", "overlap_dn_genes",	"dn_pvalue")			 
dn.pval.vec <- as.numeric(overlap.infor[,"dn_pvalue"])
dn.pval.order <- order(dn.pval.vec)
overlap.dn.order <- overlap.infor[dn.pval.order,]

write.table(overlap.dn.order, file= paste0(Data.Fold, "drug.all.DESeq.txt"),quote=F, col.names = NA, sep='\t');	

# 上调药物
overlap.up.infor <- t(sapply(names(drug.up.gene), function(current.drug)
{
	drug_genes = unique(toupper(setdiff(drug.up.gene[[current.drug]],"")));
	a = length(intersect(drug_genes,gene_symbols));
	b = length(intersect(drug_genes, top_up_genes));
	d = length(intersect(drug_genes, top_down_genes));
	
	up_sig_pvalue =  1- phyper(b-1,a,g-a, length(top_up_genes),lower.tail=T)
	dn_sig_pvalue =  1- phyper(d-1,a,g-a, length(top_down_genes),lower.tail=T)
	
	top.string <- paste(intersect(drug_genes, top_up_genes), collapse = " ")
	down.string <- paste(intersect(drug_genes, top_down_genes), collapse = " ")
	
	return(c(g, length(drug_genes), a, b, top.string, up_sig_pvalue, d, down.string, dn_sig_pvalue))
}))

colnames(overlap.up.infor) <- c("background_genes", "drug_genes", "overlap", "overlap_up", "overlap_up_genes", "up_pvalue", "overlap_dn", "overlap_dn_genes",	"dn_pvalue")
dn.pval.vec <- as.numeric(overlap.up.infor[,"dn_pvalue"])
dn.pval.order <- order(dn.pval.vec)
overlap.dn.order <- overlap.up.infor[dn.pval.order,]

write.table(overlap.dn.order, file= paste0(Data.Fold, "drug.up.DESeq.txt"),quote=F, col.names = NA, sep='\t');

# 下调药物
overlap.dn.infor <- t(sapply(names(drug.down.gene), function(current.drug)
{
	drug_genes = unique(toupper(setdiff(drug.down.gene[[current.drug]],"")));
	a = length(intersect(drug_genes,gene_symbols));
	b = length(intersect(drug_genes, top_up_genes));
	d = length(intersect(drug_genes, top_down_genes));
	
	up_sig_pvalue =  1- phyper(b-1,a,g-a, length(top_up_genes),lower.tail=T)
	dn_sig_pvalue =  1- phyper(d-1,a,g-a, length(top_down_genes),lower.tail=T)
	
	top.string <- paste(intersect(drug_genes, top_up_genes), collapse = " ")
	down.string <- paste(intersect(drug_genes, top_down_genes), collapse = " ")
	
	return(c(g, length(drug_genes), a, b, top.string, up_sig_pvalue, d, down.string, dn_sig_pvalue))
}))

colnames(overlap.dn.infor) <- c("background_genes", "drug_genes", "overlap", "overlap_up", "overlap_up_genes", "up_pvalue", "overlap_dn", "overlap_dn_genes",	"dn_pvalue")
dn.pval.vec <- as.numeric(overlap.dn.infor[,"dn_pvalue"])
dn.pval.order <- order(dn.pval.vec)
overlap.dn.order <- overlap.dn.infor[dn.pval.order,]

write.table(overlap.dn.order, file= paste0(Data.Fold, "drug.dn.DESeq.txt"),quote=F, col.names = NA, sep='\t');
