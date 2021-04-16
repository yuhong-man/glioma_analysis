#!/data2/wangb/anaconda2/bin/Rscript
args=commandArgs(T)
library(DESeq2)
load(file = args[2]);
# args[1] 为项目名称 如LUSC,BRCA
# args[2] 为上一步的缓存文件 04.step1.DESeq2.RData | 04.step2.DESeq2.RData


#############################################################################################
### DESeq2 数据处理
#############################################################################################
dds <- DESeqDataSetFromMatrix(countData = inputData, colData, design = ~ condition)
keep <- rowSums(counts(dds)) >= 20 # 对counts数进行初步的过滤
dds <- dds[keep,]
dds <- DESeq(dds) #三步走：文库大小估计；离散程度估计；统计检验
res <- results(dds,name= "condition_T_vs_N") #结果呈现
res <- res[order(res$padj),]
resDF = as.data.frame(res)
resDF$gene_id = rownames(resDF)
resDF <- resDF[,c(7,1,2,3,4,5,6)]
# 输出结果文件
write.table(resDF,file=paste(args[1],'DESeq2_DEG.list', sep = "_", collapse = NULL),sep="\t",quote= FALSE,row.names = FALSE)




#############################################################################################
### inputData colData 变量保存
#############################################################################################
inputData <- as.data.frame(counts(dds, normalized=TRUE)) ##表达量归一化
save(inputData,colData, file=paste(args[1],'DESeq2_inputData_colData.RData', sep = "_", collapse = NULL))




#############################################################################################
### MA 统计图
#############################################################################################
pdf(paste(args[1],'DESeq2_MA.pdf', sep = "_", collapse = NULL))
plotMA(res,ylim=c(-3,3))
dev.off()