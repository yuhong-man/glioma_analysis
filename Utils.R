library(DESeq2)
library(ggplot2)

my_plotPCA <- function(dds,filename, nsub=1000){
	vsd <- vst(dds, blind=FALSE, nsub=nsub)
	pcaData <- plotPCA(vsd,intgroup=c("condition"),returnData=TRUE)
	percentVar <- round(100*attr(pcaData,"percentVar"))
	PcaPlotFig <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
	  geom_point(size=3) +
	  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
	  coord_fixed()
	ggsave(PcaPlotFig, filename = filename)
  	return(pcaData)
}
