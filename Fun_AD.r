### AUXILAEY FUNCTIONS FOR THE PROJECT AGING DISEASE

######################################################################################################
#### INVERSE.NORM INVERSELY NORMALIZES THE GENE BY MAPPING TO A STANDARD NORMAL DISTRIBUTION 
####              ACCORDING TO RANK IN SMAPLES
######################################################################################################
	
    inverse.norm <- function(expr.matrix, method = "min")
	{  
		#### inverse.norm inversely normalizes the gene by mapping to a standard normal distribution according to rank in smaples
		#### Author: Jialiang Yang  
		#### Date: 6/2/2015 version: 1
		##################################################################
		### Input: 
		###       expr.matrix: matrix to be normalized, format: row-gene colum-sample
		###       method: method to deal with ties, c("average", "first", "random", "max", "min")
		### Output:
        ###       matrix.norm: the matrix after normalization		
		###***************************************************************** 

	    sample.size <- dim(expr.matrix)[2]
	    print(sample.size)  ##82
		
		ref <-  qnorm(seq(0, 1, length.out = (sample.size+2)), mean = 0, sd =1 )
		ref <- ref[2:(sample.size+1)]
		
		cat('ref mean: ',mean(ref),' sd: ',sd(ref),'\n',sep='')   ### ref mean: 1.265925e-16 sd: 0.9598929
		print(shapiro.test(ref))    ### W = 0.9973, p-value = 0.9999
	   
	    ### reversely map the samples in a gene to ref accoridng to its rank
		   normalize.matrix <- t(apply(expr.matrix,1,function(x){ref[rank(x,ties.method = method)]}))
	}
	
######################################################################################################
#### LOG2TRANS PERFORMS THE LOG2 TRANSFORMATION OF GENE EXPRESSION PROFILE
######################################################################################################

    log2trans <- function(dat.matrix, add.value = 1)
	{  
		#### log2trans performs the log2 transformation of gene expression profile
		####           and replaces 0 with global min/ 10
		#### Author: Jialiang Yang  
		#### Date: 4/6/2016 version: 1
		##################################################################
		### Input:
		###        dat.matrix: row-gene col-sample
		### Output: 
		###        dat.log2: the log2 transformed data
		###***************************************************************** 
		    dat.matrix <- dat.matrix + add.value
			
			dat.log2 <- log2(dat.matrix)
			
		    return(dat.log2)
	}

######################################################################################################
#### PLOT.PCA PLOTS THE PCA GRAPH OF THE DATA MATRIX
######################################################################################################

	plot.pca <- function(dat.matrix, figure.fname, color.vec = NULL, label = FALSE, PC.scale = FALSE)
	{  
		#### plot.pca plots the pca graph of the dat.matrix
		#### Author: Jialiang Yang  
		#### Date: 3/18/2016 version: 1
		##################################################################
		### Input: 
		###       dat.matrix: format: row-gene colum-sample
		###     figure.fname: a pdf file to output the figure
		###       color.vec:  must be the same size as col(expr.matrix) to specify color for each sample
		###            labe: decide to print the node label alongside the node or not
		###***************************************************************** 
		
			library("preprocessCore")
			library(ggplot2)
			
			pca <- prcomp(t(dat.matrix), scale.= PC.scale)
			
			### construct the data to plot
				if(length(color.vec) == 0)
				{
					data.matrix <- data.frame(pca$x[,1], pca$x[,2])
					names(data.matrix) <- c("PC1", "PC2")
					
					P <- ggplot(data.matrix, aes(x = PC1, y = PC2)) + geom_point(size = 3)+ 
					         scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
							 labs(title =  paste0("PCA Plot"), x= "PC1", y = "PC2")+
							 theme(plot.title = element_text(size = 18),axis.title=element_text(size=15),axis.text = element_text(size = 15))
					
					ggsave(plot=P,height=7,width=7,dpi=300, filename=figure.fname, useDingbats=FALSE)	
					
				}else{
					try(if(length(color.vec) != dim(dat.matrix)[2]) stop("color.vec and expression matrix are not of the same sample size!"))
					data.matrix <- data.frame(pca$x[,1], pca$x[,2], color.vec)
					names(data.matrix) <- c("PC1", "PC2", "Group")
					
					if (label){
					
						P <- ggplot(data.matrix, aes(x = PC1, y = PC2,color = Group)) + geom_point(size = 15)+ 
						    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
							geom_text(aes(label=Group,x = PC1,y = PC2))+
							labs(title =  paste0("PCA Plot"), x= "PC1", y = "PC2")+
							theme(plot.title = element_text(size = 18),axis.title=element_text(size=15),axis.text = element_text(size = 15))
						
						ggsave(plot=P,height=7,width=7,dpi=300, filename=figure.fname, useDingbats=FALSE)	
						
					}else{
					    P <- ggplot(data.matrix, aes(x = PC1, y = PC2,color = Group)) + geom_point(size = 3)+ 
					   	 scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
							labs(title =  paste0("PCA Plot"), x= "PC1", y = "PC2")+
							theme(plot.title = element_text(size = 18),axis.title=element_text(size=15),axis.text = element_text(size = 15))
						
						ggsave(plot=P,height=7,width=7,dpi=300, filename=figure.fname, useDingbats=FALSE)
					}
				}
	}
	
######################################################################################################
### PLOT.PCA.GROUP PLOTS THE PCA GRAPH FOR A SET OF GENE EXPRESSION PROFILES
######################################################################################################

    plot.pca.group <- function(group.fname.list, group.name.list, figure.fname)
	{  
		#### plot.pca.group plots the pca graph for a set of gene expression profiles
		#### Author: Jialiang Yang  
		#### Date: 3/18/2016 version: 1
		##################################################################
		### Input: 
		###       group.fname.list: a string vector to store the path of each expression files
        ###		              Format: expression matrix (\t deliminated), row-gene colum-sample
		###       group.name.list: group names foe each file in group.fname.list
        ###       figure.fname: the figure file name		
		###***************************************************************** 
	    	library(ggplot2)

		    try(if(length(group.fname.list) != length(group.name.list)) stop("group.fname.list and group.name.list are not consistent"))
			
			### (1) Identify overlapping genes
			### ----------------------------------------------
				current.fname <- group.fname.list[1]
				current.matrix <- read.delim(current.fname, header=TRUE, stringsAsFactors=FALSE, row.names = 1)
				merge.gene.ID <- rownames(current.matrix) 
				
				for (i in 2:length(group.fname.list))
				{
					current.fname <- group.fname.list[i]
					current.matrix <- read.delim(current.fname, header=TRUE, stringsAsFactors=FALSE, row.names = 1)
					current.gene.ID <- rownames(current.matrix)
					
					merge.gene.ID <- intersect(merge.gene.ID, current.gene.ID)
				}
				
				length(merge.gene.ID)  ### 17378
			
			### (2) Merge the gene expression profile for all the tissues
			### ----------------------------------------------
			    current.group <- group.name.list[1]
				current.fname <- group.fname.list[1]
				current.matrix <- read.delim(current.fname, header=TRUE, stringsAsFactors=FALSE, row.names = 1)
				merge.expr.matrix <- current.matrix[merge.gene.ID,]
				group.full.vec <- rep(current.group, dim(current.matrix)[2])
				
				for (i in 2:length(group.fname.list))
				{
					current.group <- group.name.list[i]
					current.fname <- group.fname.list[i]
					current.matrix <- read.delim(current.fname, header=TRUE, stringsAsFactors=FALSE, row.names = 1)
					merge.expr.matrix <- cbind(merge.expr.matrix, current.matrix[merge.gene.ID,])
					
					group.full.vec <- c(group.full.vec, rep(current.group, dim(current.matrix)[2]))
				}
				
				merge.expr.matrix[1:5,1:5]
				dim(merge.expr.matrix)          ###  17378  2651
				length(group.full.vec)         ### 2651
				
			### (3) Plot the graph
			### ----------------------------------------------
			    plot.pca(merge.expr.matrix, figure.fname, group.full.vec)
	}
	
######################################################################################################
#### read.VCF reads the SNP file in vcf format 
######################################################################################################
    
	read.VCF.genotype <- function(VCF.fname, FILTER.criteria = "PASS", RD.THRESH = 10)
	{
	    #### read.VCF reads the SNP information called by GATK
		##################################################################
		### Input:
		###       VCF.fname: vcf file with the first 123 lines annotations and the following lines with format 
		###                CHROM	POS	   ID	REF	ALT	QUAL	FILTER	   INFO	  FORMAT	sample
	    ###                 chr1	14574	.	A	G	49.77	SnpCluster	AC=1;AF=0.500;AN=2;BaseQRankSum=0.736;ClippingRankSum=0.736;DP=3;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=60.00;
		### Output:
		###       SNP.chrpos.matrix: a matrix with the first column bing chr and pos and the second column genotype
		###*****************************************************************
		
			library(VariantAnnotation)
			
			VCF.matrix <- readVcf(VCF.fname, "hg19")
			FILTER.vec <- rowRanges(VCF.matrix)$FILTER
			FILTER.ind <- which(FILTER.vec == "PASS")
			RD.vec <- info(VCF.matrix)$DP
			RD.ind <- which(RD.vec >= RD.THRESH)
			KEEP.ind <- intersect(FILTER.ind, RD.ind)
			
			GT.vec <- geno(VCF.matrix)$GT[KEEP.ind,]
			
			### convert GT.vec into regular format, e.g. C/T 
				VCF.locus.str <- names(GT.vec)
				
				VCF.locus.matrix <- t(sapply(VCF.locus.str, function(current.str)
				{
					str.vec <- unlist(strsplit(current.str, split = "_"))
					
					if(length(str.vec) > 2)
					{
					   str1.vec <- paste0(str.vec[1: (length(str.vec)-1)], collapse = "_")
					   
					   str.vec <- c(str1.vec, str.vec[length(str.vec)])
					}
					
					return(str.vec)
				}))
				
            GT.matrix <- cbind(VCF.locus.matrix, as.character(GT.vec))
			
			### construct genotype matrix
				genotype.matrix <- t(apply(GT.matrix, 1, function(current.locus)
				{
					if(current.locus[3] == "0/1"){current.GT <- current.locus[2]}else{
					   locus.vec <- unlist(strsplit(current.locus[2], split = "/"))
					   current.GT <- paste0(locus.vec[2], "/", locus.vec[2])
					} 
					
					return(c(current.locus[1], current.GT))
				}))
				
		return(genotype.matrix)		
	}
	
######################################################################################################
#### CONSTRUCT.COEXPR CONSTRUCTS THE CO-EXPRESSION NETWORK
######################################################################################################

    ###----------------------------------------------------------------------------------------------
	#### Programs for updated WGCNA
	###----------------------------------------------------------------------------------------------

		plotTOMHeatmap2 <- function(coexppClusters, geneModules, samplingThreshold = 1000,...)
		{
			#input are calculated from coexpressionAnalysis in the coexpp pkg
				if (length(coexppClusters) > samplingThreshold) {
					warning("Too many observations to construct heatmap directly. Sampling to tractable number of observations.")
					samples <- clusters(coexppClusters)$order[sort(sample(length(coexppClusters),samplingThreshold))]
					m=sampleMatrix(coexppClusters, samples, kind = "tom")
					Colors=geneModules[samples]
				}else {
					samples <- clusters(coexppClusters)$order
					m=tom(coexppClusters)[samples, samples]
					Colors=geneModules[samples]
				}
				diag(m) = NA
				Colors=as.character(Colors)
				heatmap(m, Rowv=NA, Colv=NA, scale="none", revC=TRUE, symm=TRUE, labRow="", labCol="", ColSideColors=Colors, RowSideColors=Colors, ...)
		}
		
		plotRsq=function(x,threshold=0.9){
			x=x[floor(x[,1])==x[,1],]
			plot(x[,1],-sign(x[,3])*x[,2],main='Scale Independence',xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit, sign R^2', type='n')
			text(x[,1],-sign(x[,3])*x[,2],labels=x[,1],cex=0.9,col='red')
			abline(h=threshold,col='red')
		}

    ###----------------------------------------------------------------------------------------------
	#### construct.coexpr
	###----------------------------------------------------------------------------------------------
	
		construct.coexpr <- function(expr.fname, out.fold, sft.Power = NULL)
		{  
			#### construct.coexpr constructs the co-expression network
			#### works for R version 3.0.3. might have problem in higher version 
			#### Author: Jialiang Yang  
			#### Date: 3/20/2016 version: 1
			##################################################################
			### Input: 
			###       expr.fname: the full path of the expression file 
			###                   format: row-gene, column-sample, separated by tab "\t"
			###       out.fold: a fold to store the TOM file, scale free plot and dendegram	
			###       sft.Power: the soft power to be set, by default NULL guess by the program
			###***************************************************************** 

				if(!file.exists(out.fold)){dir.create(out.fold)}
				
				if (length(sft.Power) == 0){out1.fold <- paste0(out.fold, "NULL/") }else{out1.fold <- paste0(out.fold,sft.Power,"/")}
				if(!file.exists(out1.fold)){dir.create(out1.fold)}
				
				out.fold <- out1.fold
				
				expr.matrix <- read.delim(expr.fname,as.is=TRUE, row.names = 1, header=TRUE)
				
				gene.name <- rownames(expr.matrix) 
				
				### Read in gene expression matrix; samples in rows and genes in columns.
				###--------------------------------------------------------------------------
					dat <- as.matrix(t(expr.matrix))
				
					### Specify the soft power coefficient; if NULL, the power coefficient will be computed automatically.
						sftPower=sft.Power     
						
					### Set the number of CPU cores to be recruited.
						ncores=10
					
				### Calculate modules
				###--------------------------------------------------------------------------
					library(coexpp)
					library(flashClust)
					library(WGCNA)
					coexppSetThreads(ncores)
					obj=coexpressionAnalysis(dat,beta=sftPower)
					disableWGCNAThreads()

				### Output modules
				###--------------------------------------------------------------------------
					modules=data.frame(gene=obj$clusters@clusters$labels,obj$intraModularStatistics,module=obj$geneModules)
					print(dim(modules))
					print(modules[1:5,1:5])
					
					write.csv(modules,file=paste0(out.fold,"modules_Normalize.csv"),row.names=FALSE)

				### Plot scale free fit
				###--------------------------------------------------------------------------
					pdf(paste0(out.fold,"scaleFree_Normalize.pdf"),width=14,height=14)
						plotScaleFree(obj$clusters)
					dev.off()
					
				###plot soft power figure
				###--------------------------------------------------------------------------
				    if (is.null(sft.Power))
					{
						pdf(paste0(out.fold,"softPower.pdf"),width=14,height=14)
						plotRsq(obj$clusters@sftStatistics,0.9)
						dev.off()
					}
		
				### Plot cluster dendrogram
				###--------------------------------------------------------------------------
					pdf(paste0(out.fold, "den.cluster_Normalize.pdf"),width=14,height=14)
					plotClustering(obj$clusters, obj$geneModules,'Dynamic Tree Cut',
						dendroLabels=FALSE,hang=0.03,addGuide=TRUE,guideHang=0.05,
						main='Gene dendrogram and module colors')
					dev.off()
					
				### Plot TOM matrix
				###--------------------------------------------------------------------------
					pdf(paste0(out.fold, "den.tom_Normalize.pdf"),width=14,height=14)
					plotTOMHeatmap2(obj$clusters, obj$geneModules)
					dev.off()
					
				### TOM Matrix
				###--------------------------------------------------------------------------
				   if(is.null(sft.Power)){sftPower = slot(obj$clusters, "beta")}else{sftPower = sft.Power}
				   TOM = TOMsimilarityFromExpr(dat, power = sftPower);
				   rownames(TOM) <- gene.name
				   colnames(TOM) <- gene.name
				   
				   save(TOM, file = paste0(out.fold,"TOM.RData"))		  
		}
	
######################################################################################################
#### GET.MODULE.SYM RETRIEVES THE GENE SYBMOLS IN EACH MODULE IN MODULE.FNAME GENERATED BY WGCNA
######################################################################################################

	get.module.sym <- function(module.fname, out.fold, genesym.fname)
	{  
		#### get.module.sym retrieves the gene sybmols in each module in module.fname
		#### Author: Jialiang Yang  
		#### Date: 9/28/2016 version: 1
		##################################################################
		### Input:
		###        module.fname: module file generated by WGCNA
		###        out.fold: the output file folder
		###        genesym.fname: a file to correpond Ensemble ID and gene sysmbol
		###***************************************************************** 

			module.matrix <- read.csv(module.fname, header=TRUE, stringsAsFactors=FALSE, row.names = 1)
			module.matrix[1:5,]
			
			gene.sym.matrix <- read.delim(genesym.fname, header=TRUE, stringsAsFactors=FALSE, row.names = 1)
			gene.sym.matrix[1:5,]
			gene.ID <- gene.sym.matrix[, "GeneID"]
			gene.symbol <- gene.sym.matrix[, "GeneSymbol"]
			
			module.vec <- as.character(module.matrix[, "module"])
			Ensemble.ID.vec <- as.character(rownames(module.matrix))
			
			module.unique <- unique(module.vec)
			
			### construct symbol txt file for each module
			sapply(module.unique, function(current.module)
			{
				current.ind <- which(module.vec == current.module)
				current.Ensemble <- Ensemble.ID.vec[current.ind]
				
				sym.ind <- match(current.Ensemble, gene.ID)
				
				module.sym <- gene.symbol[sym.ind]
				
				module.out.fname <- paste0(out.fold, current.module)
				
				write(module.sym, module.out.fname)
			})
	}
	
	
######################################################################################################
#### MDC Module Differential Connectivity analysis
######################################################################################################

	MDC <- function(work.fold, inputfnameB, inputfnameBModule, inputfname, shortnames, headCol2, headCol, corrlpower)
	{  
		#### get.module.sym retrieves the gene sybmols in each module in module.fname
		#### Author: Bin Zhang adpated by Jialiang Yang 
		#### Date: 9/28/2016 version: 1
		##################################################################
		### #    Input:
		#         1) inputfnameB       -- mRNA data (tab-delimited text file) for one disease state B
		#         2) inputfnameBModule -- module assignment file as output from WINA or WGCNA
		#         3) inputfname        -- mRNA data (tab-delimited text file) for another disease state A
		#         4) shortnames        -- shortnames for B and A, respectively. Order is important
		#         5) headCol           -- # of gene information columns in inputfname
		#         6) headCol2          -- # of gene information columns in inputfnameB
		#         7) corrlpower        -- the exponent of the power function, determined by WINA or WGCNA
		#
		#   Output:
		#         1) "*_wFDR.xls"      -- MDC with false discovery rates
		#         2) "*.png"           -- MDC plot
		###***************************************************************** 

		setwd(work.fold)		
		
		library(lattice) # require is design for use inside functions 

		inputfnameB       = inputfnameB; headCol2 = headCol2;
		inputfnameBModule = inputfnameBModule
		inputfname        = inputfname; headCol =headCol;
		shortnames        = shortnames
		corrlpower        = corrlpower

		#
		# -----------------------------End of Parameters to be changed --------------------------------------

		imgwid=600; imghei=400

		# specify the directory for holding analysis results, you need make such a sub-directory 
		#    under your working directory if it doesn't exist
		outputDir1  ="DiffConnectvty/"
		dir.create(outputDir1)

		outputDir  ="DiffConnectvty/Random/"
		dir.create(outputDir)

		# image type
		imagetype="png" #"pdf"


		fname      =getFileNameNopath(inputfname)
		fnameB     =getFileNameNopath(inputfnameBModule)
		extname    =getFileExtension(inputfname)

		fname = paste(fnameB, "_VS_", shortnames[2], "_MDC", sep="")
		flog  = paste(outputDir1, fname, "_Randoms", ".xls",       sep='')
		flog2 = paste(outputDir1, fname, "_wFDR", ".xls",       sep='')
		fimgRt= paste(outputDir1, fname, ".png",sep='')  #png only

		  #------------------------- second file --------------------------------------
		  allMatrix2 <- read.delim(inputfnameB,sep="\t", header=T)
		  dim(allMatrix2)
		  genesInfor2 <- cbind( allMatrix2[,c(1:headCol2)])
		  rowTitles2 <- as.character(genesInfor2[,1])

		  #These are the expression values
		  datExpr2 <-t(allMatrix2[,-c(1:headCol2)])
		  dim(datExpr2)

		  no.genes <- dim(datExpr2)[2]

		  # ------ modules from 2nd file --------------------------------
		  allMatrix2modinfo <- read.delim(inputfnameBModule,sep="\t", header=T)
		  dim(allMatrix2modinfo)
		  ncols2m = dim(allMatrix2modinfo)[2]

		  allMatrix2modinfo <- as.matrix(allMatrix2modinfo)[,c(1,ncols2m)]
		  genes2Idx = cbind(as.matrix(genesInfor2)[,1], c(1:no.genes) )
		  merged = merge(genes2Idx, allMatrix2modinfo, by.x=1, by.y=1, all.x=T)
		  merged = as.matrix(merged)
		  dim(merged)
		  morder = order(as.integer(merged[,2]))
		  allMatrix2modinfoXX = merged[morder,] # aligned to the other gene infos
		  modulescolorB = as.character(merged[morder,3])
		  modulescolorB = ifelse(is.na(modulescolorB), "grey", modulescolorB)

		  modtb = table(modulescolorB)
		  modulenames = names(modtb)
		  modtb = modtb[modulenames!="grey"]

		  # by modules  size
		  #mo = order(-as.integer(modtb))
		  #modtb = modtb[mo]
		  modulenames = names(modtb)

		  umodulesizes = union(sort(as.integer(modtb)), NULL)

		  no.modules  = length(modulenames)
		  no.nets     = 2

		  #------------------------- first file --------------------------------------
		  #
		  allMatrix <- read.delim(inputfname,sep="\t", header=T)
		  dim(allMatrix)
		  gIDa = as.character(as.matrix(allMatrix[,1]))
		  
		  # align two datasets
		  matchedIdx = getMatchedIndex2way(as.character(rowTitles2), gIDa)
		  allMatrix= allMatrix[matchedIdx[,2],]

		  genesInfor <- cbind(allMatrix[,c(1:headCol)])
		  rowTitles <- as.character(genesInfor[,1])

		  #These are the expression values
		  datExpr <-t(allMatrix[,-c(1:headCol)])
		  dim(datExpr)

		  #corhelp<- abs(corhelp)

		#*-------------------------------------------------------------------------------------
		#* initilization

		uno = length(umodulesizes)

		meanPermodule   = matrix(0, no.modules, no.nets)
		sdPermodule     = matrix(0, no.modules, no.nets)
		ks.pvalues      = rep(1, no.modules)
		linkCor         = matrix(-1, no.modules, no.nets)

		meanPermoduleTop   = matrix(0, no.modules, no.nets)
		sdPermoduleTop     = matrix(0, no.modules, no.nets)
		ks.pvaluesTop      = rep(1, no.modules)

		#############################################################################
		#-------------------- permute samples ---------------------------------------
		#
		no.perms = 50
		meanFoldChange = NULL
		for( y in c(1:(no.perms+1)) ) {

		  print(paste("********* random samples ", y, "**************"))

		  if(y!=(no.perms+1) ){ # permuated values
			 XdatExpr =apply(datExpr,  2, permuteVect)
			 XdatExpr2=apply(datExpr2, 2, permuteVect)

		  } else { # true value
			 XdatExpr =datExpr
			 XdatExpr2=datExpr2
		  }

		  for ( x in c(1:no.modules) ) {

			  print(paste(x, "/", no.modules, ":", modulenames[x], sep="") )

			  #xsel = sample(c(1:no.genes), umodulesizes[x], replace=F)
			  xsel = modulescolorB == modulenames[x]

			  corhelp <- cor(XdatExpr[, xsel],  use = "pairwise.complete.obs")
			  diag(corhelp) <- 0
			  corhelp = corhelp^corrlpower
			  links   <- apply(abs(corhelp), 1, sum, na.rm=T)

			  lpanel  <- lower.tri(corhelp)
			  corhelp <- abs(corhelp[lpanel])
			  isel    <- !is.na(corhelp)
			  no.corrls = length(corhelp)

			  corhelp <- corhelp[isel]
			  
			  # ---------------------------------------------------------------
			  corhelp2 <- cor(XdatExpr2[, xsel], use = "pairwise.complete.obs") 
			  corhelp2 <- corhelp2^corrlpower
			  diag(corhelp2) <- 0
			  links2   <- apply(abs(corhelp2), 1, sum, na.rm=T)
			  
			  corhelp2 <- abs(corhelp2[lpanel])
			  isel     <- !is.na(corhelp2)
			  corhelp2 <- corhelp2[isel]
			  
			meanPermodule[x, 1] = mean(corhelp2, na.rm=T)
			meanPermodule[x, 2] = mean(corhelp, na.rm=T)

			rm(corhelp2)
			rm(corhelp)
			collect_garbage()
		  }

		  rm(XdatExpr)
		  rm(XdatExpr2)
		  collect_garbage()

		  yfold = meanPermodule[,1]/meanPermodule[,2]
		  if(y!=(no.perms+1) ){ # permuated values
			meanFoldChange = cbind(meanFoldChange, yfold)
		  }

		  print(yfold)
		}

		xfinal = cbind(modulenames, meanFoldChange)
		colnames(xfinal) <- c("module", paste("MDC_random_samples_", c(1:no.perms), sep="_") )
		xtitle= colnames(xfinal)
		#write.table(xfinal, flog, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

		#############################################################################
		#-------------------- permute genes ---------------------------------------
		#
		no.perms = 50
		GlobalyRandomMDC = NULL
		for( y in c(1:no.perms) ) {
		  print(paste("********* random genes", y, "**************"))

		  GRmeanPermodule   = matrix(0, no.modules, no.nets)
		  for ( x in c(1:no.modules) ) {

			  #print(paste("random genes", umodulesizes[x], "/", uno, sep="") )
			  #xsel = sample(c(1:no.genes), umodulesizes[x], replace=F)

			  print(paste("randomize genes:", modulenames[x], " ", x, "/", no.modules, sep="") )
			  xsel = sample(c(1:no.genes), as.integer(modtb[x]), replace=F)

			  corhelp <- cor(datExpr[, xsel],  use = "pairwise.complete.obs")
			  diag(corhelp) <- 0
			  corhelp = corhelp^corrlpower
			  links   <- apply(abs(corhelp), 1, sum, na.rm=T)

			  lpanel  <- lower.tri(corhelp)
			  corhelp <- abs(corhelp[lpanel])
			  isel    <- !is.na(corhelp)
			  no.corrls = length(corhelp)

			  corhelp <- corhelp[isel]
			  
			  # ---------------------------------------------------------------
			  corhelp2 <- cor(datExpr2[, xsel], use = "pairwise.complete.obs") 
			  corhelp2 <- corhelp2^corrlpower
			  diag(corhelp2) <- 0
			  links2   <- apply(abs(corhelp2), 1, sum, na.rm=T)
			  
			  corhelp2 <- abs(corhelp2[lpanel])
			  isel     <- !is.na(corhelp2)
			  corhelp2 <- corhelp2[isel]
			  
			  GRmeanPermodule[x, 1] = mean(corhelp2, na.rm=T)
			  GRmeanPermodule[x, 2] = mean(corhelp, na.rm=T)

			  rm(corhelp2)
			  rm(corhelp)
			  collect_garbage()
		  }
		  GlobalyRandomMDC = cbind(GlobalyRandomMDC, GRmeanPermodule[,1]/GRmeanPermodule[,2])
		}

		GRMDC_mean = apply(GlobalyRandomMDC, 1, mean)
		GRMDC_sd = apply(GlobalyRandomMDC, 1, sd)
		GR_Dat = cbind(yfold, GRMDC_mean, GRMDC_sd)
		GR_FDR_N = apply(GR_Dat, 1, MDC_FDR_by_normal_distr)

		GR_FDR = apply(cbind(yfold, GlobalyRandomMDC), 1, MDC_FDR)

		xfinal = cbind(modulenames, meanFoldChange, GlobalyRandomMDC)
		colnames(xfinal) <- c("module", paste("MDC_random_samples_", c(1:no.perms), sep=""), 
							  paste("MDC_random_genes_", c(1:no.perms), sep="") )

		write.table(xfinal, flog, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

		#--------------- compute FDR --------------------
		#
		rmean = apply(meanFoldChange, 1, mean)
		rsd   = apply(meanFoldChange, 1, sd)
		mdcDat = cbind(yfold, rmean, rsd)
		mdcFDR_N = apply(mdcDat, 1, MDC_FDR_by_normal_distr) # FDR by distribution

		mdcFDR = apply(cbind(yfold, meanFoldChange), 1, MDC_FDR)

		comFDR = apply(cbind(GR_FDR, mdcFDR), 1, max)

		final = cbind(modulenames, round(yfold,2),  signif(comFDR,2), signif(mdcFDR,2), signif(GR_FDR,2))
		colnames(final) = c("module", "MDC",  "FDR", "FDR_random_samples", "FDR_random_genes")
		final

		write.table(final, flog2, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)


		#################################################################################
		#---------------------- plot ratio of means -----------------------------------
		#
		rylab = paste("MDC: k_", shortnames[1]," / k_", shortnames[2], sep="")

		imgwid=1400; imghei=600;
		ignore_small_modules=F; module_sizecut=100

		hlines = c(2,1.5,1)
		if (max(yfold)>20){ hlines = c(20, 10, 5) }

		rylab = paste("MDC(", shortnames[1],", ", shortnames[2], ")", sep="")
		od = order(-yfold)
		dispDat = cbind(yfold[od]); dispNames=modulenames[od];
		rownames(dispDat) <- dispNames; dispModuleSize = modtb[od]
		dispSign = ifelse(comFDR[od]<= 0.10, "*", "")

		yymax = closest_integer_bound(max(yfold))
		barplot_by_oneGroup(datMatrix=dispDat, maskMatrix=cbind(dispSign), 
						  rowcol=dispNames, vlines=hlines,
						  title="", #"Modular Differential Connectivity", 
						  ylab=rylab, 
						  imgfile=fimgRt, imgwid=2000, imghei=800, 
						  mcex=0.5, xlab_cex=0.65, xaxisLabel_rotate=45, 
						  legendcols=1, legendxDelta=0, legendCex=1, 
						  ipointsize=12, iunits="px", ires=300, icompression="lzw",
						  show_value=FALSE, showLegend=FALSE, show_xaxisLabel=TRUE, show_mask=TRUE,
						  xmargin=3, ymargin=5, zmargin=0.5, xlabel_orientation=0, ValueXdelta=0, 
						  ValueYdelta=2*yymax/100, valuecex=1, yMAX=yymax, yMIN=0)
	}
	
	# to split "abc|123", use sep="\\|", "abc.123" use "\\."
	splitString =function(mystring, separator="; "){
	  splitted = NULL
	  for (each in mystring){
		 if (is.na(each) | is.null(each)){
			next
		 }
		 a=unlist( strsplit(each, separator) )
		 splitted =c(splitted, a)
	  }
	  #a=unlist( strsplit(mystring, separator) )
	  return(splitted )
	}


	#get the filename without extension
	#
	getFileExtension=function(fullfname){
		splitted=unlist( strsplit(fullfname, "\\.") )
		
		if( length(splitted) >1){
		  return (splitted[length(splitted)])
		} else{
		  return ("")
		}
	}

	#get the filename without extension
	getFileName=function(fullfname){
		ext=getFileExtension(fullfname)
		if(ext ==""){
		   return (fullfname)
		}
		extd = paste(".", ext, sep="")
		splitted=splitString(fullfname, extd)

		splitted[1]
	}

	#get the filename without extension
	getFileNames=function(fullfnames){

	  final = NULL
	  for(ef in fullfnames) {
		 fn = getFileName(ef)
		 final = c(final, fn)
	  }
	  return (final)
	}

	#get the filename without extension and path information
	getFileNameNopath=function(fullfname){
	   res = NULL
	   for(each in fullfname) {
		myfilename = getFileName(each)
		splitted=unlist( strsplit(myfilename, "/") )
		 res= c(res, splitted[length(splitted) ])
	   }
	   return (res)
	}

	# assume that find indices of the common components in the two sets
	#
	getMatchedIndex2way=function(cvector, dvect){
	  fullindex = c(1:length(cvector) )
	  orgIdx    = cbind(cvector, fullindex)

	  index2    = c(1:length(dvect))
	  subIdex   = cbind(dvect, index2)

	  merged    = merge(orgIdx, subIdex, by.x=1, by.y=1, all=F)
	  merged    = as.matrix(merged)
	  
	  outIndex  = cbind(as.integer(merged[,2]), as.integer(merged[,3]) )
	  mo = order(outIndex[,1])

	  outIndex = outIndex[mo,]

	  return (outIndex)
	}

	permuteVect = function(myvect){
	  pidx = sample(c(1:length(myvect)), length(myvect), replace=F)
	  return(myvect[pidx])
	}

	# a matrix with three cols: MDC random_MDC_mean random_MDC_sd
	#
	MDC_FDR = function(mdcvect){
	  if(mdcvect[1] < 1) {
		 fdr = sum(mdcvect[1]>mdcvect[-1])/length(mdcvect[-1])
	  } else {
		 fdr = sum(mdcvect[1]<mdcvect[-1])/length(mdcvect[-1])
	  }
	  return (fdr)
	}

	# a matrix with three cols: MDC random_MDC_mean random_MDC_sd
	#
	MDC_FDR_by_normal_distr = function(mdcvect){
	  if(mdcvect[1] < 1) {
		 fdr = pnorm(mdcvect[1], mean=mdcvect[2], sd=mdcvect[3], lower.tail=TRUE)
	  } else {
		 fdr = pnorm(mdcvect[1], mean=mdcvect[2], sd=mdcvect[3], lower.tail=FALSE)
	  }
	  return (fdr)
	}

	# xaxisLabel_rotate = 45 # slant x-labels
	#
	barplot_by_oneGroup = function(datMatrix, maskMatrix=NULL, 
					  rowcol=c("red", "green", "blue", "grey"), vlines=NULL,
					  title="", ylab="", imgfile=NULL, imgwid=800, imghei=400, 
					  mcex=1.6, xlab_cex=1, xaxisLabel_rotate=NULL, usr_delta=2, 
					  legendcols=1, legendxDelta=0, legendCex=0.8, 
					  ipointsize=12, iunits="px", ires=72, icompression="lzw",
					  show_value=T, showLegend=T, show_xaxisLabel=T, show_mask=T,
					  xmargin=4, ymargin=4, zmargin=1, xlabel_orientation=0, ValueXdelta=0, 
					  ValueYdelta=NULL, valuecex=0.8, yMAX=NULL, yMIN=0) {

	   xmatrixMj   = datMatrix
	   xmaskMatrix = maskMatrix
	 
	   if(!is.null(maskMatrix)){
		  bar_y = paste(datMatrix[,1], maskMatrix[,1] , sep="")
	   }

	   tab_names  = rownames(xmatrixMj)
	   no.groupsX = dim(xmatrixMj)[1]; no.tabs=dim(xmatrixMj)[2]

		  interv = 1; space4grp = interv+1
		  xcoord =  0.5+ seq(1, space4grp*no.groupsX, space4grp)
		  
	   gene_max = max(xmatrixMj)
	   which_min= as.integer(which.min(gene_max))

	   if( is.null(yMAX) ) {
		  ymax = (as.integer(max(xmatrixMj)/200)+1)*200; #max(ymeans);
		  ymax = (as.integer(max(xmatrixMj)/20)+1)*20; #max(ymeans);
	   } else{
		  ymax = yMAX
	   }

	   if(is.null(yMIN) ) {ymin=min(xmatrixMj)
	   } else {ymin=yMIN}

	   mycolor = rowcol[c(1:no.groupsX)]

	   if(!is.null(imgfile)) {
		 openImgDev(imgfile, iwidth =imgwid, iheight = imghei, ipointsize=ipointsize, iunits=iunits, ires=ires, icompression=icompression)
		 par(mfrow=c(1,1), mar=c(xmargin, ymargin, zmargin, 0) + 0.1, cex=mcex,las=1)#, srt=90)
	   }

	   barplot(datMatrix[,1], beside = TRUE, space=interv,
			col = mycolor, axisnames=F, ylab=ylab, ylim = c(ymin, ymax),
			legend =F)#shortnames) #legend =F

	   if(title!="") {
		   title(main = title, font.main = 1)
	   }

	   #err.bp(t(ymeans), t(ysds), two.side=F) 
	   if(show_xaxisLabel) {
		   if(is.null(xaxisLabel_rotate) ) {
			  axis(1, at =xcoord, labels =rownames(xmatrixMj), las=xlabel_orientation, 
				   col = "black", col.axis="black", lwd = 0, tick = T, line =1)
		   } else{
			  axis(1, at =xcoord, labels =rep("", nrow(datMatrix)), las=xlabel_orientation, 
				   col = "black", col.axis="black", lwd = 0, tick = T, line =1)
			  text(xcoord, par("usr")[3] - usr_delta, srt = xaxisLabel_rotate, adj = 1, labels = rownames(xmatrixMj), 
			  xpd = TRUE, cex=xlab_cex)
			  par(srt=0) 
		   }
	   }

	   ValueYdelta2 = ValueYdelta
	   if(is.null(ValueYdelta)){
		  ValueYdelta2 = ymax/20
	   }

	   if(show_value) {
		  text(x=xcoord, y=datMatrix[,1]+ValueYdelta2, srt=0, labels=bar_y, cex=valuecex, col="black")
	   }

	   if(show_mask) {
		  text(x=xcoord, y=datMatrix[,1]+ValueYdelta2, srt=0, labels=maskMatrix[,1], cex=valuecex, col="black")
	   }


	   if(!is.null(vlines)){
			for(xx in vlines){abline(h=xx,col="gray", lty=3)}
	   }

	   if(!is.null(imgfile)) {
		 par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1, srt=0)
		 dev.off()
	   }

	}

	closest_integer_bound = function(val){
	  logv = log10(val)
	  ii   = as.integer(logv)
	  idec = val/(10^ii) # 2.50
	  iremaind = idec - as.integer(idec)

	  ix = ifelse(iremaind<0.5, 0.5, 1)
	  return ( (as.integer(idec)+ix)*(10^ii) )
	}

	#iunits: "px", "cm", "mm"
	#
	# Letter Size 8.5 x 11 inches
	# images widths: 8.3cm ~ 3.26in, 12.35cm ~ 4.86in, 17.35cm ~ 6.83in
	#        heights: max 23.35cm ~ 9.19in
	#
	openImgDev=function(imgname, iwidth = 1024, iheight=1024, ipointsize=12, iunits="px", ires=72, icompression="lzw")
	{
	  imgtype = getFileExtension(imgname)
	  
	  if (imgtype=="ps"){
		 postscript(imgname,width=iwidth, height=iheight, pointsize=ipointsize)
	  }else if (imgtype=="png"){
		 png(imgname, width=iwidth, height=iheight, pointsize=ipointsize, units=iunits, res=ires)
	  }else if ( (imgtype=="jpg") | (imgtype=="jpeg") ){
		 jpeg(imgname, width=iwidth, height=iheight, pointsize=ipointsize,units=iunits, res=ires,quality =100)
	  }else if ( (imgtype=="tif")|(imgtype=="tiff") ){
		 tiff(imgname, width=iwidth, height=iheight, pointsize=ipointsize,units=iunits, res=ires,compression=icompression)
	  }else if ( (imgtype=="pdf") | (imgtype=="PDF") ){
		 pdf(imgname)
		 return
	  }else{
		 png(imgname, width=iwidth, height=iheight, pointsize=ipointsize, units=iunits, res=ires)
	  }
	  trellis.device(new = FALSE, col = TRUE) 
	}

	collect_garbage=function(){
		while (gc()[2,4] != gc()[2,4]){}
	}

	SetTextContrastColor <- function(color)
	{
	  ifelse( mean(col2rgb(color)) > 127, "black", "white")
	}

	repmat <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

######################################################################################################
#### GO.ENRICH PERFORMS GENE ONTOLOGY ENRICHMENT ANALYSIS BY WGCNA
######################################################################################################

	GO.enrich <- function(moduleColors, module.symbol, GO.enrich.fname, nBestP = 10, organism = "human")
	{  
		#### GO.enrich performs Gene ontology enrichment analysis by WGCNA
		#### Author: Jialiang Yang  
		#### Date: 6/2/2015 version: 1
		##################################################################
		### Input:
		###        moduleColors: color of the gene indicates the module a gene belongs to
		###        module.symbol: the gene symbol (moduleColors and module.symbol must be of the same length)
		###             example:  moduleColors <- c("black", "black", "blue", "blue")
		###                       module.symbol <- c("ZNF593", "ZNF13", "INS", "TP53")
		###     GO.enrich.fname: a file to store the GO enrichment result
		###              nBestP:  top enriched GOs to outupt, default 10 
		###            organism: the species to choose, default "human"
		### Output: top nBestP GO enrichment results for each module (e.g., black and blue)
		###***************************************************************** 

		### Convert gene symbol to entrez ID
		### ----------------------------------------------
			library(org.Hs.eg.db)

			xx <- as.list(org.Hs.egALIAS2EG)
			
			# Remove pathway identifiers that do not map to any entrez gene id
			xx <- xx[!is.na(xx)]

			if(length(xx) > 0){
			# The entrez gene identifiers for the first two elements of XX
			xx[1:2]
			# Get the first one
			xx[[1]]
			}

			module.entrezID.map <- xx[module.symbol]
			names(module.entrezID.map) <- module.symbol
			length(module.entrezID.map)        ### 4

			module.entrezID.list <- sapply(module.entrezID.map, function(current.map)
			{
				if(length(current.map) == 0)
				{
					return("NA")
				}else
				{
					return(current.map[1])
				}
			})

			module.entrezID <- as.character(unlist(module.entrezID.list))
			print(length(module.entrezID))   ### 4
			
		### GO enrichment analysis
		### ----------------------------------------------
			library(WGCNA)		
			library(org.Mm.eg.db)
			library(GO.db)

			GOenr = GOenrichmentAnalysis(moduleColors, module.entrezID, organism = organism, nBestP = nBestP);
			tab = GOenr$bestPTerms[[4]]$enrichment
			names(tab)

			write.table(tab, file =  GO.enrich.fname, sep = ",", quote = TRUE, row.names = FALSE)
	}


######################################################################################################
#### KEGG.ENRICH PERFORMS KEGG ENRICHMENG ANALYSIS BY HYPER-GEOMETRIC TEST
######################################################################################################

	KEGG.enrich <- function(moduleColors, module.symbol, KEGG.enrich.fname, nBestP = 10, organism = "hsa")
	{  
		#### KEGG.enrich performs KEGG enrichment analysis by hyper geomtric test
		#### Author: Jialiang Yang  
		#### Date: 3/24/2016 version: 1
		##################################################################
		### Input:
		###        moduleColors: color of the gene indicates the module a gene belongs to
		###        module.symbol: the gene symbol (moduleColors and module.symbol must be of the same length)
		###             example:  moduleColors <- c("black", "black", "blue", "blue")
		###                       module.symbol <- c("ZNF593", "ZNF13", "INS", "TP53")
		###     KEGG.enrich.fname: a file to store the KEGG enrichment result
		###              nBestP:  top enriched GOs to outupt, default 10 
		###            organism: the species to choose, default "hsa"(e.g. human)
		###                       refer to http://www.genome.jp/kegg/catalog/org_list.html for organism codes
		### Output: top nBestP GO enrichment results for each module (e.g., black and blue)
		###***************************************************************** 

		### Convert gene symbol to entrez ID
		### ----------------------------------------------
			library(org.Hs.eg.db)

			xx <- as.list(org.Hs.egALIAS2EG)
			# Remove pathway identifiers that do not map to any entrez gene id
			xx <- xx[!is.na(xx)]

			if(length(xx) > 0){
			# The entrez gene identifiers for the first two elements of XX
			xx[1:2]
			# Get the first one
			xx[[1]]
			}

			module.entrezID.map <- xx[module.symbol]
			names(module.entrezID.map) <- module.symbol
			length(module.entrezID.map)        ### 4
			
		### retrieve most updated KEGG pathways from KEGG.db
		### Website: http://www.bioconductor.org/packages/release/data/annotation/html/KEGG.db.html
		### ----------------------------------------------
			library(KEGG.db)
			
			### mapps KEGG pathway identifers to Entrez gene or open reading frame identifiers
				xx <- as.list(KEGGPATHID2EXTID)
				if(length(xx) > 0){
				# Get the value of the first key
				xx[[1]]
				# Get the values for multiget for a few keys
				if(length(xx) >= 3){
				xx[1:3]
				}
				}
			
			### retrieve the identifiers only specific to the orgamism we compare, e.g. hsa for human
			    xx.name <- names(xx)
				
				xx.org.char <- sapply(xx.name, function(current.name)
				{
				     return(gsub("[0-9]", "", current.name))
				})
				
				xx.org.ind <- which(xx.org.char == organism)
				
				xx.org <- xx[xx.org.ind]
			
			### mapps KEGG pathway identifers to KEGG pathway names
				yy <- as.list(KEGGPATHID2NAME)
				if(length(yy) > 0){
				# get the value for the first key
				yy[[1]]
				# Get the values for a few keys
				if(length(yy) >= 3){
				yy[1:3]
				}
				}
			
			### combine the KEGG pathway identifers and names
				xx.org.name <- names(xx.org)
				
				### only retrieve the numbers
				xx.org.name.num <- sapply(xx.org.name, function(current.name)
				{
				     return(gsub("[^0-9]", "", current.name))
				})
				
				xx.function.name <- yy[xx.org.name.num]
				xx.org.num <- length(xx.function.name)
				print(xx.org.num)   ### 229
				
				xx.org.ident.func <- sapply(1:xx.org.num, function(current.num)
				{
					return(paste0(xx.org.name[current.num], "_", xx.function.name[current.num]))
				})
				
				names(xx.org) <- xx.org.ident.func
				
				print(xx.org[1:2])
	
		### generate background universal genes: we use a combination of genes in GO and KEGG
		### ----------------------------------------------
			gene.KEGG <- unique(unlist(xx.org))
			print(length(gene.KEGG))  ### 5894
		
			gene.universal <- gene.KEGG
			
	    ### Enrichment analysis by hyper-geometric test
		### ----------------------------------------------	
		
			### retrive the genes in each module
			module.names <- unique(moduleColors)
			
			module.enrich <- lapply(module.names, function(current.module)
			{
				current.module.gene <- as.character(module.entrezID.map)[which(moduleColors == current.module)]	
				drawn.gene <- intersect(current.module.gene,gene.universal)
                numDrawn<- length(drawn.gene)				
				
				### hyper geometric test for each KEGG pathways and rank according to significance
				enrich.vec<- sapply(xx.org, function(current.KEGG)
				{
					white.gene<- current.KEGG
					black.gene<-setdiff(gene.universal,white.gene)
					white.drawn.gene<- intersect(white.gene,drawn.gene)
					numW<-length(white.gene)
					numB<-length(black.gene)
					numWdrawn<-length(white.drawn.gene)
					phyper(numWdrawn - 1L, numW, numB, numDrawn, lower.tail=FALSE)
				})
				
				enrich.FDR <- p.adjust(enrich.vec, method = "bonferroni")
				
				### return top nBestP KEGG pathway according to FDR
				module <- rep(current.module, nBestP)
				modSize <- rep(length(current.module.gene), nBestP)
				bkgrModSize <- rep(length(drawn.gene), nBestP)
				rank.vec <- 1:nBestP
				
				top.KEGG <- sort(enrich.FDR)[1:nBestP]
				top.KEGG.name <- names(top.KEGG)
				
				enrich.matrix <- cbind(module, modSize, bkgrModSize, rank.vec, enrich.vec[top.KEGG.name], top.KEGG, top.KEGG.name)
			})
			
			module.enrich.matrix <- do.call(rbind, module.enrich)
			
			colnames(module.enrich.matrix) <- c("module", "modSize", "bkgrModSize", "rank", "enrihmentP", "BonferoniP", "KEGGTerm")	
			
			write.table(module.enrich.matrix, file = KEGG.enrich.fname, sep = ",", quote = TRUE, row.names = FALSE)
	}
	
######################################################################################################
#### WRITE.GENESET WRITES A LIST OF GENE SET INTO A FILE
######################################################################################################	

	### returns string w/o leading or trailing whitespace
	### funciton trim removes spaces at the begining and end of a sring or list
	  trim <- function (x) gsub("^\\s+|\\s+$", "", x)

	write.geneset <- function(gene.set.list, out.fname)
	{
	  #### write.geneset writes a list of gene set into a file
	  #### Author: Jialiang Yang  
	  #### Date: 12/8/2015 version: 1
	  ##################################################################
	  ### Input:
	  ###        gene.set.list: gene set list
	  ###        out.fname: output file name
	  ###***************************************************************** 
	  
	  set.num <- length(gene.set.list)
	  set.name <- names(gene.set.list)
	  
	  sink(out.fname)
	  
	  for (i in 1: set.num)
	  {
		current.name <- set.name[i]
		current.geneset <- trim(gene.set.list[[i]])
		unique.geneset <- unique(current.geneset)
		gene.num <- length(unique.geneset)
			
		if (gene.num == 0) 
		{
			cat(current.name,'\t genes_0 \n')
		}else
		{
			if(gene.num == 1)
			{
				cat(current.name,'\t genes_1', unique.geneset,  '\n')
			}else
			{
				cat(current.name,'\t',paste('genes',gene.num,sep='_'),'\t',sep = '')

				for(j in 1:(gene.num-1))
				{
					cat(unique.geneset[j],'\t',sep='')
				}
						
				cat(unique.geneset[gene.num],'\n', sep= '') 
			}
			
		}
	   }

	   sink()
	 }
	 
######################################################################################################
#### READ.GENESET READS THE GENE SET FILE GENERATED BY CAT
######################################################################################################

	read.geneset<- function(gene.set.file)
	{
	  #### read.geneset reads the gene set file generated by cat
	  #### Author: Jialiang Yang  
	  #### Date: 3/19/2015
	  ##################################################################
	  ### Input: 
	  ###        gene.set.file: the gene set file name with format setname \t anything \t gene1 \t gene2 \t ... in each role
	  ### Output: 
	  ###        set.gene: a list file stores the unique genes in each set, i.e. names(set.gene) = set name set.gene[[1]] <- c(gene1, gene2, ...)
	  ###*****************************************************************
	  
	  pathway.string <- readLines(gene.set.file, -1)
	  
	  pathway.list <- strsplit(pathway.string,'\t')
	  
	  ### retrive pathway names
	  set.name <- unlist(lapply(pathway.list, function(y)
	  {
		 return(y[1])
	  }))
	  
	  ### retrieve pathway genes
	  set.gene <- lapply(pathway.list, function(y)
	  {
		 y.len <- length(y)
		 y.vec <- y[3:y.len]
		 y.vec <- unique(setdiff(y.vec, ""))
	  })
	  
	  names(set.gene)<- set.name
	  return(set.gene)
	}

######################################################################################################
#### read.CMAP READS THE GENE SET FILE GENERATED BY CMAP
######################################################################################################

	read.CMAP<- function(gene.set.file)
	{
	  #### read.geneset reads the gene set file generated by cMap
	  #### Author: Jialiang Yang  
	  #### Date: 8/15/2017
	  ##################################################################
	  ### Input: 
	  ###        gene.set.file: the gene set file name with format setname \t gene1 \t gene2 \t ... in each role
	  ### Output: 
	  ###        set.gene: a list file stores the unique genes in each set, i.e. names(set.gene) = set name set.gene[[1]] <- c(gene1, gene2, ...)
	  ###*****************************************************************
	  
	  pathway.string <- readLines(gene.set.file, -1)
	  
	  pathway.list <- strsplit(pathway.string,'\t')
	  
	  ### retrive pathway names
		  set.name <- unlist(lapply(pathway.list, function(y)
		  {
			 return(y[1])
		  }))
	  
	  ### retrieve pathway genes
		  set.gene <- lapply(pathway.list, function(y)
		  {
			 y.len <- length(y)
			 y.vec <- y[2:y.len]
			 y.vec <- unique(setdiff(y.vec, ""))
		  })
		  
	  names(set.gene)<- set.name
	  return(set.gene)
	}
	
######################################################################################################
#### read.L1000 READS THE GENE SET FILE GENERATED BY L1000
######################################################################################################

	read.L1000 <- function(gene.set.file)
	{
	  #### read.geneset reads the gene set file generated by L1000
	  #### Author: Jialiang Yang  
	  #### Date: 8/15/2017
	  ##################################################################
	  ### Input: 
	  ###        gene.set.file: the gene set file name with format setname \t gene1,1.0 \t gene2,1.0 \t ... in each role
	  ### Output: 
	  ###        set.gene: a list file stores the unique genes in each set, i.e. names(set.gene) = set name set.gene[[1]] <- c(gene1, gene2, ...)
	  ###*****************************************************************
	  
	  pathway.string <- readLines(gene.set.file, -1)
	  
	  pathway.list <- strsplit(pathway.string,'\t')
	  
	  ### retrive pathway names
		  set.name <- unlist(lapply(pathway.list, function(y)
		  {
			 return(y[1])
		  }))
	  
	  ### retrieve pathway genes
		  set.gene <- lapply(pathway.list, function(y)
		  {
			 y.len <- length(y)
			 y.vec <- gsub(",1.0", "", y[2:y.len])
			 y.vec <- unique(setdiff(y.vec, ""))
		  })
		  
	  names(set.gene)<- set.name
	  return(set.gene)
	}
	
######################################################################################################
#### FILTER.GENESET FILTERS THE GENESET
######################################################################################################

	filter.geneset<- function(geneset,genesymbol,size= 5)
	{
	  ##### filter.geneset filters the geneset with size less than size
	  ####################################################################
	  ### Input:
	  
	  ### Output:
	  ####################################################################
	  
	  new.geneset<- lapply(1:length(geneset),function(x){ intersect(geneset[[x]],genesymbol) })
	  names(new.geneset)<- names(geneset)
	  new.geneset.size<- unlist(lapply(1:length(new.geneset),function(x){ length(new.geneset[[x]]) }))
	  new.geneset[new.geneset.size>=size]
	}

	read.filter.geneset<- function(gene.set.file,genesymbol,size = 5)
	{
	   geneset<- read.geneset(gene.set.file)
	   filter.geneset(geneset,genesymbol,size)
	}

######################################################################################################
### spanning.network 
######################################################################################################	

    spanning.network <- function(gene.set, network.matrix, method = "SP")
	{
	  ##### map.network maps the aging ang disease genes into the modularized network
	  ####################################################################
	  ### Input:
	  ###        gene.set: the set of genes to be kept in the spanning network (in gene symbol)
      ###        network.matrix: the whole network(adjacency matrix) with row and column names the genes(symbol)
	  ###        method: 
	  ###               "SP": add all genes in the shortest path between genes in the gene set into the network
	  ###               "Steiner": the steiner tree spanning network
	  ###               "overlap": the subnetwork containing only overlapped geens (with edges) between gene.set and network.matrix 
      ### Output:
      ###        network.span: the spanning network	  
	  ####################################################################
	  
		  gene.symbol <- rownames(network.matrix)
		  gene.mapped <- intersect(gene.set, gene.symbol)
		  print(length(gene.mapped))
		  
		  ### method 1: shortest path
		  if (method == "SP")
		  {
			### generate graph
				library(igraph)
				network.graph <- graph_from_adjacency_matrix(network.matrix, mode ="undirected", weighted = NULL, diag = FALSE)
				
				network.vertices.ID <- V(network.graph)$name
				
				from.vec <- match(gene.mapped, network.vertices.ID)
				to.vec <- match(gene.mapped, network.vertices.ID)
				
				nodes.expand.list <- sapply(from.vec, function(current.node)
				{
				    current.path <- shortest_paths(network.graph, from = current.node, to = to.vec, output = "vpath")
					current.path.vec <- unique(names(unlist(current.path$vpath)))
					return(current.path.vec)
				})
				
				node.expand.vec <- unique(unlist(nodes.expand.list))
		  }
		  
		  network.span <- network.matrix[gene.symbol %in% node.expand.vec, gene.symbol %in% node.expand.vec]
		  network.span.final <- network.span[rowSums(network.span)!=0,colSums(network.span)!=0]
		  
		return( network.span.final)
	}
	
######################################################################################################
### MAP.NETWORK MAPS THE AGING ANG DISEASE GENES INTO THE MODULARIZED NETWORK WITHOUT EXPANDING
######################################################################################################

	map.network <- function(network.data.fname, GENAGE.fname, GENE.THRESH, DISEASE.fname, geneset.save.fname, module.fname, 
	                        PERMUTE.N, network.span.method = "SP", protien.coding.symbol)
	{
	  ##### map.network maps the aging ang disease genes into the modularized network without expanding
	  ####################################################################
	  ### Input:
	  ###        network.data.fname: an RData file to store the whole netwrok 
	  ###                            contain "intM" format: ITGA7 PPP1R9A SRGN GRB7 PAK1
	  ###                                           ITGA7       0       0    0    0    0
	  ###                                           PPP1R9A     0       0    0    0    0
	  ###                                           SRGN        0       0    0    0    0
	  ###                                           GRB7        0       0    0    0    0
	  ###                                           PAK1        0       0    0    0    0
      ###              GENAGE.fname: aging gene file in csv format 
	  ###                                           GenAge ID	symbol	aliases	   name	                    why	entrez gene id	...
      ###                                                  1	GHR	      GHBP	 growth hormone receptor	mammal	2690	   ...
	  ###               GENE.THRESH: miminum numbers of genes mapped
	  ###             DISEASE.fname: file\t deliminated) to save disease genes format: 
	  ###                                                       Bipolar disorder	genes_35	ADCY2	ANK3	ODZ4	TRANK1	MIR2113	
      ###                                                         Thyroid cancer	genes_44	DIRC3	FOXE1	NRG1	MBIP	NKX2-1
	  ###        geneset.save.fname: a folder to save results
	  ###              module.fname:   the module (\t deliminated) file the network restricted to format: black	genes_944	PLEKHN1	C1orf159	TTLL10	ACAP3 ...
	  ###                 PERMUTE.N: number of permutation times
	  ###        network.span.method: 
	  ###               "SP": add all genes in the shortest path between genes in the gene set into the network
	  ###               "Steiner": the steiner tree spanning network
	  ###               "overlap": the subnetwork containing only overlapped geens (with edges) between gene.set and network.matrix 
	  ###          protien.coding.symbol: a vector to store protien coding gene symbols
	  ### Output:
	  #################################################################### 
	  
	  ### deal with big network and aging and diesease genes
	  ### ----------------------------------------------------------------
	      prefix.full.name.vec <- unlist(strsplit(network.data.fname, split = "/"))
		  prefix.full.name <- prefix.full.name.vec[length(prefix.full.name.vec)]
		  prefix.name.vec <- unlist(strsplit(prefix.full.name, split = "[.]"))
		  prefix.name <- paste0(prefix.name.vec[1:(length(prefix.name.vec)-1)], collapse = ".")
		  print(prefix.name)
		  
		  load(network.data.fname)   ### Network.RData

		  print(dim(intM))   #  10687 10687
		  print(sum(intM)/2) # 206674

		  genesymbol<- rownames(intM)
		  print(genesymbol[1:5])

		  GenAge.matrix <- read.csv(GENAGE.fname, header=TRUE, stringsAsFactors=FALSE)
		  print(GenAge.matrix[1:5,])
		  
		  GenAge<- intersect(GenAge.matrix[,"symbol"],genesymbol)
		  print(length(GenAge))  ### 292

		  disease.geneset<- read.filter.geneset(DISEASE.fname,genesymbol, GENE.THRESH)
		  print(disease.geneset[1:5])
		  print(length(disease.geneset)) ###  224
		  
		  actual.geneset<- c(list(GenAge=GenAge),disease.geneset)
		  length(actual.geneset)  ###  225

		  geneset.save.vector <- unlist(strsplit(geneset.save.fname, split = "/"))
		  fold.name <- paste(geneset.save.vector[1:(length(geneset.save.vector)-1)], collapse = "/")

		  if(!file.exists(fold.name))
		  {
			 dir.create(fold.name)
		  }
		  
		  save(intM,genesymbol,actual.geneset,file = geneset.save.fname)   ### STRING.cut0.900_geneset.RData
		  
	   ### creat folders to save intermidiate results
	   ### ----------------------------------------------------------------
		  work.fold.vec <- unlist(strsplit(geneset.save.fname, "/"))
		  work.fold <- paste(work.fold.vec[1:(length(work.fold.vec)-2)], collapse = "/")
		  
		  input.fold.name <- paste0(work.fold, "/input/")
		  
		  if(!file.exists(input.fold.name))
		  {
			 dir.create(input.fold.name)
		  }

	  ### deal with modules: map the selected modules into the network and do permutation
	  ### ----------------------------------------------------------------
		subnetwork.geneset <- read.geneset(module.fname)
		library(igraph)
		  
		module.size.matrix <- t(sapply(1:length(subnetwork.geneset), function(i)
		{
		    gene.set.i <- subnetwork.geneset[[i]]
			print(length( gene.set.i))
			
			gene.set.protein.coding <- intersect(gene.set.i, protien.coding.symbol)
			
			overlap.gene.symbol <- intersect(gene.set.i, genesymbol)
			print(length(overlap.gene.symbol))
			
			if (length(overlap.gene.symbol) <= 1)
			{
				return(c(length(gene.set.i), length(gene.set.protein.coding), length(overlap.gene.symbol), 0, 0)) 
			}else
			{
			    network.small <- intM[rownames(intM) %in% gene.set.i,colnames(intM) %in% gene.set.i]
				network.small.shrink <- network.small[rowSums(network.small)!=0,colSums(network.small)!=0]
				size.small <- dim(network.small.shrink)[1]
				
				### construct spannning network from the mapped gene
				   intM.i <- spanning.network(gene.set.i, intM, method = network.span.method) 
				   print(dim(intM.i))
				
			    ### retrieve modules if larger than GENE.THRESH genes (in the module) could be mapped
				   gene.i<- rownames(intM.i)
				   actual.geneset.i <- filter.geneset(geneset=actual.geneset,genesymbol=gene.i,size= GENE.THRESH)

				### if at least geneAge and another disease is enriched in thepathway, do the permutation for gene age 
				   if(('GenAge' %in% names(actual.geneset.i))&(length(actual.geneset.i)>=2))
				   {
					  geneset.size<- length(actual.geneset.i[['GenAge']])
					  
					  permute.geneset.i<- lapply(1:PERMUTE.N,function(i){
					  sort(sample(gene.i,geneset.size,replace=F))})
					  
					  names(permute.geneset.i)<- paste0('Size',geneset.size,'_',1:PERMUTE.N)
				   
					  geneset.i<- c(actual.geneset.i,permute.geneset.i)
					  # length(geneset.i)

					  save(gene.i,intM.i,actual.geneset.i,permute.geneset.i,geneset.i,file=paste0(input.fold.name,prefix.name,
						'_',gsub(':','',fixed=T,strsplit(names(subnetwork.geneset)[i],'_')[[1]][1]),'.RData'))
				   }
				   
				return(c(length(gene.set.i), length(gene.set.protein.coding), length(overlap.gene.symbol),size.small, dim(intM.i)[1])) 
			}
		}))
		
		rownames(module.size.matrix) <- names(subnetwork.geneset)
		colnames(module.size.matrix) <- c("Gene Num", "Protein coding","Overlap with PPI", "GeroNet Min", "GeroNet SP")
	   
	   return(module.size.matrix)
	}
	
######################################################################################################
### RWR PERFORMS RANDOM NETWORK WITH RESTART
######################################################################################################

	rwr <- function(W,P0,gamma,converge.delta) 
	{
	   #####  rwr performs random network with restart
	   ####################################################################
	   ### Input:
	  
	   ### Output:
	   #################################################################### 
	   
		PT <- P0
		k <- 0
		delta <- 1
		while  (delta > converge.delta) {
		  PT1 <- (1-gamma)*W
		  PT2 <- PT1 %*% t(PT)
		  PT3 <- (gamma*P0)
		  PT4 <- t(PT2) + PT3
		  delta <- sum(abs(PT4 - PT))
		  PT <- PT4
		  k <- k + 1
		}
		return(PT)
	}

	randomWalk.matrix <- function(intM, queryGenes, p0=p0, gamma=NA, converge.delta=1e-6) 
	{
	  #####  randomWalk.matrix performs random network with restart until converge
	  ####################################################################
	  ### Input:
	  
	  ### Output:
	  #################################################################### 
	   
		p0[queryGenes] <- 1
		p0 <- p0/sum(p0)
		res <- rwr(intM,t(p0),gamma,converge.delta)	
		return(drop(res))
	}

######################################################################################################
### ASSOCIATION.PATHWAY CALCULATE AGE-DISEASE ASSOCIATION ON EACH PATHWAY
######################################################################################################

	association.pathway <- function(subnetwork.data.fanme, gammavalue, output.folder, prefix) 
	{
	   #####  association.pathway calculate age-disease association on each pathway
	   ####################################################################
	   ### Input:
	  
	   ### Output:
	   #################################################################### 
	   
	   load(subnetwork.data.fanme)
	   genesymbol <- gene.i
	   geneset <- geneset.i
	   intM <- intM.i
	   
	   score.matrix<- matrix(NA,nrow=length(genesymbol),ncol=length(geneset),dimnames=list(genesymbol,names(geneset)))
	   
	   Ng <- dim(intM)[1]   
	  
	   for (i in 1:Ng) 
	   {
		 intM[,i] <- intM[,i]/sum(intM[,i])
	   }

	   p0 <- numeric(length=Ng)
	   names(p0) <- row.names(intM)
	   
	   print(dim(score.matrix))
	   
	   ### perform random walk for aging and all disease genes
		   for(i in 1:ncol(score.matrix))
		   {
			 print(i)
			 queryGenes<- geneset[[i]]
			 try(yy<- randomWalk.matrix(intM, queryGenes, p0=p0, gamma=gammavalue))
			 try(score.matrix[,i]<- yy[genesymbol])
		   }
	   
	   ### calculate enrichment score
		   actual.geneset<- actual.geneset.i
		   GenAge<- actual.geneset[['GenAge']]
	   
	   ### permutation
		   permute.GenAge<- names(permute.geneset.i)
		   actual.permute.GenAge<- c('GenAge',permute.GenAge)
		   acutal.disease<- setdiff(names(actual.geneset),'GenAge')
		   print(length(acutal.disease))
	   
	   
	   ### disease to geneAge: geneAge as seeded gene set
		   disease2GenAge<- matrix(NA,nrow=length(actual.permute.GenAge),ncol=length(acutal.disease),dimnames=list(actual.permute.GenAge,acutal.disease))

		   N<- nrow(score.matrix)
		   G<- length(GenAge)

		   minus<- sqrt( G/(N - G) )
		   add<- sqrt( (N - G)/G )
		   
		   for(j in 1:ncol(disease2GenAge))
		   {
			  # print(j)
			  score<- sort(score.matrix[,acutal.disease[j]],decreasing=T)
			  all.ranked.gene<- names(score)

			  index<- c(which(!duplicated(score))[-1] -1,length(score))
			  ES.NA<- rep(NA,length=length(index))

			  disease2GenAge[,j]<- unlist(lapply(1:length(actual.permute.GenAge),function(i)
			  {
				 # cat(i,' : ',j,'\n',sep='')
				 ES.change<- rep(-minus,length(all.ranked.gene))
				 ES.change[all.ranked.gene %in% geneset[[actual.permute.GenAge[i]]]]<- +add

				 ES<- ES.NA
				 ES[1]<- sum(ES.change[1:index[1]])
				 for(ix in 2:length(index)){ ES[ix]<- ES[ix-1] + sum(ES.change[(index[ix-1]+1):index[ix]]) }

				 max(ES)
			  }))
		   }
	   
	   ### GeneAge to disease: disease as seeded gene set
		   GenAge2disease <- matrix(NA,nrow=length(acutal.disease),ncol=length(actual.permute.GenAge),dimnames=list(acutal.disease,actual.permute.GenAge))
		   
		   N<- nrow(score.matrix)
		   
		   for(j in 1:ncol(GenAge2disease))
		   {
			  print(j)
			  score<- sort(score.matrix[,actual.permute.GenAge[j]],decreasing=T)
			  all.ranked.gene<- names(score)
			  
			  index<- c(which(!duplicated(score))[-1] -1,length(score))
			  ES.NA<- rep(NA,length=length(index))
			  
			  GenAge2disease[,j]<- unlist(lapply(1:length(acutal.disease),function(i)
			  {
				 # cat(i,' : ',j,'\n',sep='')
				 disease.gene<- geneset[[acutal.disease[i]]]

				 G <- length(disease.gene)
				 minus<- sqrt( G/(N - G) )
				 add <- sqrt( (N - G)/G )
				 
				 ES.change<- rep(-minus,length(all.ranked.gene))
				 ES.change[all.ranked.gene %in% disease.gene]<- +add
				 
				 ES<- ES.NA
				 ES[1]<- sum(ES.change[1:index[1]])
				 for(ix in 2:length(index)){ ES[ix]<- ES[ix-1] + sum(ES.change[(index[ix-1]+1):index[ix]]) }
			  
				 max(ES)
			  }))
		   }
	   
	   ### calculate the weighted sum of the scores, e.g. 
		   beta.list<- seq(from=0,to=1,by=0.1)
		   
		   pvalue.list<- vector('list',length(beta.list))
		   names(pvalue.list)<- beta.list

		   for(j in 1:length(beta.list))
		   {
			  beta<- beta.list[j]
			  
			  pvalue<- unlist(lapply(acutal.disease,function(disease){
			  
			  ### the combined score
			  ES.actual<- beta*disease2GenAge['GenAge',disease] + (1 - beta)*GenAge2disease[disease,'GenAge']
			  ES.permute<- beta*disease2GenAge[permute.GenAge,disease] + (1 - beta)*GenAge2disease[disease,permute.GenAge]

			  z<- (ES.actual - mean(ES.permute))/sd(ES.permute)
			  pnorm(q=z,mean=0,sd=1,lower.tail=F,log.p=F)
			  }))

			  fdr<- p.adjust(pvalue,method='fdr')
			  
			  pvalue.list[[j]]<- cbind(disease=acutal.disease,pvalue,fdr)
		   }
		   
		   if(!file.exists(output.folder))
		   {
			  dir.create(output.folder)
		   }
		   
		   save(score.matrix,disease2GenAge,GenAge2disease,pvalue.list,file=paste0(output.folder,prefix,'_gammavalue',gammavalue,'.RData'))
	}
	
######################################################################################################
### SUMMARY.MNA SUMMARIZES THE ASSOCIATION P-VLAUES FOR EACH SUBNETWORK
######################################################################################################

	summary.MNA <- function(output.fold, prefix.list, module.fname, summary.fold, beta.list) 
	{
	   ##### summary.MNA summarizes the association p-vlaues for each subnetwork
	   ####################################################################
	   ### Input:
	  
	   ### Output:
	   #################################################################### 
	   
	   file.list<- sort(dir(output.fold,pattern='[.]RData'))    ### result from 2 template.R
	   print(prefix.list)
	   
	   if(!file.exists(summary.fold))
	   {
		  dir.create(summary.fold)
	   }
	   
	   ### calculate one by one
	   for(i in 1:length(prefix.list))
	   {
		   pvalue.list.xx<- vector('list',length(beta.list))
		   names(pvalue.list.xx)<- beta.list
		   
		   subnetwork.geneset <- read.geneset(module.fname)
		   
		   select.geneset.name <- names(subnetwork.geneset)
		   print(select.geneset.name)
		   
		   names(select.geneset.name)<- paste0(prefix.list[i],'_',gsub(':','',fixed=T,
			   unlist(lapply(select.geneset.name,function(x){ strsplit(x,'_')[[1]][1] }))),'_gammavalue0.7.RData')
			   
		   file.list.i<- file.list[grep(prefix.list[i],file.list)]
		   print(length(file.list.i))
		   
		   for(j in 1:length(file.list.i))
		   {
			  print(j)
			  load(paste0(output.fold,file.list.i[j]))
			  
			  for(k in 1:length(beta.list))
			  {
				  beta<- beta.list[k]
				  
				  if(beta %in% seq(from=0,to=1,by=0.1))
				  {
					 tab.x<- pvalue.list[[as.character(beta)]]
				  } else {
					  acutal.disease<- colnames(disease2GenAge)
					  actual.permute.GenAge<- rownames(disease2GenAge)
					  permute.GenAge<- actual.permute.GenAge[actual.permute.GenAge!='GenAge']
					  
					  pvalue<- unlist(lapply(acutal.disease,function(disease){
						  ES.actual<- beta*disease2GenAge['GenAge',disease] + (1 - beta)*GenAge2disease[disease,'GenAge']
						  ES.permute<- beta*disease2GenAge[permute.GenAge,disease] + (1 - beta)*GenAge2disease[disease,permute.GenAge]
						  
						  z<- (ES.actual - mean(ES.permute))/sd(ES.permute)
						  pnorm(q=z,mean=0,sd=1,lower.tail=F,log.p=F)
					  }))
					  
					  fdr<- p.adjust(pvalue,method='fdr')
					  tab.x<- cbind(disease=acutal.disease,pvalue,fdr)
					}
					
				 pvalue.list.xx[[as.character(beta)]]<- rbind(pvalue.list.xx[[as.character(beta)]],
								  cbind(subnetwork=select.geneset.name[file.list.i[j]],tab.x))
			  }
		   }
		   
		   for(k in 1:length(beta.list))
		   { 
			  beta<- as.character(beta.list[k])
			  tab<- pvalue.list.xx[[beta]]
			  pvalue<- as.numeric(tab[,'pvalue'])
			  tab[,'fdr']<- p.adjust(pvalue,method='fdr')
			  write.table(tab[order(pvalue),],file=paste0(summary.fold,prefix.list[i],'_beta',beta,'_pvalue.csv'),row.names=F,sep=',',quote=T)
			}
		}
	}

######################################################################################################
#### DESEQ2 PIPELINE
######################################################################################################
	
	run.DESeq2 <- function(sampleTable, DATA.FOLD, design, DIFF.FOLD, cond.vec, alpha)
	{  
	  #### run.DESeq.simple: performs the differential gene analysis from raw reads count
	  ####            by DEseq using simple design (two groups)
	  ##################################################################
	  ### Input:
	  ###       meta.table.fname: the metadata table file
	  ###       DATA.FOLD: the folder to store count files
	  ###       DIFF.FOLD: the folder stores all the results
	  ###       design: the design formula, the differential parameter should be in the end
	  ###       cond.vec: the condition to compare c("Young", "Old")
	  ###       alpha: p-value cutoff
	  ### Output: 
	  ###       error.vec: the summary of read depth and no consistent reads for each quality score
	  ###***************************************************************** 
	  
	  library("DESeq2")
	  library("BiocParallel")
	  library("ggplot2")
	  ###library("pheatmap")
	  ###library("vsn")
	  
	  register(MulticoreParam(4))
	  
	  setwd(DATA.FOLD)
	  
	  dds <- DESeqDataSetFromHTSeqCount(sampleTable, DATA.FOLD, design)
	  dds <- dds[ rowSums(counts(dds)) > 1, ]
		  
	  dds$condition <- factor(dds$condition, levels=cond.vec)
	  
	  ### Call differential genes
		  dds <- DESeq(dds, parallel=TRUE)
		  res <- results(dds, parallel=TRUE, alpha = alpha)
		  
		  resOrdered <- res[order(res$padj),]
		  
		  write.table(resOrdered, file = paste0(DIFF.FOLD,"res_DESeq", "_FDR_", alpha,".txt"),quote=F,col.names=NA,sep='\t')
		  
		  print(summary(res))
		  
		  print(sum(res$padj < alpha, na.rm=TRUE))
	  
	  ### plot the log2 fold changes attributable to a given variable over the mean of normalized counts
		  pdf(file = paste0(DIFF.FOLD,"log2FoldChange.pdf"), width = 7, height = 7)
			plotMA(res, main="DESeq2", ylim=c(-2,2))
			idx <- identify(res$baseMean, res$log2FoldChange)
			resMLE <- results(dds, addMLE=TRUE, parallel=TRUE)
			plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))
		  dev.off()
		  
	  ### plot Counts
	  	  pdf(file = paste0(DIFF.FOLD,"plotCounts.pdf"), width = 7, height = 7)
			plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
			
			d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)
			ggplot(d, aes(x=condition, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(25,100,400))
		  dev.off()
		
		sink(file = paste0(DIFF.FOLD,"MethodDescription.txt"))
		     mcols(res)$description
		sink()
		
	   ## Extracting tranformed values
		# rld <- rlog(dds, blind=FALSE)
		# vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
		# vsd.fast <- vst(dds, blind=FALSE)
		# head(assay(rld), 3)

		
	  ### Effects of transformations on the variance
	    # pdf(file = paste0(DIFF.FOLD,"variance.pdf"), width = 7, height = 7)
			# notAllZero <- (rowSums(counts(dds))>0)
			# meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
			# meanSdPlot(assay(rld[notAllZero,]))
			# meanSdPlot(assay(vsd[notAllZero,]))
		# dev.off()
		
	  ### Heatmap of the count matrix
	     # pdf(file = paste0(DIFF.FOLD,"heatmapCount20.pdf"), width = 7, height = 7)
			# select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
			# nt <- normTransform(dds) # defaults to log2(x+1)
			# log2.norm.counts <- assay(nt)[select,]
			# df <- as.data.frame(colData(dds)[,c("condition", "type")])
			
			# pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)
			
			# pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
			
			# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)
			
		# dev.off()
	
	}
######################################################################################################
#### DESEQ SIMPLE PIPELINE
######################################################################################################

    run.DESeq.simple <- function(Sample.fname,  DIFF.FOLD, condA, condB)
	{  
	  #### run.DESeq.simple: performs the differential gene analysis from raw reads count
	  ####            by DEseq using simple design (two groups)
	  ##################################################################
	  ### Input:
	  ###       Sample.fname: the metadata table file
	  ###       DIFF.FOLD: the folder stores all the related files
	  ### Output: 
	  ###       error.vec: the summary of read depth and no consistent reads for each quality score
	  ###***************************************************************** 
	  
	  setwd(DIFF.FOLD)
	  
	  ### prepare data
		  Sample.matrix <- read.delim(Sample.fname, header = T, quote = "", row.names = 1,  stringsAsFactors=FALSE)
		  Sample <- as.data.frame(Sample.matrix)
		  
		  samplesDESeq = with(Sample, data.frame(shortname = I(shortname), countf = I(countf),condition = condition))
		  
		  library("DESeq")
		  cds = newCountDataSetFromHTSeqCount(samplesDESeq)
		  
		  cds = estimateSizeFactors(cds)
		  sizeFactors(cds)
			### MOAloneN2 MOAloneN3    MOCON1    MOCON2    MOCON3
			### 0.6213529 0.6443940 0.7582592 1.4401263 2.3792545
				
	  ### To inspect sample relationships, invoke a variance-stabilizing transformation and inspect a principal component analysis (PCA) plot
		  cdsB = estimateDispersions(cds, method = "blind")
		  vsd = varianceStabilizingTransformation(cdsB)
			
		  pdf(file = "Sample.PCA.pdf", width = 7, height = 7)
			p = plotPCA(vsd, intgroup = c("condition"))
		  dev.off()
		
	  ### Use estimateDispersions to calculate dispersion values
		  cds = estimateDispersions(cds)
		  
		  pdf(file = "plot.dispEsts.pdf", width = 7, height = 7)
		  plotDispEsts(cds)
		  dev.off()
		
	  ### Perform the test for differential expression by using nbinomTest
		  res = nbinomTest(cds, condA, condB)
			
		  pdf(file = "plot.res.pdf", width = 7, height = 7)
			plotMA(res)
		  dev.off()
			
	  ### Inspect the result tables of significantly upregulated and downregulated genes, at a 10% false discovery rate (FDR)
		  resSig = res[which(res$padj < 0.1),]
		  print(head( resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ] ))
		  print(head( resSig[ order(resSig$log2FoldChange, decreasing = FALSE), ] ))
		
	  ### Count the number of genes with significant differential expression at a FDR of 10%
		  print(table( res$padj < 0.1 ))			
				### FALSE  TRUE
				### 46099   196
				
	  ### Write into a file
		  write.table(res, file = "res_DESeq.txt",quote=F,col.names=NA,sep='\t')
		  
		  pdf(file = "hist.pdf", width = 7, height = 7)
			hist(res$pval, breaks = 100)
		  dev.off()
	}

	
######################################################################################################
#### edgeR SIMPLE PIPELINE
######################################################################################################

    run.EDGER.simple <- function(Sample.fname, DIFF.FOLD, condA, condB)
	{  
	  #### run.edgeR.simple: performs the differential gene analysis from raw reads count
	  ####            by edgeR using simple design (two groups)
	  ##################################################################
	  ### Input:
	  ###       Sample.fname: the metadata table file
	  ###       DIFF.FOLD: the folder stores all the related files
	  ###***************************************************************** 
	  
	  setwd(DIFF.FOLD)
	  
	  ### prepare data
		  Sample.matrix <- read.delim(Sample.fname, header = T, quote = "", row.names = 1,  stringsAsFactors=FALSE)
		  Sample <- as.data.frame(Sample.matrix)
		    
		  library("edgeR")
		  counts = readDGE(Sample$countf)$counts
		  
	  ### Filter weakly expressed and noninformative (e.g., non-aligned) features
		  noint = rownames(counts) %in% c("no_feature","ambiguous","too_low_aQual", "not_aligned","alignment_not_unique")
		  cpms = cpm(counts)
		  keep = rowSums(cpms > 1) >= 2 & !noint
		  ###In edgeR, it is recommended to remove features without at least 1 read per million in n of the
          #### samples, where n is the size of the smallest group of replicates (here, n = 2 for the Alone group)
		  counts = counts[keep,]

	  ### Visualize and inspect the count table as follows
          colnames(counts) = Sample$shortname	  
		  head( counts[,order(Sample$condition)], 5 )
		  
							 # ADAloneN2 ADAloneN3 ADCON1 ADCON2 ADCON3
			# ENSG00000000419       604       326    581   1152    585
			# ENSG00000000457       138        98    102    258    133
			# ENSG00000000460        43        28     20     92     46
			# ENSG00000000971      4462      2492   5789  13328   6698
			# ENSG00000001036       753       476    865   1353    791

	   ### Create a DGEList object (edgeR s container for RNA-seq count data)  
		   d = DGEList(counts = counts, group = Sample$condition)
		   d = calcNormFactors(d)   ### Estimate normalization factors
				
	   ### Figure to inspect the relationships between samples using a multidimensional Scaling (MDS) plot
	        pdf(file = "plotMDS.pdf", width = 7, height = 7)
				plotMDS(d, labels = Sample$shortname, col = c("darkgreen","blue")[factor(Sample$condition)])
			dev.off()
			
		### Estimate tagwise dispersion (simple design)
		    d = estimateCommonDisp(d)
			d = estimateTagwiseDisp(d)
		
        ### Create a visual representation of the mean-variance relationship using the plotMeanVar		
			pdf(file = "meanVar.pdf", width = 7, height = 7)
				plotMeanVar(d, show.tagwise.vars = TRUE, NBline = TRUE)
				plotBCV(d)
			dev.off()
	   
	    ### Test for differential expression
			de = exactTest(d, pair = c(condA,condB))
			res = de$table
			
			de.fdr <- p.adjust(as.numeric(res[,"PValue"]))
			res.full <- cbind(res, de.fdr)
			
			colnames(res.full) <- c(colnames(res), "fdr")
				
	  ### Write into a file
		  write.table(res.full, file = "res_EDGER.txt",quote=F,col.names=NA,sep='\t')
	}

######################################################################################################
#### ENSEMBLETOSYMBOL: ENSEMBLE ID TO GENE ID TO GENE SYMBOL
######################################################################################################

    EnsembleToSymbol <- function(Ensemble.ID)
	{  
	  #### EnsembleToSymbol converts Ensemble ID to gene symbol by biomart
	  ##################################################################
	  ### Input:
	  ###       Ensemble.ID: ENSG00000280795
	  ### Output: 
	  ###       Gene.Symbol: MPRIPP1
	  
	  ### note: need to go to headnode to access biomart ssh -X login1
	  ###***************************************************************** 
	  
	  ### remove the . before the Ensemble.ID
		Ensemble.ID.clean <- sapply(Ensemble.ID, function(current.ID)
		{
			current.ID.vec <- unlist(strsplit(current.ID, split = "[.]"))
			return(current.ID.vec[1])
		})
	    
		print(Ensemble.ID.clean[1:5])
	  
	  ### convert the ID to gene symbole
		ensembl_genes <- as.character(Ensemble.ID.clean)
		print(ensembl_genes[1:5])
		
		library(biomaRt) 
		mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
		
		gene.data <- getBM( filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),
	                 values= ensembl_genes, mart= mart)
					 
	  ### further map the IDs to symbol using gene.data as reference
	  ### we need this step since gene.data and ensembl_genes are not of the same dimension
		match.ind <- match(ensembl_genes, as.character(gene.data[, "ensembl_gene_id"]))
		
		ID.symbol <- cbind(Ensemble.ID, gene.data[match.ind, "hgnc_symbol"])
	  
	  return(ID.symbol)
	}

######################################################################################################
### FISHER.OVERLAP PERFORMS FISHER'S EXACT TEST ON OVERLAPPING OF SET1 AND SET2 ON WHOLE.SET
######################################################################################################

	fisher.overlap <- function(set1, set2, whole.set)
	{
	  #### fisher.overlap performs fisher's exact test on overlapping of set1 and set2 on whole.set
	  #### Author: Jialiang Yang  
	  #### Date: 2/24/2015 version: 1
	  ##################################################################
	  ### Input:
	  ###        set1: first set to compare
	  ###        set2: second set to compare
	  ###        whole.set: whole set
	  ### Output:
	  ###        fish.sig: significance of the overlapping between set1 and set2
	  ###***************************************************************** 
	  
	  set1.vec <- rep(0, length(whole.set))
	  set1.ind <- which(match(whole.set, set1) != "NA")
	  
	  if(length(set1.ind) == 0){return(c(1, length(whole.set), 0, 0))}
	  set1.vec[set1.ind] <- 1
	  
	  set2.vec <- rep(0, length(whole.set))
	  set2.ind <- which(match(whole.set, set2) != "NA")
	  if(length(set2.ind) == 0){return(c(1, length(whole.set), 0, 0))}
	  set2.vec[set2.ind] <- 1

	  pvalue <- fisher.test(set1.vec, set2.vec , alternative = "greater")$p.value
	  
	  return(c(pvalue, length(whole.set), length(set1.ind), length(set2.ind)))
	}


	
######################################################################################################
### KDA PERFORMS THE KEY DRIVER ANALYSIS
######################################################################################################	
	
	findNLayerNeighborsLinkPairs <- function( linkpairs , subnetNodes , nlayers = 1 , directed = FALSE )
	{
	  #linkpairs=linkpairs; subnetNodes=ineighbors; nlayers=nlayers-1; directed=directed
	   merged <- merge( linkpairs , subnetNodes , by.x = 1 , by.y = 1 , all = FALSE )
	   merged <- as.matrix( merged )

	   if ( isTRUE( all.equal( dim( merged )[1] , 0 ) ) )
	   {
		   return( NULL )
	   }

	   if ( isTRUE( all.equal( dim( merged )[1] , 1 ) ) &
			isTRUE( all.equal( merged[1,1] , merged[1,2] ) ) )
	   {
		   return( merged )
	   }

	   if ( !directed )
	   {
	# undirected networks
		 mergeleft <- merge( linkpairs , subnetNodes , by.x = 2 , by.y = 1 , all = FALSE )
		 mergeleft <- as.matrix( mergeleft )
		 mergeleft <- mergeleft[,c( 2 , 1 )] # keep the original link direction
		 merged <- rbind( merged , mergeleft )
	   }
	   dim1 <- dim( merged )[1]
	   if ( isTRUE( all.equal( dim1 , 0 ) ) )
	   {
	# no links
			 return( NULL )
	   }else if ( is.null( dim1 ) )
	   {
	# only one link
		   merged <- rbind( merged )
	   }
	   ineighbors <- union( merged[,1] , merged[,2] )

	   if ( isTRUE( all.equal( nlayers , 1 ) ) )
	   {
		   res <- getSubnetworkLinkPairs( linkpairs , subnetNodes = ineighbors )
		   return( res )
	   }

	   # stop earlier if no change
	   #
	   common <- intersect( ineighbors , subnetNodes )
	   if ( length( common ) == length( ineighbors ) )
	   {
		   res <- getSubnetworkLinkPairs( linkpairs , subnetNodes = ineighbors )
		   return( res )
	   }

	   return( findNLayerNeighborsLinkPairs( linkpairs , ineighbors , nlayers - 1 , directed ) )
	}
	
	getSubnetworkLinkPairs <- function( linkpairs , subnetNodes )
	{
	   mergeright <- merge( linkpairs , subnetNodes , by.x = 2 , by.y = 1 , all = FALSE )
	   if ( isTRUE( all.equal( dim( mergeright )[1] , 0 ) ) )
	   {
		 return( NULL )
	   }

	   mergeleft <- merge( mergeright , subnetNodes , by.x = 2 , by.y = 1 , all = FALSE )

	   #mergeright = merge(linkpairs, subnetNodes, by.x=1, by.y=1, all=F)
	   #mergeleft2 = merge(mergeright, subnetNodes, by.x=2, by.y=1, all=F)

	   if ( isTRUE( all.equal( dim( mergeleft )[1] , 0 ) ) )
	   {
		 return( NULL )
	   }

	   return( as.matrix( mergeleft ) )   
	}
	
	downStreamGenes <- function( netpairs , seednodes , N = 100 , directed = TRUE )
	{
	   prenodes <- seednodes
	   cnt <- N
	   pcdiff <- 1
	   while( length( pcdiff ) > 0 && cnt > 0 )
	   {
		  retlinks <- findNLayerNeighborsLinkPairs( linkpairs = netpairs , subnetNodes = prenodes ,
				   nlayers = 1 , directed = directed )
		  if( is.null( retlinks ) )
		  {
			  return( NULL )
		  }
		  curnodes <- union( retlinks[,1] , retlinks[,2] ) 
		  pcdiff <- setdiff( curnodes , prenodes )
		  prenodes <- curnodes
		  cnt <- cnt - 1
	   }

	   if ( is.null( retlinks ) )
	   {
		   return( NULL )
	   }else{
		  return( curnodes )
	   }
	}
	
	setInSets <- function( setC , setlist )
	{
	   for ( i in c( 1:length( setlist ) ) )
	   {
		   isSub <- setsub( setC , setlist[[i]] )
		   if ( isSub )
		   { #setC is a subset of setlist[[i]]
			  return( TRUE )
		   }
	   }
	   return( FALSE )
	}
	
	setsub <- function( setA , setB )
	{
		if ( length( setA ) > length( setB ) )
		{
			return( FALSE )
		}
		setAB <- union( setA , setB )
		return ( setequal( setAB , setB ) )
	}
	
	concatenate <- function( myvect , mysep="" )
	{
	  if ( is.null( myvect ) )
	  {
		return ( "" )
	  }
	  else if ( isTRUE( all.equal( length( myvect ) , 1 ) ) )
	  {
		return ( as.character( myvect ) )
	  }
	  return( paste( as.character( myvect ), sep = "" , collapse = mysep ) )
	}
	
	keyDriverAnalysis <- function( inputnetwork , signature , directed = TRUE , nlayerExpansion = 1 ,
		nlayerSearch = 6 , enrichedNodesPercentCut = -1 , boostHubs = TRUE , dynamicSearch = TRUE ,
		FETpValueCut = 0.05 , useCorrectedpValue = TRUE , outputfile = NULL )
	{
		if ( !is.null( outputfile ) )
		{
			#onetFname <- paste( outputfile , outputDir , key2 , ".pair" , sep = '' )
			snpFname <- paste( outputfile , ".snp" ,  sep = '' )
			kdFname <- paste( outputfile , "_keydriver.xls" , sep = '' )
		}

		# overlap between network & signature
		wholenodes <- union( inputnetwork[,1] , inputnetwork[,2] )
		no.wholenodes <- length( wholenodes )
		wholeOvlp <- intersect( wholenodes , signature )
		no.wholeOvlp <- length( wholeOvlp )

		if ( nlayerExpansion >= 1 )
		{
			# expand network by n-layer nearest neighbors
			expandNet <- findNLayerNeighborsLinkPairs( linkpairs = inputnetwork ,
					subnetNodes = signature , nlayers = nlayerExpansion , directed = directed )
		}
		else if ( isTRUE( all.equal( nlayerExpansion , 0 ) ) )
		{
			# no expansion
			expandNet <- getSubnetworkLinkPairs( linkpairs = inputnetwork ,
												 subnetNodes = signature )
		}
		else
		{
			expandNet <- inputnetwork
		}

		print( paste( "dim(expandNet): " , dim( expandNet ) ) )

		allnodes <- sort( union( expandNet[,1] , expandNet[,2] ) )
		no.nodes <- length( allnodes )
		#write.table(expandNet, onetFname, sep="\t",quote=FALSE, col.names=T, row.names=FALSE)

		# convert IDs into indices
		netIdxSrc <- getMatchedIndexFast( allnodes , expandNet[,1] )
		netIdxDst <- getMatchedIndexFast( allnodes , expandNet[,2] )
		signatIdx <- getMatchedIndexFast( allnodes , intersect( allnodes , signature ) )
		expandNetIdx <- cbind( netIdxSrc , netIdxDst )

		################################################################################################
		# 4. keydriver for a given network
		#
		#linkpairs=expandNet;signature=genes;directed=directed; nlayers=6; min_downstreamnodes=min_ds_cut; FETpValueCut=0.05; boostHubs=T; dynamicSearch=dynamicSearch
		
		if ( directed )
		{
			ret <- keydriverInSubnetwork( linkpairs = expandNetIdx , signature=signatIdx ,
					background = c( no.wholenodes , no.wholeOvlp ) , directed = directed ,
					nlayers = nlayerSearch , enrichedNodesPercentCut = enrichedNodesPercentCut ,
					FETpValueCut = FETpValueCut , boostHubs = boostHubs ,
					dynamicSearch = dynamicSearch , bonferroniCorrection = useCorrectedpValue )
		}
		else
		{
			ret <- keydriverInSubnetwork( linkpairs = expandNetIdx , signature = signatIdx ,
					background = c( no.wholenodes , no.wholeOvlp ) , directed = directed ,
					nlayers = nlayerSearch , enrichedNodesPercentCut = enrichedNodesPercentCut ,
					FETpValueCut = FETpValueCut , boostHubs = boostHubs ,
					dynamicSearch = dynamicSearch , bonferroniCorrection = useCorrectedpValue )
		}

		if ( is.null( ret ) )
		{
			return( NULL )
		}

		# retrieve results
		#
		fkd <- ret[[1]]
		parameters <- ret[[2]]

		print(fkd)
		fkd[,1] <- allnodes[as.integer( fkd[,1] )]   

		if ( !is.null( outputfile ) )
		{
			write.table( fkd , kdFname , sep="\t" , quote = FALSE , col.names = TRUE ,
						 row.names = FALSE )

			################################################################################################
			#  output networks & key drivers for visualization
			#
			#     Cytoscape output: 1) network file - *_cys.txt 2) node property file: *_cys-nodes.txt
			#

			nodeprop <- configureNodeVisualization( allnodes = allnodes , signature = genes ,
													kdaMatrix = fkd )

			hnList <- nodeprop[[1]] # node subcategpries
			listprop <- nodeprop[[2]] # visual properties for each subcategory
			legend <- nodeprop[[3]] # legend table for visual properties

			resf <- makeSNP( netpairsWtype = expandNet , edgecolorlevels = c( "grey" ) ,
					  highlightNodes = hnList , normColor = "grey" , highColor = listprop[,1] ,
					  normShape = "circle" , highShape = listprop[,2] , normNodeSize = "40" ,
					  highNodeSize = listprop[,3] , normFontSize = "12" ,
					  highFontSize = listprop[,4] , legendtable = legend , snafile = snpFname )

			result <- list( expandNet , fkd , ret[[2]] , getFileFullNameNopath( resf ) )
	#		result <- list( expandNet , fkd , ret[[2]] , getFileName( snpFname ) )
			names( result ) <- c( "subnetwork" , "keydrivers" , "parameters" , "files" )
		}
		else
		{
			result <- list( expandNet , fkd , ret[[2]] )
			names( result ) <- c( "subnetwork" , "keydrivers" , "parameters" )
		}

		return( result )
	}
	
	replaceString <- function( fullfnames , oldstr , newstr )
	{
	  no.files <- length( fullfnames )
	  res <- NULL
	  for ( each in fullfnames )
	  {
		#print(paste(i, "/", no.files, ":", each) )
		each2 <- paste( each , oldstr , sep = "" )
		splitted <- splitString( each2 , oldstr )

		neweach <- concatenate( splitted , newstr )

	# How could this ever be run?  There is no way for this conditional to evaluate to TRUE!
	#    if ( FALSE )
	#	{
	#      neweach <- ""
	#      for ( is in splitted )
	#	  {
	#        neweach <- paste( neweach , newstr , sep = is )
	#      }
	#    }

		#oldeach  = paste(pathnet, each,   sep="")
		#neweach  = paste(pathnet, newstr, splitted[2], sep="")
		#a=file.rename(from=oldeach, to=neweach)
		#print(a)

		res <- c( res , neweach )
	  }
	  return( res )
	}
	
	removeDuplicatedLinks <- function( linkpairs , directed = FALSE )
	{
		if ( isTRUE( all.equal( dim( linkpairs )[1] , 1 ) ) )
		{
		   return( linkpairs )
		}
		links <- paste( linkpairs[,1] , linkpairs[,2] , sep = "\t" )
		# 1. remove duplications 
		#
		cleanedlinkMatrix <- union( links , NULL )
		# 2. remove inversed duplications
		#
		linkMatrix  <- as.matrix( getAllParts( cleanedlinkMatrix , "\t" ) )
		if ( directed || isTRUE( all.equal( dim( linkMatrix )[1] , 1 ) ) )
		{
			return( linkMatrix )
		}
		#  first, remove self-self interactions
		#
		selfSelfLinks <- linkMatrix[,1] == linkMatrix[,2]
		linkMatrix <- linkMatrix[!selfSelfLinks,]
		cleanedlinkMatrix <- cleanedlinkMatrix[!selfSelfLinks]
		# Now, create reverse links
		reversedLinks <- cbind( paste( linkMatrix[,2] , linkMatrix[,1] , sep = "\t" ) , c( 1:length( cleanedlinkMatrix ) ) )
		removedCols <- as.integer( merge( cleanedlinkMatrix , reversedLinks , by.x = 1 , by.y = 1 , all = FALSE )[,2] )
		if ( length( removedCols ) > 0 )
		{
			# construct non-duplicated interactions
			#
			dupLinks <- cleanedlinkMatrix[removedCols]
			dupLinksRev <- reversedLinks[removedCols]
			uniques <- NULL
			for ( i in c( 1:length( dupLinks ) ) )
			{
			   if ( !( is.element( dupLinks[i] , uniques ) | is.element( dupLinks[i] , uniques ) ) )
			   {
				   uniques <- c( uniques , dupLinks[i] )
			   }
			}
			xlinkMatrix <- c( cleanedlinkMatrix[-removedCols] , uniques )
		}
		else
		{
			xlinkMatrix <- cleanedlinkMatrix
		}
		return( getAllParts( xlinkMatrix , "\t" ) )
	}

	
	mergeTwoMatricesByKeepAllPrimary <- function( primaryMatrix , minorMatrix , missinglabel = "" ,
		                            keepAllPrimary = TRUE , keepPrimaryOrder = TRUE ,
									keepAll = FALSE )
	{
	  no.promarycols <- dim( primaryMatrix )[2]
	  no.mustbegenes <- dim( primaryMatrix )[1]

	  # we add in one more column to indicate which genes are mustbeincluded after being merged with mcg
	  keyword <- "mustbeused"
	  mustbeGenesMatrix <- cbind( primaryMatrix , c( 1:no.mustbegenes ) , rep( keyword , no.mustbegenes ) )

	  if ( is.null( colnames( primaryMatrix ) ) )
	  {
		colnames( mustbeGenesMatrix ) <- c( c( 1:no.promarycols ) , "primorder" , keyword )
	  }
	  else
	  {
		colnames( mustbeGenesMatrix ) <- c( colnames( primaryMatrix ) , "primorder" , keyword )
	  }
	# Why is this uncommented?
	#  dim( mustbeGenesMatrix )

	  if ( is.null( keepAllPrimary ) )
	  { #normal merge: to have the common elements
		myMatrix <- merge( mustbeGenesMatrix , minorMatrix , by.x = 1 , by.y = 1 , all.x = FALSE ,
						   sort = FALSE , all = FALSE )
	  }
	  else
	  {
		myMatrix <- merge( mustbeGenesMatrix , minorMatrix , by.x = 1 , by.y = 1 , all.x = TRUE ,
						   sort = FALSE , all = TRUE )
	  }
	# Again, why is this left uncommented?
	#  dim( myMatrix )
	  nocols.mymatrix <- dim( myMatrix )[2]

	  #the mustbeused genes which are not included in minor have NAs in the column $mustbeused
	  #so we can use this information to figure out which mustbeused genes missing in minorMatrix
	  myMatrix[,nocols.mymatrix] <- ifelse( is.na( myMatrix[,nocols.mymatrix] ) , missinglabel ,
									as.character( myMatrix[,nocols.mymatrix] ) )

	  orders <- order( as.numeric( as.matrix( myMatrix[,no.promarycols + 1] ) ) )
	  if ( keepPrimaryOrder )
		  myMatrix <- myMatrix[orders,]

	  if ( is.null( keepAllPrimary ) )
	  {
		 selected <- rep( T , dim( myMatrix )[1] )
	  }
	  else
	  {
		 if ( keepAllPrimary )
		 {
		   selected <- !( is.na( myMatrix[,no.promarycols + 2] ) )
		 }
		 else #return the row-elements in minor which are missed in primary
		 {
		   selected <- is.na( myMatrix[,no.promarycols + 2] )
		 }
	  }
	# Why are these here?
	#  sum(selected)

	  #keep the primary matrix and remove the mustbeused column
	  return( myMatrix[selected,-c( no.promarycols + 1 , no.promarycols + 2 )] )
	}
	
	configureNodeVisualization <- function( allnodes , signature , kdaMatrix , bNodeSz = 40 , bFontSz = 12 )
	{
		
		# SIG--signature; NSIG--not signature; GKD--Global KeyDriver; LKD--Local KeyDriver; NKD--Not KeyDriver
		#
		xcategories <- c( "SIG_GKD" , "SIG_LKD" , "SIG_NKD" , "NSIG_GKD" , "NSIG_LKD" , "NSIG_NKD" )
		xcolors <- c( "red" , "blue" , "lightgreen" , "red" , "blue" , "grey" )
		names( xcolors ) <- xcategories
		xshapes <- c( "square" , "square" , "circle" , "circle" , "circle", "circle" )
		names( xshapes ) <- xcategories
		xsizes <- c( 3 * bNodeSz , 2 * bNodeSz , bNodeSz , 3 * bNodeSz , 2 * bNodeSz , bNodeSz )
		names(xsizes) <- xcategories
		xfontsz <- c( 3 * bFontSz , 2 * bFontSz , bFontSz , 3 * bFontSz , 2 * bFontSz , bFontSz )
		names( xfontsz ) <- xcategories
		
		no.nodes <- length( allnodes )
		
		# legend table 
		legendtb <- cbind( xcategories , xshapes , xcolors , xcolors , xsizes , xfontsz )
		colnames( legendtb ) <- c( "label" , "shape" , "color" , "border" , "node_size" , "font_size" )
		
		sigInNet <- intersect( allnodes , signature )
		sigStatus <- rep( "NSIG" , no.nodes )
		names( sigStatus ) <- allnodes
		sigStatus[sigInNet] <- "SIG"
		kdrStatus <- rep( "NKD" ,  no.nodes )
		names( kdrStatus ) <- allnodes

		nf.cols <- dim( kdaMatrix )[2]
		nf.rows <- dim( kdaMatrix )[1]
		keydrvNames <- NULL
		if ( nf.rows > 0 )
		{
			keydrv <- as.integer( kdaMatrix[,nf.cols] )
			# global driver
			keysel <- which( keydrv == 1 )
			keydrvNames <- kdaMatrix[keysel,1]
			kdrStatus[keydrvNames] <- "GKD"

			# local driver
			if ( any( keydrv == 0 ) )
			{
				keysel <- which( keydrv == 0 )
				keydrvNames <- kdaMatrix[keysel,1]
				kdrStatus[keydrvNames] <- "LKD"
			}

			# combined signature-keydriver status
			#
			sigkdrStatus <- paste( sigStatus , kdrStatus , sep = "_" )
			hnList <- tapply( allnodes , sigkdrStatus , list ) # make a list for each category
			sigkdrNames <- names( hnList )

			isNonSig <- intersect( xcategories[4:6] , sigkdrNames ) # if all nodes are signatures, we use only circle for display
			if ( isTRUE( all.equal( length( isNonSig ) , 0 ) ) )
			{
				xshapes <- c( "circle" , "circle" , "circle" , "circle" , "circle" , "circle" )
				names( xshapes ) <- xcategories
			}

			# set up actual visualization properties
			yHighColor <- xcolors[sigkdrNames]
			yHighShape <- xshapes[sigkdrNames]
			yHighSize <- xsizes[sigkdrNames]
			yHighFontSZ <- xfontsz[sigkdrNames]
		}else{
			hnList <- list( sigInNet ) # highlight only signature
			yHighColor <- c( "brown" )
			yHighShape <- c( "circle" )
			yHighSize <- c( "1" )
			yHighFontSZ <- c( "1" )
		}
		return( list( hnList , cbind( yHighColor , yHighShape , yHighSize , yHighFontSZ ) , legendtb ) )
	}
	
	degreeByLinkPairs <- function( linkpairs , directed = FALSE , cleangarbage = FALSE )
	{
		codepair <- c( 0 , 1 )  #[1] for no connection, [2] for connection

		edgesInNet <- dim( linkpairs )[1]
		
		# consider both columns
		# Need to find out why we build the object this way
		#  Could be a cleaner way, but there may be a reason why
		#  it's currently done this way
		allnodenames <- NULL
		allnodenames <- c( allnodenames , as.character( linkpairs[,1] ) )
		allnodenames <- c( allnodenames , as.character( linkpairs[,2] ) )
		 
		nametable <- table( allnodenames )
		# Not sure why this was here, uncommented out
		# length( nametable )
		
		uniquenames <- names( nametable )
		no.uniquenames <- length( uniquenames )

		totallinks <- as.integer( nametable ) # no of links for each node
		totalmatrix <- cbind( names( nametable ) ,  totallinks )

		if ( directed )
		{
			# outlines
			dnodenames <- as.character( linkpairs[,1] )
			dnametable <- table( dnodenames )
			duniquenames <- names( dnametable )
			dmatrix <- cbind( names( dnametable ) , as.integer( dnametable ) )
			colnames( dmatrix ) <- c( "node" , "links" )

			iolinks <- mergeTwoMatricesByKeepAllPrimary( primaryMatrix = cbind( uniquenames ) ,
					minorMatrix = dmatrix , missinglabel = "0" , keepAllPrimary = TRUE ,
					keepPrimaryOrder = TRUE , keepAll = FALSE )
			outlinks <- as.integer( as.matrix( iolinks[,2] ) )


			# inlines
			dnodenames <- as.character( linkpairs[,2] )
			dnametable <- table( dnodenames )
			duniquenames <- names( dnametable )
			dmatrix <- cbind( names( dnametable ) , as.integer( dnametable ) )
			colnames( dmatrix ) <- c( "node" , "links" )

			iolinks <- mergeTwoMatricesByKeepAllPrimary( primaryMatrix = cbind( uniquenames ) ,
					minorMatrix = dmatrix , missinglabel = "0" , keepAllPrimary = TRUE ,
					keepPrimaryOrder = TRUE , keepAll = FALSE )
			inlinks <- as.integer( as.matrix( iolinks[,2] ) ) 

		}else{
			inlinks <- totallinks
			outlinks <- totallinks
		}

		#hubidx    = order(-totallinks)

		# output in/out links for each gene
		#
		linksMatrix <- cbind( inlinks , outlinks , totallinks )
		colnames( linksMatrix ) <- c( "inlinks" , "outlinks" , "totallinks" )
		rownames( linksMatrix ) <- uniquenames

		rm( inlinks , outlinks , totallinks )

		if ( cleangarbage )
		{
			collect_garbage()
		}

		return( data.frame( linksMatrix ) )
	}
	
	getAllParts <- function( fullfnames , sep = "-" , retLen = FALSE )
	{
	  unSplit <- function( x ) { unlist( strsplit( x , sep ) ) }
	  ReturnMatrix <- function( x , returnMatrixEnv , unSplit )
	  {
		  usplit <- unSplit( x )
		  returnMatrixEnv$returnMatrixList[[returnMatrixEnv$i]] <- usplit
		  returnMatrixEnv$i <- returnMatrixEnv$i + 1
		  return( length( usplit ) )
	  }
	  makeMatrix <- function( x , returnMatrixEnv , rowLengths )
	  {
		  returnMatrixEnv$returnMatrix[returnMatrixEnv$i,1:rowLengths[returnMatrixEnv$i]] <- x
		  returnMatrixEnv$i <- returnMatrixEnv$i + 1
		  return( rowLengths[returnMatrixEnv$i-1] )
	  }
	  returnLengths <- function( x , unSplit ) { return( length( unSplit( x ) ) ) }
	  if ( retLen )
	  {
		return( unlist( lapply( fullfnames , returnLengths , unSplit ) ) )
	  }else{
		returnMatrixEnv <- new.env()
		returnMatrixEnv$returnMatrixList <- vector( "list" , length = length( fullfnames ) )
		returnMatrixEnv$i <- 1
		rowLengths <- unlist( lapply( fullfnames , ReturnMatrix , returnMatrixEnv , unSplit ) )
		returnMatrixEnv$returnMatrix <- matrix( "" , nrow = length( rowLengths ) , ncol = rowLengths[which.max( rowLengths )] )
		returnMatrixEnv$i <- 1
		rowLengths <- unlist( lapply( returnMatrixEnv$returnMatrixList , makeMatrix , returnMatrixEnv , rowLengths ) )
		return( returnMatrixEnv$returnMatrix )
	  }
	}
	
	getFileExtension <- function( fullfname )
	{
		splitted <- unlist( strsplit( fullfname , "\\." ) )

		if ( length( splitted ) > 1 )
		{
		  return( splitted[length( splitted )] )
		}else{
		  return( "" )
		}
	}

	#get the filename without path information
	getFileFullNameNopath <- function( fullfnames )
	{
		res <- NULL
		for ( each in fullfnames )
		{
			splitted <- unlist( strsplit( each , "/" ) )
			res <- c( res , splitted[length( splitted )] )
		}
		return( res )
	}
	
	getFileName <- function( fullfname )
	{
		ext <- getFileExtension( fullfname )
		if (ext == "" )
		{
		   return( fullfname )
		}
		extd <- paste( "." , ext , sep = "" )
		return( splitString( fullfname , extd )[1] )
	}
	
	getMatchedIndexFast <- function( cvector , subvect )
	{
		orgIndex <- cbind( cvector , c( 1:length( cvector ) ) )
		subIndex <- cbind( subvect , c( 1:length( subvect ) ) )
		merged <- as.data.frame( merge( subIndex , orgIndex , by.x = 1 , by.y = 1 , all.x = TRUE ) )
		if ( dim( merged )[1] > 1 )
		{
			od <- order( merged[,2] )  # restore the original order of subvect
			merged <- merged[od,]
		}
		return( merged[,3] )
	}

		
	keydriverInSubnetwork <- function( linkpairs, signature, background = NULL, directed = TRUE, nlayers = 6, 
	     enrichedNodesPercentCut = -1, FETpValueCut = 0.05, boostHubs = TRUE, dynamicSearch = TRUE, bonferroniCorrection = TRUE)
	{
		allnodes <- union( linkpairs[,1] , linkpairs[,2] )
		no.subnetsize <- length( allnodes )

		# whole network nodes as the signature
		networkAsSignature <- isTRUE( all.equal( length( setdiff( allnodes , signature ) ) , 0 ) )
		
		overlapped <- intersect( allnodes , signature )
		no.overlapped <- length( overlapped ) # within the subnetwork

		if ( is.null( background ) )
		{
			background2 <- c( no.subnetsize , no.overlapped ) 
		}else{
			background2 <- background
		}

		keydrivers <- NULL
		kdMatrix <- NULL
		kdIndex <- NULL # indices of keydrivers in dsnodesList

		dsnodesList <- as.list( rep( 0 , no.subnetsize ) )
		no.dsnodes <- rep( 0 , no.subnetsize )
		cnt <- 1

		intv <- as.integer( no.subnetsize / 10 )
		print( "find downstream genes" )

		# set up searching range
		if ( dynamicSearch )
		{ # dynamic search for optimal layer
			layers4Search <- c( 1:nlayers )
		}else{  # fixed layer for optimal layer
			layers4Search <- c( nlayers )
		}
		# if the network itself is the signature, no need for dynamic search
		if ( networkAsSignature )
		{
			layers4Search <- c( nlayers )
		}

		for ( i in c( 1:no.subnetsize ) )
		{
			if ( isTRUE( all.equal( i%%intv , 0 ) ) )
			{
				print( paste( i , "/" , no.subnetsize ) )
			}

			# initialization
			minpv <- 1
			minNoHits <- 0
			minNoIdn <- 0
			minLayer <- 0
			minDn <- 0
			minFc <- 0
			minpvW <- 1
			minFcW <- 0
			for ( y in layers4Search )
			{
				#netpairs=linkpairs; seednodes=allnodes[i]; N=nlayers; directed=directed
				idn <- downStreamGenes( netpairs = linkpairs , seednodes = allnodes[i] ,
						N = y , directed = directed )
				idn <- setdiff(idn, allnodes[i])
				no.idn <- length(idn)

				if ( !networkAsSignature )
				{# do enrichment test for only expanded subnetwork
					hits <- intersect( idn , overlapped )
					no.hits <- length( hits )

					if ( isTRUE( all.equal( no.hits , 0 ) ) )
					{
						next
					}

					foldchg <- ( no.hits / no.idn ) / ( no.overlapped / no.subnetsize )
					pv <- phyper( no.hits - 1 , no.idn , no.subnetsize - no.idn ,
							no.overlapped , lower.tail = FALSE )

					foldchgW <- ( no.hits / no.idn ) / ( background2[2] / background2[1] )
					pvW <- phyper( no.hits - 1 , no.idn , background2[1] - no.idn ,
							background2[2] , lower.tail = FALSE )

					if ( pv < minpv )
					{
						minpv <- pv
						minNoHits <- no.hits
						minNoIdn <- no.idn
						minLayer <- y
						minFc <- foldchg
						minpvW <- pvW
						minFcW <- foldchgW
					}
				}else{ # for non-expanded subnetwork
					no.hits <- no.idn
					minpv <- 0
					minNoHits <- no.idn
					minNoIdn <- no.idn
					minLayer <- y
					minFc <- 1
				}
			} #y
			
			# record the down stream genes for the biggest layer
			if ( no.idn > 0 )
			{
				dsnodesList[[i]] <- idn
				no.dsnodes[i] <- no.idn
			}

			correctMinPv <- minpv * no.subnetsize
			correctMinPv <- ifelse( correctMinPv > 1 , 1 , correctMinPv )
			
			res <- c( minNoHits , minNoIdn , no.overlapped , no.subnetsize , background2[2] ,
					background2[1] , length( signature ) , minLayer , minFcW , minpvW , minFc ,
					minpv , correctMinPv )
			kdMatrix <- rbind( kdMatrix , res )
			#print(res)
		}
		
		mymincut <- enrichedNodesPercentCut * no.overlapped
		if ( enrichedNodesPercentCut <= 0 )
		{
			mymincut = mean(no.dsnodes) + sd(no.dsnodes)
		}
		cutmatrix <- c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut )

		# pick up key drivers by pvalue and no. of downstream genes
		ncols <- dim( kdMatrix )[2]
		
		if ( bonferroniCorrection )
		{ # use corrected pvalue
			kdSel <- ( kdMatrix[,ncols] < FETpValueCut ) & ( kdMatrix[,2] >= mymincut )
		}else{
			kdSel <- ( kdMatrix[,ncols-1] < FETpValueCut) & (kdMatrix[,2] >= mymincut )
		}

		if ( sum( kdSel ) > 0 )
		{
			keydrivers <- allnodes[kdSel]
			kdIndex <- c( 1:no.subnetsize )[kdSel]
			n.drivers <- length( keydrivers )

			#******************* local driver or not **************************************
			#
			# check whether a driver is in the downstream of other drivers
			keydrv <- rep( 0 , no.subnetsize )
			#if (!networkAsSignature) {
			for ( i in c( 1:n.drivers ) )
			{
				# Note that kdIndex[i] is the index of ith keydriver in kdMatrix  
				# restrict to only candidate drivers 
				iselA <- ( kdMatrix[,2] > kdMatrix[kdIndex[i],2] ) & kdSel
				isel <- c( 1:no.subnetsize )[iselA]
				if ( sum(isel) > 0 )
				{
					if ( directed )
					{
						ilocal <- setInSets( setC = allnodes[kdIndex[i]] ,
								setlist = dsnodesList[isel] )
					}else{
						ilocal <- setInSets( setC = dsnodesList[[kdIndex[i]]] ,
								setlist = dsnodesList[isel] )
					}
					keydrv[kdIndex[i]] <- !ilocal + 0
				}else{
					keydrv[kdIndex[i]] <- TRUE
				}
			}
			#}
		}
		
		# promote genes with many direct links to be key drivers
		#
		#              inlinks outlinks totallinks
		#0610031J06Rik       2        0          2
		#1110001J03Rik       0        1          1
		#
		if ( boostHubs )
		{
			
			if ( !networkAsSignature )
			{
				# for expanded network, restrict the boosted nodes to the key driver candidates
				kdSelB <- rep( FALSE , no.subnetsize )
				kdSelB[kdIndex] <- TRUE
				psel <- kdMatrix[,ncols-3] * no.subnetsize < 0.05
				kdPsd <-1
				kdPmean <- 1
				kdpvalues <- -log10( kdMatrix[,ncols-3] )
				kdpvalues <- ifelse( is.na( kdpvalues ) , 0 , kdpvalues )
				#histogram(kdpvalues)
				if ( sum( psel ) > 0 )
				{
					kdPmean <- mean( kdpvalues[psel] )
					kdPsd <- sd( kdpvalues[psel] )
					#kdPmean= median(kdpvalues[psel]); kdPsd= mad(kdpvalues[psel])
					print( as.numeric( signif( kdpvalues[psel] , 2 ) ) )
					directSel <- ( kdpvalues > ( kdPmean + kdPsd ) )
					directSel <- ifelse( is.na( directSel ) , FALSE , directSel )
					#print(directSel)
					if ( sum( directSel ) > 0 )
					{
						kdSel <- kdSel | directSel
						dIndex <- c( 1:no.subnetsize )[kdSel]
						keydrv <- rep( FALSE , no.subnetsize )
						keydrv[dIndex] <- TRUE
					}
				}
				cutmatrix <- rbind( c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut ,
								kdPmean , kdPsd , kdPmean + kdPsd ) )

				colnames( cutmatrix ) <- c( "mean_downstream" , "sd_downstream" ,
						"enrichedNodes_cut" , "mean_logP" , "sd_logP" , "cut_logP" )
			}else{
				# for non-expanded network, consider all the nodes in the subnetwork
				kdSelB <- rep( TRUE , no.subnetsize )
				# align the degree with allnodes
				mydegree <- degreeByLinkPairs( linkpairs = linkpairs , directed = directed ,
						cleangarbage = FALSE )
				mIdx <- getMatchedIndexFast( rownames( mydegree ) , allnodes )
				mydegree <- mydegree[mIdx,]
				if ( directed )
				{
					directSel <- mydegree[,2] > mean( mydegree[,2] ) + 2 * sd( mydegree[,2] )
					cutmatrix <- rbind( c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut ,
									mean( mydegree[,2] ) , sd( mydegree[,2] ) ,
									mean( mydegree[,2] ) + 2 * sd( mydegree[,2] ) ) )
				}else{
					directSel <- mydegree[,3] > mean( mydegree[,3] ) + 2 * sd( mydegree[,3] )
					cutmatrix <- rbind( c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut ,
									mean( mydegree[,3] ) , sd( mydegree[,3] ) ,
									mean( mydegree[,3] ) + 2 * sd( mydegree[,3] ) ) )
				}
				directSel <- directSel & kdSelB

				directeHub <- rownames( mydegree )[directSel]
				isDirectHub <- setElementInSet( allnodes , directeHub )

				keydrv[isDirectHub] <- TRUE
				kdSel <- kdSel | isDirectHub
				colnames( cutmatrix ) <- c( "mean_downstream" , "sd_downstream" ,
						"cut_downstream" , "mean_degree" , "sd_degree" , "cut_degree" )
			}
		}else{
			cutmatrix <- rbind( c( mean( no.dsnodes ) , sd( no.dsnodes ) , mymincut , "F" ) )
			colnames( cutmatrix ) <- c( "mean_downstream" , "sd_downstream" , "cut_downstream" ,
					"boost_directhubs" )
		}
		
		if ( isTRUE( all.equal( sum( kdSel ) , 0 ) ) )
		{
			return( NULL )
		}
		
		##
		# in this case, signature is the network nodes themselves, so pvalue will be 0 for all nodes
		# so the driver will be the ones with most downstream genes
		#
		isSignature <- rep( 0 , no.subnetsize )
		names( isSignature ) <- allnodes
		isSignature[overlapped] <- 1

		fkd <- cbind( allnodes , isSignature , kdMatrix , keydrv + 0 )[kdSel,]

		if ( sum( kdSel ) > 1 )
		{
			nf.cols <- dim( fkd )[2]
			if ( networkAsSignature )
			{
				mo <- order( -as.integer( fkd[,3] ) )
			}else{
				mo <- order( as.numeric( fkd[,nf.cols - 1] ) )
			}

			fkd <- fkd[mo,]
			# put key driver on the top
			mo <- order( -as.integer( fkd[,nf.cols] ) )
			fkd <- fkd[mo,]
		}else{
			fkd <- rbind( fkd )
		}
		
		colnames( fkd ) <- c( "keydrivers" , "isSignature" , "hits" , "downstream" ,
				"signature_in_subnetwork" , "subnetwork_size" , "signature_in_network" ,
				"network_size" , "signature" , "optimal_layer" , "fold_change_whole" ,
				"pvalue_whole" , "fold_change_subnet" , "pvalue_subnet" ,
				"pvalue_corrected_subnet" , "keydriver" )

		print( fkd )

		ret <- as.list( c( 1:2 ) )
		ret[[1]] <- fkd
		ret[[2]] <- cutmatrix

		return( ret )
	}
	
	configureNodeVisualization <- function( allnodes , signature , kdaMatrix , bNodeSz = 40 , bFontSz = 12 )
	{
		
		# SIG--signature; NSIG--not signature; GKD--Global KeyDriver; LKD--Local KeyDriver; NKD--Not KeyDriver
		#
		xcategories <- c( "SIG_GKD" , "SIG_LKD" , "SIG_NKD" , "NSIG_GKD" , "NSIG_LKD" , "NSIG_NKD" )
		xcolors <- c( "red" , "blue" , "lightgreen" , "red" , "blue" , "grey" )
		names( xcolors ) <- xcategories
		xshapes <- c( "square" , "square" , "circle" , "circle" , "circle", "circle" )
		names( xshapes ) <- xcategories
		xsizes <- c( 3 * bNodeSz , 2 * bNodeSz , bNodeSz , 3 * bNodeSz , 2 * bNodeSz , bNodeSz )
		names(xsizes) <- xcategories
		xfontsz <- c( 3 * bFontSz , 2 * bFontSz , bFontSz , 3 * bFontSz , 2 * bFontSz , bFontSz )
		names( xfontsz ) <- xcategories
		
		no.nodes <- length( allnodes )
		
		# legend table 
		legendtb <- cbind( xcategories , xshapes , xcolors , xcolors , xsizes , xfontsz )
		colnames( legendtb ) <- c( "label" , "shape" , "color" , "border" , "node_size" , "font_size" )
		
		sigInNet <- intersect( allnodes , signature )
		sigStatus <- rep( "NSIG" , no.nodes )
		names( sigStatus ) <- allnodes
		sigStatus[sigInNet] <- "SIG"
		kdrStatus <- rep( "NKD" ,  no.nodes )
		names( kdrStatus ) <- allnodes

		nf.cols <- dim( kdaMatrix )[2]
		nf.rows <- dim( kdaMatrix )[1]
		keydrvNames <- NULL
		if ( nf.rows > 0 )
		{
			keydrv <- as.integer( kdaMatrix[,nf.cols] )
			# global driver
			keysel <- which( keydrv == 1 )
			keydrvNames <- kdaMatrix[keysel,1]
			kdrStatus[keydrvNames] <- "GKD"

			# local driver
			if ( any( keydrv == 0 ) )
			{
				keysel <- which( keydrv == 0 )
				keydrvNames <- kdaMatrix[keysel,1]
				kdrStatus[keydrvNames] <- "LKD"
			}

			# combined signature-keydriver status
			#
			sigkdrStatus <- paste( sigStatus , kdrStatus , sep = "_" )
			hnList <- tapply( allnodes , sigkdrStatus , list ) # make a list for each category
			sigkdrNames <- names( hnList )

			isNonSig <- intersect( xcategories[4:6] , sigkdrNames ) # if all nodes are signatures, we use only circle for display
			if ( isTRUE( all.equal( length( isNonSig ) , 0 ) ) )
			{
				xshapes <- c( "circle" , "circle" , "circle" , "circle" , "circle" , "circle" )
				names( xshapes ) <- xcategories
			}

			# set up actual visualization properties
			yHighColor <- xcolors[sigkdrNames]
			yHighShape <- xshapes[sigkdrNames]
			yHighSize <- xsizes[sigkdrNames]
			yHighFontSZ <- xfontsz[sigkdrNames]
		}
		else
		{
			hnList <- list( sigInNet ) # highlight only signature
			yHighColor <- c( "brown" )
			yHighShape <- c( "circle" )
			yHighSize <- c( "1" )
			yHighFontSZ <- c( "1" )
		}
		return( list( hnList , cbind( yHighColor , yHighShape , yHighSize , yHighFontSZ ) , legendtb ) )
	}
	
	setElementInSet <- function( setA , setB )
	{
		found <- rep( FALSE , length( setA ) )
		for ( i in c( 1:length( setA ) ) )
		{
		   idiff <- setdiff( setA[i] , setB )
		   found[i] <- isTRUE( all.equal( length( idiff ) , 0 ) )
		}
		return ( found )
	}

	KDA<- function(fname,genes,layer,cnet,directed=F)
	{
	  ##### KDA performs the key driver analysis
	  ####################################################################
	  ### Input:
	  
	  ### Output:
	  #################################################################### 
	  
		  outputDir<- paste0(fname,'_layer',layer,'/')
		  fname.vec <- unlist(strsplit(fname, split = "/"))
		  fname <- fname.vec[length(fname.vec)]
		  
		  if(layer >=1){
			expandNet<- findNLayerNeighborsLinkPairs(linkpairs = cnet , subnetNodes = genes , nlayers = layer , directed = FALSE)
		   } else {
		   expandNet<- getSubnetworkLinkPairs(linkpairs = cnet , subnetNodes = genes)
		   }
		   
		   print(dim(expandNet))
		   
		   allnodes<- union(expandNet[,1] , expandNet[,2])
		   
		   if (directed){
			ret<- keydriverInSubnetwork(linkpairs = expandNet , signature = genes, background=NULL, directed = directed , nlayers = 6 , enrichedNodesPercentCut=-1, 
			FETpValueCut=0.05, boostHubs=T, dynamicSearc=T, bonferroniCorrection=T)
		   } else {
			 ret<- keydriverInSubnetwork(linkpairs = expandNet, signature = genes, directed = directed, nlayers = 2, enrichedNodesPercentCut=-1, FETpValueCut=0.05,
			 boostHubs=T, dynamicSearch=T, bonferroniCorrection=T)
			}
			
			if (!is.null(ret)) 
			{
			   dir.create(outputDir)
			   
			   fkd<- ret[[1]]
			   
			   write.table(fkd , file=paste(outputDir, fname , '_keydriver.xls' , sep = '') , sep = '\t' , quote = FALSE , col.names = TRUE , row.names = FALSE)
			  
			   paraMatrix<- ret[[2]]	   
			   write.table(paraMatrix , file=paste(outputDir, fname , '_parameters.xls' , sep = '') , sep = '\t' , quote = FALSE , col.names = TRUE , row.names = FALSE)
			
			  nodeprop = configureNodeVisualization(allnodes=allnodes, signature=genes, kdaMatrix=fkd)
			  
			  hnList     = nodeprop[[1]] # node subcategpries
			  listprop   = nodeprop[[2]] # visual properties for each subcategory
			  legend     = nodeprop[[3]] # legend table for visual propertie
			  
			  suppressWarnings(makeSNP(netpairsWtype   = expandNet, 
			  edgecolorlevels = c('grey'),
			  highlightNodes  = hnList,
			  normColor='grey',   highColor=listprop[,1],
			  normShape='circle', highShape=listprop[,2],
			  normNodeSize ='40',  highNodeSize =listprop[,3],
			  normFontSize ='12',  highFontSize =listprop[,4],
			  legendtable=legend, snafile=paste(outputDir, fname , '_net.txt' , sep = '')))
			  }
	}
	
	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SNP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# the third column of netpairsWtype is the link type 
#  
# highlightNodes can be a vector of genes or a list of gene sets with
#    highColor and highShape for each gene set, i.e.,
#
#  highlightNodes[[i]] <-  highColor[i] and highShape[i]
#
	makeSNP <- function( netpairsWtype , highlightNodes = NULL , edgecolorlevels ,
			normColor = "grey" , highColor = "red" , normShape = "50" , highShape = "50" , 
			normNodeSize = "1" , highNodeSize = "2" , normFontSize = "1" , highFontSize = "2" ,
			directed = TRUE , legendtable = NA , snafile = "tmp.sna" )
	{
		xfname <- getFileName( snafile )
		fcys <- paste( xfname , "_cys.txt" , sep = "" )
		fcysn <- paste( xfname , "_cys-nodes.txt" , sep = "" )
		write.table( netpairsWtype , fcys , sep = "\t" , quote = FALSE , col.names = TRUE ,
				row.names = FALSE )
		pcols <- dim( netpairsWtype )[2]
		# consider both columns
		uniquenames <- union( as.character( netpairsWtype[,1] ) , as.character( netpairsWtype[,2] ) )
		uniquenames <- sort( uniquenames )
		no.uniquenames <- length( uniquenames )

		# initialize matrix with the size equal to no.uniquenames
		#
		name2idxMatrix <- cbind( uniquenames , c( 1:no.uniquenames ) )
		# 0. make link index: A1 A2 T & A I
		# A1 A2 T & A I ==> A1 A2 T I1
		#
		leftIdx <- merge( netpairsWtype , name2idxMatrix , by.x = 1 , by.y = 1 , all = FALSE ) 
		#  A1 A2 T I1 & A I ==> A2 A1 T I1 I2: I1 and I2 are the indices of A1 and A2 respectively
		#
		allIdx <- merge( leftIdx , name2idxMatrix , by.x = 2 , by.y = 1 , all = FALSE )
		no.pairs <- dim( allIdx )[1]
		if ( isTRUE( all.equal( pcols , 2 ) ) )
		{
			no.elevels <- length( edgecolorlevels )
			greyIdx <- c( 1:no.elevels )[edgecolorlevels=="grey"]
			linksIdxMatrix <- cbind( allIdx[,c( 3 , 4 )] , rep( greyIdx , no.pairs ) )
		}
		else
		{
			linksIdxMatrix <- allIdx[,c(4,5,3)]
		}
		# 1. head string
		header <- paste( "*vertices " , as.character( no.uniquenames ) , sep = "" )
		# 2. vertices matrix 
		if ( is.null( highlightNodes ) )
		{
			verticesMatrix <- cbind( uniquenames ,
					rep( normColor , no.uniquenames ) ,
					rep( normShape , no.uniquenames ) ,
					rep( normNodeSize , no.uniquenames ) ,
					rep( normFontSize , no.uniquenames ) )
		}
		else
		{
			verticesMatrix <- matrix( "" , no.uniquenames , 5 )
			verticesMatrix[,1] <- uniquenames
			verticesMatrix[,2] <- rep( normColor , no.uniquenames )
			verticesMatrix[,3] <- rep( normShape , no.uniquenames )
			verticesMatrix[,4] <- rep( normNodeSize , no.uniquenames )
			verticesMatrix[,5] <- rep( normFontSize , no.uniquenames )
			xcolor <- rep( normColor , no.uniquenames )
			names( xcolor ) <- uniquenames
			xshape <- rep( normShape , no.uniquenames )
			names( xshape ) <- uniquenames
			xNsize <- rep( normNodeSize , no.uniquenames )
			names( xNsize ) <- uniquenames
			xSsize <- rep( normFontSize , no.uniquenames )
			names( xSsize ) <- uniquenames
			# set color and shape for highlighted nodes
			#
			if ( !is.list( highlightNodes ) )
			{
				highlightNodes2 <- intersect( highlightNodes , uniquenames )
				xcolor[highlightNodes2] <- highColor[1]
				xshape[highlightNodes2] <- highShape[1]
				xNsize[highlightNodes2] <- highNodeSize[1]
				xSsize[highlightNodes2] <- highFontSize[1]
			}
			else
			{
				no.highnodeSets <- length( highlightNodes )
				for ( il in c( 1:no.highnodeSets ) )
				{
					highlightNodes2 <- intersect( highlightNodes[[il]] , uniquenames )
					if ( isTRUE( all.equal( length( highlightNodes2 ) , 0 ) ) )
					{
						next
					}
					xcolor[highlightNodes2] <- highColor[il]
					xshape[highlightNodes2] <- highShape[il]
					xNsize[highlightNodes2] <- highNodeSize[il]
					xSsize[highlightNodes2] <- highFontSize[il]
				}
			}
			verticesMatrix[,2] <- as.character( xcolor )
			verticesMatrix[,3] <- as.character( xshape )
			verticesMatrix[,4] <- as.character( xNsize )
			verticesMatrix[,5] <- as.character( xSsize )
		}
		#verticesMatrix <- as.matrix(verticesMatrix)
		colnames( verticesMatrix ) <- c( "nodename" , "color" , "shape" , "size" , "font_size" )
		write.table( verticesMatrix , fcysn , sep = "\t" , quote = FALSE , col.names = TRUE , row.names = FALSE )
		#**************************************************************************
		#
		# 3. output indexed netpairs
		#
		# Legend
		if ( !any( is.na( legendtable ) ) )
		{
			mhead <- paste( "Legend" , dim( legendtable )[1] , sep = " " )
			write.table( as.matrix( mhead ) , snafile , sep = "\t" , quote = FALSE , col.names = FALSE , row.names = FALSE , append = FALSE )
			cat( legendtable , file = snafile , sep = "\t" , labels = colnames( legendtable ) , append= TRUE )
			# vertex
			write.table( as.matrix( header ) , snafile , sep = "\t" , quote = FALSE , col.names = FALSE , row.names = FALSE , append = TRUE )
		}
		else
		{
			# vertex
			write.table( as.matrix( header ) , snafile , sep = "\t" , quote = FALSE , col.names = FALSE , row.names = FALSE , append = FALSE )
		}
		cat( verticesMatrix , file = snafile , sep = "\t" , labels = colnames( verticesMatrix ) , append = TRUE )
		#edge color
		write.table( t( as.matrix( c( "edge_color_levels" , edgecolorlevels ) ) ) , snafile , sep = "\t" , quote = FALSE , col.names = FALSE , row.names = FALSE , append = TRUE )
		# link pairs based index   
		#
		write.table( rbind( c( "src" , "dst" , "type" ) ) , snafile , sep = "\t" , quote = FALSE , col.names = FALSE , row.names = FALSE , append = TRUE )
		write.table( linksIdxMatrix , snafile , sep = "\t" , quote = FALSE , col.names = FALSE , row.names = FALSE , append = TRUE )
		return( c( fcys , fcysn ) )
	}
	
	splitString <- function( mystring , separator = "; " )
	{
	  splitted <- NULL
	  for ( each in mystring )
	  {
		 if ( is.na( each ) | is.null( each ) )
		 {
			next
		 }
		 splitted <- c(splitted , unlist( strsplit( each , separator ) ) )
	  }
	  #a=unlist( strsplit(mystring, separator) )
	  return( splitted )
	}

	
    keydriver.plot <- function(keydriver.fname, TOP.NUM, edge.raw.fname, node.raw.fname, out.kda.FOLD)
	{
	  #### keydriver.plot generate the edge and node file to plot key drivers
	  #### Author: Jialiang Yang  
	  #### Date: 3/19/2015
	  ##################################################################
	  ### Input: 
	  ###        keydriver.fname: the key driver file provided by KDA
	  ###        TOP.NUM: the top key drivers to plot
	  ###        edge.raw.fname: the raw edge file
	  ###        node.raw.fname: the raw node file specify properties of all nodes
	  ###        out.kda.FOLD: the output file name
	  ###*****************************************************************
	  
		  driver.matrix <- read.delim(keydriver.fname, header=TRUE, stringsAsFactors=FALSE)
		  driver.gene <- driver.matrix[, "keydrivers"]
		  driver.log <- -log10(driver.matrix[, "pvalue_corrected_subnet"])
		  
		  Top.driver.gene <- driver.gene[1:TOP.NUM]
		  Top.driver.log <- driver.log[1:TOP.NUM]
		  
		  ### deal with edge
		  ### ------------------------------------------------------------------------------
		      edge.matrix <-  read.delim(edge.raw.fname, header=TRUE, stringsAsFactors=FALSE)
		      print(edge.matrix[1:5,])
		     
			  ### find genes connected to the key drivers
			  edge.keep.judge <- apply(edge.matrix, 1, function(current.edge)
			  {
				if((current.edge[1] %in% Top.driver.gene) || (current.edge[2] %in% Top.driver.gene))
				{return(TRUE)}else{return(FALSE)}
			  }) 
		      
			  edge.keep.matrix <- edge.matrix[edge.keep.judge,]
			  colnames(edge.keep.matrix) <-c("protein1", "protein2")
			  node.keep.vec <- unique(c(edge.keep.matrix[,1], edge.keep.matrix[,2]))
			  
		  ### deal with node
		  ### ------------------------------------------------------------------------------
			  node.matrix <- read.delim(node.raw.fname, header=TRUE, stringsAsFactors=FALSE)
			  print(node.matrix[1:5,])
			  
			  node.ind <- match(node.keep.vec, as.character(node.matrix[, "node"]))
			  node.select.matrix <- node.matrix[node.ind,]
			  
			  ### decide node color
			  color.vec <- as.character(unlist(apply(node.select.matrix, 1, function(current.node)
			  {
			    if((current.node[5] == 1) && (current.node[6] == 1)){return("light green")}
				else if((current.node[5] == 1) && (current.node[6] == 0)){return("red")}
				else if((current.node[5] == 0) && (current.node[6] == 1)){return("blue")}
				else{return("grey")}
			  })))
			  
			  ### remove grey nodes
              node.keep.ind <- union(which(color.vec != "grey"), match(Top.driver.gene, as.character(node.select.matrix[, "node"])))	
			  
              node.keep.matrix <- node.select.matrix[node.keep.ind,]	
              node.final <- node.keep.matrix[, "node"]
              node.color.final <- color.vec[node.keep.ind]			  

			  ### node shape
			  node.keydriver.ind <- match(Top.driver.gene, node.final)
			  driver.gene <- rep(0, length(node.final))
			  shape.vec <- rep("circle", length(node.final))
			  driver.gene[node.keydriver.ind] <- 1
			  shape.vec[node.keydriver.ind] <- "square"
			  
			  
			  ### node size
			  size.vec <- rep(50, length(node.final))
			  size.vec[which(shape.vec=="square")] <- 100
			  
			  font.vec <- size.vec - 10
			  
			  node.new.matrix <- cbind(node.final, shape.vec, node.color.final, size.vec, font.vec, driver.gene, node.keep.matrix[, 2:9])
			  colnames(node.new.matrix)[1:6] <- c("node", "shape", "color", "size", "font", "driver")
	   
	          write.table(node.new.matrix,file=paste0(out.kda.FOLD, "KDA.node.agingdisease.txt"),sep='\t',quote=F,row.names=F)
			  
			### keep only edges not involve grey nodes
			### ------------------------------------------------------------------
                edge.grey.judge <- apply(edge.keep.matrix, 1, function(current.edge)
				{
					if((current.edge[1] %in% node.final) && (current.edge[2] %in% node.final))
					{return(TRUE)}else{return(FALSE)}
			    }) 

			    edge.matrix.final <- edge.keep.matrix[edge.grey.judge,]
			  
			    write.table(edge.matrix.final,file=paste0(out.kda.FOLD, "KDA.edge.agingdisease.txt"),sep='\t',quote=F,row.names=F)
	}
	
	 