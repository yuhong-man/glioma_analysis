#!/data2/wangb/anaconda2/bin/Rscript
args <- commandArgs(TRUE)
path_exp <- args[1]
# expression matrix, with row names as genes and column names as samples. 

path_design <- args[2] 
# path_design should be like: 
# sample_id       condition
# L129_B03        N
# L281_A06        N
# L130_A10        N
# L129_F04        N
# L129_F03        T
# L226_H02        T
# L128_G01        T
# L141_A05        T

nPerms = 5

data1 <- read.table(path_exp, sep = '\t',header = TRUE, row.names = 1,stringsAsFactors=FALSE)
list_for_design <- read.table(path_design,sep = '\t',header =  TRUE, row.names = 1,stringsAsFactors=FALSE)
#data1 <- read.csv("/data2/wangb/lymphoma/merge.csv",sep = ',',header = TRUE, row.names = 1,stringsAsFactors=FALSE)
#list_for_design <- read.csv("/data2/wangb/lymphoma/list_for_design.csv",sep = ',',header =  TRUE, row.names = 1,stringsAsFactors=FALSE)
library(GO.db)
library(impute)
library(preprocessCore)
library(DGCA)

# names(list_for_design)[1] <- group
# normal <- list_for_design[,1]
# new_list <- data.frame(list_for_design,normal,stringsAsFactors=FALSE)
# new_list[which(new_list[,group]$lymphoma=="T"),1] <- 1
# new_list[which(new_list[]=="N"),1] <- 0
# new_list[which(new_list$normal=="T"),2] <- 0
# new_list[which(new_list$normal=="N"),2] <- 1
# new_list <- apply(new_list,2,as.numeric)
# quit()


column_condition = 'condition'
matrix_for_design = data.frame()
for(i in rownames(list_for_design)){
	for (c in unique(list_for_design[,column_condition])){ 
		if (list_for_design[i,column_condition] == c){			
			matrix_for_design[i,c] <- 1
		}else{
			matrix_for_design[i,c] <- 0
		}
	}
}
matrix_for_design <- as.matrix(matrix_for_design)

print(dim(data1))
print(dim(matrix_for_design))

ddcor_res = ddcorAll(inputMat = data1, design=matrix_for_design, compare = c("N", "T"), nPerms=nPerms)
# head(ddcor_res, 3)
