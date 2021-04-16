#!/usr/bin/env Rscript

plotGoForDF <- function(df_in, path_output, input_gene_type='SYMBOL', path_network_txt='/data2/wangb/projects/20200708_PPI/protein.actions.v11.0.txt'){
	df_in = df_in[,'gene_id']
	print(head(df_in))
	# Remove the .suffix after ENSxxxxxxxxxxxx. 
	df_in = gsub('[.].{1,10}', '', df_in)
	print(head(df_in))
	temp_file = path_output
	write.table(df_in, temp_file, sep='\t', col.names=F, row.names=F, quote=F)
	system('sleep 2')
        cmd_enrich = paste('plotGO.R', temp_file, input_gene_type, '&')
	cmd_kda = paste('KeyDriverAnalysis.R', path_network_txt, temp_file, '&')
	cmd_to_run = paste(cmd_enrich, cmd_kda)
        return (cmd_to_run)
}


splitToUpAndDownDataFrameAndEnrichAndKDA <- function(path_deg_list, threshold_for_padj, threshold_for_log2fc, input_gene_type='SYMBOL'){  
	df_deg = read.table(path_deg_list, sep='\t', header=T, skipNul=T)
	print('original dim')
	print(dim(df_deg))
	# Drop rows containing NA. https://bbs.pinggu.org/thread-1170824-1-1.html
	df_deg =  df_deg[apply(df_deg, 1, function(x)!any(is.na(x))), , drop=F]
	# print(tail(df_deg))
	print('after dropping na dim')
	print(dim(df_deg))
	df_deg = df_deg[df_deg['padj'] <= threshold_for_padj, ] 
	print('after filtering by padj by ')
	print(threshold_for_padj)
	print(dim(df_deg))
	print('threshold_for_log2fc and dim of up and down ')
	print(threshold_for_log2fc)
	df_deg_up = df_deg[df_deg['log2FoldChange'] > threshold_for_log2fc, ]
	print(dim(df_deg_up))
	df_deg_down = df_deg[df_deg['log2FoldChange'] < -threshold_for_log2fc, ]
	print(dim(df_deg_down))
	all_cmds_to_run = c()	
	print(class(df_deg_up))
	for(direction_and_df_in in list(list('up', df_deg_up), list('down', df_deg_down))){
		direction = direction_and_df_in[1]
		df_in = as.data.frame(direction_and_df_in[2])
		# print(df_in)
		# print(class(df_in))
		# print(df_in[,'gene_id'])
		# current_time_txt = format(Sys.time(),format="%Y%m%d_%H%M%S.txt") 
		path_output = paste0(path_deg_list, '.', direction, '.txt')
		cmd_to_run = plotGoForDF(df_in, path_output=path_output, input_gene_type=input_gene_type) 
                all_cmds_to_run = paste(all_cmds_to_run, cmd_to_run)
	}
    all_cmds_to_run = paste(all_cmds_to_run, 'wait') 
	print(all_cmds_to_run)
	system(all_cmds_to_run)
}

if(sys.nframe() == 0){
	print('args: path_deg_list(from the deseq2 step), threshold')
	args <- commandArgs(T)
	path_deg_list = args[1]
	# The header should be: gene_id baseMean        log2FoldChange  lfcSE   stat    pvalue  padj
	threshold_for_padj = as.numeric(args[2])
	threshold_for_log2fc = as.numeric(args[3]) # Absolute value. 
	input_gene_type = args[4]
	# For both positive and negative.  
	splitToUpAndDownDataFrameAndEnrichAndKDA(path_deg_list, threshold_for_padj, threshold_for_log2fc, input_gene_type=input_gene_type)
}
