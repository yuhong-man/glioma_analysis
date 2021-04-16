#!/data2/wangb/bin/python3
#coding=utf-8
import sys 
import os 
import pandas as pd 
'''
ExtractMatrixColumns.py  $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.exp.txt $path_required_samples  $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.exp.required_samples.txt 
PlotHeatmapBypheatmap.R  $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.exp.required_samples.txt  10 $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.exp.required_samples.pdf F T T F row
'''

BIN_PLOT_HEATMAP = 'PlotHeatmapBypheatmap.R'

def extractTopNFromDEGList(path_deg_list, topN=10):
    df_in = pd.read_csv(path_deg_list, sep='\t')
    dfs_gt_or_lt = [] 
    for gt_or_lt in ['__gt__','__lt__']:
        df_this_direction = df_in[getattr(df_in['log2FoldChange'], gt_or_lt)(0)].sort_values('padj', ascending=True).iloc[:topN]
        dfs_gt_or_lt.append(df_this_direction)
    df_topN_bidirectional = pd.concat(dfs_gt_or_lt)
    print(df_topN_bidirectional)
    return df_topN_bidirectional

def plotDEGHeatmapFromDEGListForTopNGenes(path_exp, path_deg_list, path_sample_order, topN=10, transpose_exp=False):
    df_deg_topN = extractTopNFromDEGList(path_deg_list, topN)
    path_deg_topN = '{}.top{}.xls'.format(path_deg_list, topN)
    df_deg_topN.to_csv(path_deg_topN, sep='\t', index=False)

    df_exp = pd.read_csv(path_exp, sep='\t', index_col=0)
    if transpose_exp:df_exp = df_exp.T 
    sample_order = pd.read_csv(path_sample_order, sep='\t', index_col=0, header=None).index
    df_exp_for_plot = df_exp.loc[df_deg_topN['gene_id'],sample_order]
    path_exp_for_plot = '{}.topDEG{}.xls'.format(path_exp, topN)
    df_exp_for_plot.to_csv(path_exp_for_plot, sep='\t')

    path_output = path_exp_for_plot+'.pdf'
    cmd_plot = [BIN_PLOT_HEATMAP, path_exp_for_plot, '3', path_output, 'F T T F row']
    cmd_plot = ' '.join(cmd_plot)
    os.system(cmd_plot)

if __name__ == "__main__":    
    import ShowHelp    
    path_exp = sys.argv[1] 
    path_deg_list = sys.argv[2]   
    path_sample_order = sys.argv[3]
    plotDEGHeatmapFromDEGListForTopNGenes(path_exp, path_deg_list, path_sample_order)
