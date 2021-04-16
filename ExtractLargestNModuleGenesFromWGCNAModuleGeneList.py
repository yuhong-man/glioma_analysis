#!/data2/wangb/bin/python3
import sys
import pandas as pd  

if __name__ == "__main__":
    '''
    Input is like: 
    name    label   color
    ENSG00000000003.13      52      salmon4
    ENSG00000000005.5       32      violet
    ENSG00000000419.11      19      lightyellow
    ENSG00000000457.12      1       turquoise
    ENSG00000000460.15      1       turquoise
    ENSG00000000938.11      13      salmon
    ENSG00000000971.14      20      royalblue
    ENSG00000001036.12      13      salmon
    ENSG00000001084.9       6       red
    ENSG00000001167.13      1       turquoise
    ENSG00000001460.16      1       turquoise
    Output is a gene list (each per line). 
    '''
    path_module_gene_list = sys.argv[1]
    largest_n = int(sys.argv[2])
    at_least_gene_num = int(sys.argv[3])
    path_output_gene_list = sys.argv[4]

    df_module_gene_list = pd.read_csv(path_module_gene_list, sep='\t')

    modules = df_module_gene_list['color'].unique()

    df_color_counts = df_module_gene_list['color'].value_counts() 
    selected_colors = list(df_color_counts[df_color_counts >= at_least_gene_num][:largest_n].index)
    df_module_gene_list_selected = df_module_gene_list[df_module_gene_list['color'].isin(selected_colors)] 
    df_module_gene_list_selected['name'].to_csv(path_output_gene_list, sep='\t', index=False, header=None)