#!/data2/wangb/bin/python3
import os, sys 
import pandas as pd 

BIN_PLOT_GO = 'plotGO.R'

def enrich(list_genes, path_output_prefix):
    with open(path_output_prefix, 'w') as fw:
        fw.write('\n'.join(list_genes))
    os.system('{} {}'.format(BIN_PLOT_GO, path_output_prefix))
    

def readModuleColorGeneList(path_in, prefix):
    '''
    name    label   color
    1/2-SBSRNA4     27      white
    A1BG    4       yellow
    A2LD1   9       magenta
    A2M     1       turquoise
    A2MP1   1       turquoise
    A4GALT  12      tan
    AAAS    6       red
    AACS    2       blue
    '''
    if not os.path.exists(prefix):
        try:
           os.makedirs(prefix)
        except:
            pass  
    df_in = pd.read_csv(path_in, sep='\t')
    groups = df_in.groupby('color')
    for color, df_this_color in groups:
        enrich(df_this_color.loc[:,'name'], path_output_prefix='{}/{}_{}.txt'.format(prefix, color, list(df_this_color['label'])[0]))


if __name__ == "__main__":
    import ShowHelp
    path_module_color_gene, prefix = sys.argv[1:3]
    readModuleColorGeneList(path_module_color_gene, prefix)
