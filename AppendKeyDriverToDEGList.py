#!/data2/wangb/bin/python3
'''
python version of: 
path_result_from_kda=$1
cut -f 1 $path_result_from_kda > keydrivers.txt
ConvertTCGASymbolToOrgDbSymbol.R ENSEMBLPROT SYMBOL keydrivers.txt  keydrivers.txt.symbols keydrivers
cut -f 1 keydrivers.txt.symbols |  grep -v "ENSP" > symb
cat symb 

For outputing keydriver symbols filtered by pvalue.  
'''
import os 
import sys 
import pandas as pd
import tempfile

def convertGeneNameFrom_To_(from_type, to_type, path_in, which_col_to_convert, pvalue_threshold, pvalue_column, remove_initial_symbol_col='T', keep_previous_col_name='T'):
    '''
    Invoke ConvertTCGASymbolToOrgDbSymbol.R for converting gene names. 
    '''
    path_out = '{}.converted_gene_name'.format(path_in)
    path_in_filtered = '{}.filtered.xls'.format(path_in)
    df_in = pd.read_csv(path_in, sep='\t')
    print(df_in.shape)
    df_in = df_in[df_in[pvalue_column] < pvalue_threshold]
    print(df_in.shape)
    df_in.to_csv(path_in_filtered, sep='\t', index=False)
    # cmd = ['ConvertTCGASymbolToOrgDbSymbol.R', from_type, to_type, path_in_filtered, path_out, which_col_to_convert, remove_initial_symbol_col, keep_previous_col_name]
    
    cmd = ['ENSPG_toGENEID.py ', path_in_filtered, which_col_to_convert, path_out]
    
    cmd = ' '.join(map(str, cmd))
    print(cmd)
    os.system(cmd)
    df_out_converted = pd.read_csv(path_out, sep='\t')
    genes = list(df_out_converted[which_col_to_convert])
    print('genes after conversion: {}'.format(genes))
    return genes 

def appendGenesToDEG(path_deg, genes_to_append, path_deg_output, gene_id_col='gene_id', \
    fill_dict=dict(log2FoldChange=10, padj=0.001, baseMean=100, lfcSE=1, stat=1, pvalue=0.001), 
    **kargs):
    df_in = pd.read_csv(path_deg, sep='\t')
    df_append = pd.DataFrame({gene_id_col:genes_to_append})
    fill_dict.update(kargs)
    for k,v in fill_dict.items():
        df_append.loc[:,k] = v 
    df_concat = pd.concat([df_in, df_append])
    df_concat = df_concat.loc[:,df_in.columns]
    df_concat.drop_duplicates('gene_id', inplace=True, keep='first')
    df_concat.to_csv(path_deg_output, sep='\t', index=0)

    
if __name__ == "__main__":
    import ShowHelp
    path_keydriver, pvalue_threshold, pvalue_column, path_deg, path_output, log2FoldChange_to_fill = sys.argv[1:7]
    pvalue_threshold = float(pvalue_threshold)
    log2FoldChange_to_fill = float(log2FoldChange_to_fill)
    # df_in = pd.read_csv(path_in, sep='\t')
    # df_in = df_in[df_in[pvalue_column] < pvalue_threshold]
    # df_in.to_csv(sep='\t') 

    genes = convertGeneNameFrom_To_('ENSEMBLPROT', 'SYMBOL', path_keydriver, 'keydrivers', pvalue_threshold=pvalue_threshold, pvalue_column=pvalue_column)
    appendGenesToDEG(path_deg, genes, path_output, log2FoldChange=log2FoldChange_to_fill)