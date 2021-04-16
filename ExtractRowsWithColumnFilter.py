#!/data2/wangb/bin/python3
import sys 
import pandas as pd 

def extractRowsWithColumn_LargerThan_(path_in, column, filter_condition):
    df_in = pd.read_csv(path_in, sep='\t')
    df_out = df_in[eval("df_in[column]{}".format(filter_condition))]
    return df_out 

if __name__ == "__main__":
    path, column, filter_condition, path_out = sys.argv[1:5] 
    df_out = extractRowsWithColumn_LargerThan_(path, column, filter_condition)
    df_out.to_csv(path_out, sep='\t', index=False)