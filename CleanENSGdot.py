#!/data2/wangb/bin/python3
description = \
'''
For removing the suffix .xxx in ENSGxxxxx.xxx. 
'''
import sys 
import pandas as pd 

if __name__ == "__main__":
    path_in, path_out = sys.argv[1:3]
    df_in = pd.read_csv(path_in, sep='\t', index_col=0) 
    index_name = df_in.index.name
    index_name = str(index_name).split('.')[0]
    df_in.index = [str(i).split('.')[0] for i in df_in.index]
    df_in.index.name = index_name
    print(df_in.head())
    df_in.to_csv(path_out, sep='\t')
