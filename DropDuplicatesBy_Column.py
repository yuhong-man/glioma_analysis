#!/data2/wangb/bin/python3
import pandas as pd 
import sys 

def drop_duplicates_by(df_in, drop_by, keep='first'):
    df_out = df_in.drop_duplicates(subset=drop_by, keep=keep)
    return df_out

if __name__ == "__main__":
    path_in, drop_by_column, path_out = sys.argv[1:4]
    try:
        keep = sys.argv[4]
    except:
        keep = 'first'
    if drop_by_column.startswith('['):
        drop_by_column = eval(drop_by_column)
    df_in = pd.read_csv(path_in, sep='\t')
    df_out = drop_duplicates_by(df_in, drop_by_column, keep=keep)
    df_out.to_csv(path_out, sep='\t', index=False)

