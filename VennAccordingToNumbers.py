#!/data2/wangb/bin/python3
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd
import argparse

if __name__ == "__main__": 
    parser = argparse.ArgumentParser()  
    parser.add_argument('-f', dest='path', help="input file path", type=str, nargs=1) 
    parser.add_argument('-o', dest='output_prefix', help="output_prefix", type=str, nargs=1)
    parser.add_argument('-t', dest='title_column', help="title_column", type=str, nargs=1)
    parser.add_argument('-p', dest='p_column', help="p value column", type=str, nargs=1)
    parser.add_argument('-P', dest='pvalue_threshold', help="p value column", type=float, nargs=1, default=0.01)
    parser.add_argument('-A', dest='A', help="column for set A number", type=str, nargs=1)
    parser.add_argument('-a', dest='a', help="column for set A name", type=str, nargs=1)
    parser.add_argument('-B', dest='B', help="column for set B number", type=str, nargs=1)
    parser.add_argument('-b', dest='b', help="column for set B name", type=str, nargs=1)
    parser.add_argument('-i', dest='intersect', help="column for intersected number", type=str, nargs=1)
    args = parser.parse_args()

    file_in = args.path[0]
    df_in = pd.read_csv(file_in, sep='\t')
    # print(df_in.columns)
    new_columns = ['index'] + list(df_in.columns[1:])
    # print(new_columns)
    df_in.columns = new_columns
    title_column = args.title_column[0]
    # if title_column in df_in.columns:
    #     pass 
    # else:
    #     try:
    #         title_column = int(title_column)
    #         title_column = df_in.columns[title_column]
    #     except:
    #         raise Exception('title_column {} not in table'.format(title_column))

    output_prefix = args.output_prefix[0] 

    for index,row in df_in.iterrows():
        pvalue = row[args.p_column[0]]
        if pvalue > args.pvalue_threshold:
            continue
        intersect_num = row[args.intersect[0]]
        A_num = row[args.A[0]]
        A_minus_B_num = A_num - intersect_num
        B_num = row[args.B[0]]
        B_minus_A_num = B_num - intersect_num
        plt.figure(pvalue)
        title_a = row[args.a[0]] if args.a[0] in row.index else args.a[0]
        title_b = row[args.b[0]] if args.b[0] in row.index else args.b[0]
        venn2(subsets=(A_minus_B_num, B_minus_A_num, intersect_num), set_labels=(title_a, title_b))
        plt.title(pvalue)
        # plt.tight_layout()
        path_output = '{}_{}.pdf'.format(output_prefix, row[title_column]).replace(' ','_').replace(':','_').replace('|','_').replace('/','_')[-254:]
        plt.savefig(path_output)
        plt.close()
        
