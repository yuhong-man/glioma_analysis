# example
# sh 02.filter_gene.sh 20160128-LUSC-RNAseqGene.txt gene_with_protein_product.txt


awk -F '\t' '{print $1}' $1 | awk -F '|' '{print $1}' | sort | uniq > A.list;
awk -F '\t' '{print $2}' $2 | sort | uniq > B.list;
comm -12 A.list B.list > filter.list;
rm A.list B.list;
