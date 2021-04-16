path_exp=$1
path_required_samples=$2
path_deg_list=$3
threshold_fold_change=$4
threshold_padj=$5
if  [ ! -n "$6" ] ;then
    fontsize=10
else
    fontsize=$6
fi

if  [ ! -n "$7" ] ;then
    p_value_col=padj
else
    p_value_col=$7
fi

# less $path_deg_list | awk '$3>"'"$threshold_fold_change"'"&&$7<"'"$threshold_padj"'" {print $1}' > $path_deg_list.up.genes
# less $path_deg_list | awk '$3<-"'"$threshold_fold_change"'"&&$7<"'"$threshold_padj"'" {print $1}' > $path_deg_list.down.genes
# cat  $path_deg_list.up.genes $path_deg_list.down.genes > $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.genes 

ExtractDegGenesFromDegList.py $path_deg_list $threshold_fold_change $threshold_padj $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.genes $p_value_col


ExtractMatrixRows.py $path_exp $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.genes $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.exp.txt 

ExtractMatrixColumns.py  $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.exp.txt $path_required_samples  $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.exp.required_samples.txt 
PlotHeatmapBypheatmap.R  $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.exp.required_samples.txt $fontsize $path_deg_list.$threshold_fold_change.$threshold_padj.up_down.exp.required_samples.pdf F T T F row
