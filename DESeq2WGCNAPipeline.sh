path_input_read_counts=$1
path_group_data=$2 # path_group_data must contain T and N only. 
path_clinical_traits=$3 # could be none 
folder_for_results=$4
prefix=$5
use_deg_for_wgcna=$6
use_wgcna_large_module_genes_for_KDA=$7
THRESHOLD_FOR_DEG_LOG2FOLDCHANGE=$8 #1
THRESHOLD_FOR_DEG_PADJ=$9 #0.05
INPUT_GENE_TYPE=${10}  # SYMBOL ENSEMBL ...
DO_DGCA=${11}

# The lines about 04.step2.DESeq2_analysis.R must be adjusted manually. 
echo Usage:
head -n 15 $0

set -euo pipefail
set -v
export PATH=/data2/wangb/pipelines/deseq2/code/:$PATH



PATH_DRUG_SIGNATURE=/data2/wangb/tools/single_drug_perturbations-p1.0.gmt

# Create a dir for storing the results. 
mkdir $folder_for_results

cp $path_input_read_counts $folder_for_results/01.inputdata.txt
cp $path_group_data $folder_for_results/01.coldata.txt
if [[ -f $path_clinical_traits ]] ; then 
cp $path_clinical_traits $folder_for_results/01.clinical.txt
fi 

cd $folder_for_results

# DGCA can be run in background
if [[ $DO_DGCA. = "T". ]];then
DGCA.R 01.inputdata.txt 01.coldata.txt > DGCA.log 2>&1 &
fi
# In most cases, this will fail due to OOM error. Ignore. 

# Core function: DESeqDataSetFromMatrix and PCA
04.step1.DESeq2_analysis.R $prefix


# export PATH=/data2/wangb/pipelines/deseq2/code/:$PATH
# Core function: 
# Filter those g or l than some threshold and do PCA again. 

04.step2.DESeq2_analysis.R $prefix F  T  g 100000 g 100000 # Not filter

# The PCA below must be adjusted manully. 
# 04.step2.DESeq2_analysis.R $prefix T  T  l -4 g -100000 
# 04.step2.DESeq2_analysis.R $prefix T  N  l -4 g -100000 


# 04.step2.DESeq2_analysis_rmoutlier_ignore_group.R  $prefix T g 5  g -100000  # For coldata

04.step2.DESeq2_analysis_rmoutlier_ignore_group.R  $prefix T  l -4 g -100000 # For coldata2
 


04.step3.DESeq2_analysis.R $prefix 04.step2.DESeq2.RData 
path_deg_list=${prefix}_DESeq2_DEG.list

# Convert to symbol.  
path_deg_list_with_symbol=$path_deg_list.symbol.xls
CleanENSGdot.py $path_deg_list $path_deg_list.withoutdot

# Allows failure. 
set +euo pipefail
ConvertTCGASymbolToOrgDbSymbol.R ENSEMBL SYMBOL $path_deg_list.withoutdot $path_deg_list_with_symbol gene_id T T
DropDuplicatesBy_Column.py $path_deg_list_with_symbol gene_id $path_deg_list_with_symbol.dedup.xls
path_deg_list_with_symbol=$path_deg_list_with_symbol.dedup.xls
set -euo pipefail
if [[ -f $path_deg_list_with_symbol ]]; then 
	deg_files=($path_deg_list $path_deg_list_with_symbol)
else
	deg_files=($path_deg_list)
fi

# Plot DEG heatmap. 
sort 01.coldata.txt -k 2 -r |sed '/sample_id/d' > 01.coldata.sort.txt
PlotHeatmapByDEGList.sh 01.inputdata.txt 01.coldata.sort.txt $path_deg_list $THRESHOLD_FOR_DEG_LOG2FOLDCHANGE $THRESHOLD_FOR_DEG_PADJ

PlotHeatmapByDEGListTopNGenesRestrictReadCount.py 01.inputdata.txt $path_deg_list 01.coldata.sort.txt
# 05.step1.WGCNA_01_analysis.R $prefix N 
# 05.step1.WGCNA_01_analysis.R $prefix T 




# KeyDriverAnalysis (KDA) is included. 
for path_deg in ${deg_files[@]}
do 
	ReadDEGListAndEnrich.R $path_deg $THRESHOLD_FOR_DEG_PADJ $THRESHOLD_FOR_DEG_LOG2FOLDCHANGE $INPUT_GENE_TYPE > $path_deg.enrich.log 2>&1 & 
done 

# drugs.run.R 
for path_deg in ${deg_files[@]}
do
	output_folder_drugs_run=drugs.run_by_`basename $path_deg`
	mkdir $output_folder_drugs_run # drugs.run_after_convert_to_symbol
	drugs.run.R $path_deg $PATH_DRUG_SIGNATURE $output_folder_drugs_run/ $THRESHOLD_FOR_DEG_PADJ $THRESHOLD_FOR_DEG_LOG2FOLDCHANGE
	cd $output_folder_drugs_run/
	mkdir venn 
	cd venn 
	VennAccordingToNumbers.py  -f ../drug.up.DESeq.txt -o  drug.up.DESeq -t index -p dn_pvalue -A background_genes -a   all_genes -B drug_genes -b  index -i overlap_dn   
	VennAccordingToNumbers.py  -f ../drug.dn.DESeq.txt -o  drug.dn.DESeq -t index -p up_pvalue -A background_genes -a   all_genes -B drug_genes -b  index -i overlap_up
	cd ../..
done

# DiffCoExPipeline can be run in background. 
DiffCoExPipeline.R 01.inputdata.txt 01.coldata.txt > DiffCoExPipeline.log 2>&1 &

if [[ x${use_deg_for_wgcna} = x"T" ]]; then  # -eq is used only when they are numbers. 
	path_deg_genes=$path_deg_list.$THRESHOLD_FOR_DEG_LOG2FOLDCHANGE.$THRESHOLD_FOR_DEG_PADJ.genes
	ExtractDegGenesFromDegList.py $path_deg_list $THRESHOLD_FOR_DEG_LOG2FOLDCHANGE $THRESHOLD_FOR_DEG_PADJ $path_deg_genes
	path_genes_for_wgcna=$path_deg_genes
else
	path_genes_for_wgcna=""
fi 


echo analysing WGCNA
# Only when clinical file exists do we need to do WGCNA. 
field_of_interest=sex # or pathology_N_stage
module=blue  # 14
for group in T N
do 
	05.step1.WGCNA_01_analysis.R $prefix $group $path_genes_for_wgcna
	# 05.step2.WGCNA_01_analysis.R $prefix N 10000000 1
	# 05.step2.WGCNA_01_analysis.R $prefix T 10000000 1
	05.step2.WGCNA_01_analysis.R $prefix $group 10000000 1 
	# 05.step3.WGCNA_01_analysis_ReadClinicalTsv.R $prefix T 2#3#4#6#7#8#11#12#13
 
	# if [[ -f $folder_for_results/01.clinical.txt ]]; then
		# 05.step3.WGCNA_01_analysis_ReadClinicalTsv.R $prefix T 1#3#4#5#6#11#12#13#15#16#17#22#24#26#29#30#31
	05.step3.WGCNA_01_analysis.R $prefix $group 1#3#4#5#6#11#12#13#15#16#17#22#24#26#29#30#31
	# fi 
	
	# 05.step4.WGCNA_02_analysis.R $prefix T 
	log_file=05.step4.WGCNA_02_analysis.R.${group}.log 
	05.step4.WGCNA_02_analysis.R $prefix $group > $log_file


	POWER_THRESHOLD=`tail -n1 $log_file | cut -d ' ' -f 2 `
	05.step5.WGCNA_02_analysis.R $prefix $group $POWER_THRESHOLD

	# args[1] 为项目名称 如LUSC,BRCA
	# args[2] 为条件类型 如T,N
	05.step6.WGCNA_03_analysis.R $prefix $group
	
	05.step8.WGCNA_04_analysis.R $prefix $group $POWER_THRESHOLD & 
	# module_color=lavenderblush3
	# args[1] 为项目名称 如LUSC,BRCA
	# args[2] 为条件类型 如T,N
	# args[3] 选择感兴趣的表现数据
	# args[4] module???  
	05.step7.WGCNA_03_analysis.R $prefix $group $field_of_interest $module  

	# args[1] 为项目名称 如LUSC,BRCA
	# args[2] 为条件类型 如T,N
	# args[3] power阈值

	# args[1] 为项目名称 如LUSC,BRCA
	# args[2] 为条件类型 如T,N
	# args[3] 临床数据
	05.step9.WGCNA_04_analysis.R $prefix $group $field_of_interest
	LARGEST_N_MODULES_FOR_KDA=3
	THRESHOLD_TO_DEFINE_LARGE_MODULE=5
	path_genes_from_largest_n_modules=${group}.largest_n_modules.genes
	ExtractLargestNModuleGenesFromWGCNAModuleGeneList.py ${group}*.module_color_gene.list $LARGEST_N_MODULES_FOR_KDA $THRESHOLD_TO_DEFINE_LARGE_MODULE $path_genes_from_largest_n_modules
	ReadModuleColorGeneListAndEnrichPerColor.py  ${group}*.module_color_gene.list ${group}_EnrichmentPerColor &
	KeyDriverAnalysis.R /data2/wangb/projects/20200708_PPI/protein.actions.v11.0.txt $path_genes_from_largest_n_modules & 
	
done


wait



path_deg_list_appended_keydriver=$path_deg_list.append_keydriver.list
# drugs.run.R after appending keydrivers to DEG list. 
AppendKeyDriverToDEGList.py $path_deg_list.up.txt.KDA_keydriver.xls 0.05 pvalue_subnet $path_deg_list $path_deg_list_appended_keydriver 10
AppendKeyDriverToDEGList.py $path_deg_list.down.txt.KDA_keydriver.xls 0.05 pvalue_subnet $path_deg_list_appended_keydriver   $path_deg_list_appended_keydriver  -10

path_deg=$path_deg_list_appended_keydriver
output_folder_drugs_run=drugs.run_by_`basename $path_deg` 
mkdir $output_folder_drugs_run
drugs.run.R $path_deg $PATH_DRUG_SIGNATURE $output_folder_drugs_run/ $THRESHOLD_FOR_DEG_PADJ $THRESHOLD_FOR_DEG_LOG2FOLDCHANGE

for p in 0.05 0.01 0.001 1e-5 1e-7 1e-10 
do 
	ExtractRowsWithColumnFilter.py $path_deg_list padj "<$p" $path_deg_list.$p.xls
	for direction in up down
	do 
		ExtractRowsWithColumnFilter.py $path_deg_list.$direction.txt.KDA_keydriver.xls pvalue_subnet "<$p" $path_deg_list.$direction.txt.KDA_keydriver.xls.$p.xls
	done
	
	ExtractRowsWithColumnFilter.py $output_folder_drugs_run/drug.up.DESeq.txt dn_pvalue "<$p" $output_folder_drugs_run/drug.up.DESeq.txt.$p.txt
	ExtractRowsWithColumnFilter.py $output_folder_drugs_run/drug.dn.DESeq.txt up_pvalue "<$p" $output_folder_drugs_run/drug.dn.DESeq.txt.$p.txt
	
done 
