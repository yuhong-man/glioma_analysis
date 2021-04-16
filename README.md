# Code, data and usage for DLBC analysis including DEG, WGCNA etc. 
## Code
Files with suffix in R, pl, sh, py.
## Data
txt.
## Usage
Here is a all-in-one code for the whole pipeline. 
DESeq2WGCNAPipeline.sh \
$path_input_read_counts
$path_group_data \
$path_clinical_traits \
$folder_for_results \
$prefix \
$use_deg_for_wgcna \
$use_wgcna_large_module_genes_for_KDA \
$THRESHOLD_FOR_DEG_LOG2FOLDCHANGE \
$THRESHOLD_FOR_DEG_PADJ \
$INPUT_GENE_TYPE \
$DO_DGCA

### Explainations for some input:  
path_input_read_counts: the read count matrix. 
path_group_data: group information.
path_clinical_traits: clinical trait information. 
use_deg_for_wgcna: whether to use only DEGs for WGCNA analysis. 
use_wgcna_large_module_genes_for_KDA: whether to use genes from large modules form WGCNA for key driver analysis. 
THRESHOLD_FOR_DEG_LOG2FOLDCHANGE: for DESeq2. 
THRESHOLD_FOR_DEG_PADJ: for DESeq2. 
INPUT_GENE_TYPE: SYMBOL or ENSEMBL
DO_DGCA: whether run DGCA or not. DGCA is too computationally expensive. 

## Example
ipynb
