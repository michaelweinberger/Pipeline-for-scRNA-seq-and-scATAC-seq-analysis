#!/bin/bash

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=3-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15
#SBATCH --partition=long



#############################################################
####                                                    #####
####              User defined variables                #####
####                                                    #####
#############################################################

#################################################################################################################################################


# specify project name
project="Farbehi"

# specify directory containing parent and child scripts
script_dir="/ceph/project/tsslab/mweinber/tmp/scripts_single_cell"

# specify output directory
# if directory does not exist, it will be generated
out_dir="/ceph/project/tsslab/mweinber/tmp/out_single_cell"

# specify species ("human", "mouse", "zebrafish")
species="mouse"






##############   scRNA-seq Mapping   ##############

# indicate if scRNA-seq mapping via Cellranger should be run ("Yes" or "No")
scRNA_mapping="Yes"

# for cellranger count and aggr, supply a sample information txt file
# should be tab-delimited, contain a header and the first 2 columns are required to be:
# "sample_id" (the sample-specific prefix of the fastq files, if samples have been re-sequenced supply all prefixes as a comma-separated list)
# "fastq_dir" (the directory containing the fastq files specified by "sample_id", if samples have been re-sequenced and are stored in multiple directories supply all directories as a comma-separated list)
# Further optional columns can be added: e.g. "sample_name", "tissue", "condition" etc.
sample_info="/project/tsslab/mweinber/2023_datasets/datasets_scRNAseq/Farbehi_et_al_10.7554eLife.43882/Farbehi_et_al_sample_setup.txt"



##############   scRNA-seq Doublet removal + Clustering   ##############

# indicate if scRNA-seq clustering should be run ("Yes" or "No")
scRNA_clustering="Yes"

# indicate if Scanpy or Seurat should be used for clustering ("scanpy" or "seurat")
scRNA_analysis="scanpy"

# supply path to a cellranger output directory containing "barcodes.tsv.gz", "features.tsv.gz" and "matrix.mtx.gz" (normally path/to/outs/count/filtered_feature_bc_matrix).
# Does not need to be supplied if scRNA-seq Mapping (above) has been run, in which case the generated scRNA-seq Mapping outputs will be used.
cellranger_out_dir=""

# supply path to a metadata .tsv file containing a "barcode" column of cell barcodes, a "sample_id" column of sample identifiers, and optional metadata columns.
# Does not need to be supplied if scRNA-seq mapping (above) has been run, in which case the generated scRNA-seq mapping outputs will be used (${out_dir}/cellranger/cellranger_aggr_cell_metadata.tsv).
metadata=""

# indicate the Leiden clustering resolution to be used
leiden_res=0.4



##############   scRNA-seq Annotation + Visualization   ##############

# indicate if scRNA-seq annotation should be run ("Yes" or "No")
scRNA_annotation="No"

# supply a .csv file containing a "cluster" column with cell cluster numbers and a "cell_type" column with cell type labels,
# optionally include an "order" column with integers indicating the order in which cell types should appear in plot legends
cluster_anno="/ceph/project/tsslab/mweinber/tmp/out_single_cell/scanpy/Farbehi_annotation.csv"

# supply path to a clustered scanpy (file ending ".h5ad") or seurat (file ending ".rds") object,
# scanpy object metadata need to contain "leiden" column with cell cluster labels, seurat object metadata need to contain "seurat_clusters" column.  
# Does not need to be supplied if scRNA-seq annotation (above) has been run, in which case the generated output object will be used.
scRNA_annotation_input=""



##############   scRNA-seq velocity analysis   ##############

# indicate if scRNA-seq clustering should be run ("Yes" or "No")
scRNA_velocity="No"

# supply path to a directory containing cellranger "outs/possorted_genome_bam.bam" scRNA-seq mapping BAM files.
# Does not need to be supplied if scRNA-seq mapping (above) has been run, in which case the generated mapping outputs will be used.
scRNA_velocity_cellranger_dir=""

# supply path to a clustered scanpy (file ending ".h5ad") or seurat (file ending ".rds") object, cell barcodes should match those in "scRNA_velocity_cellranger_dir" BAM files,
# object metadata need to contain "barcode" column with cell barcodes and "cell_type" column with cell type labels.  
# Does not need to be supplied if scRNA-seq annotation (above) has been run, in which case the generated output object will be used.
scvelo_input=""






##############   scATAC-seq Mapping   ##############

# indicate if scATAC-seq mapping via Cellranger-ATAC should be run ("Yes" or "No")
scATAC_mapping="No"

# for cellranger-atac count and aggr, supply a sample information txt file
# should be tab-delimited, contain a header and the first 2 columns are required to be:
# "sample_id" (the sample-specific prefix of the fastq files, if samples have been re-sequenced supply all prefixes as a comma-separated list)
# "fastq_dir" (the directory containing the fastq files specified by "sample_id", if samples have been re-sequenced and are stored in multiple directories supply all directories as a comma-separated list)
# Further optional columns can be added: e.g. "sample_name", "tissue", "condition" etc.
sample_info_ATAC=""



##############   scATAC-seq Clustering + Annotation  ##############

# indicate if scATAC-seq clustering should be run ("Yes" or "No")
scATAC_clustering="No"

# supply path to a cellranger-ATAC output "fragments.tsv.gz" file.
# Does not need to be supplied if scATAC-seq Mapping (above) has been run, in which case the generated mapping output will be used.
cellranger_out_ATAC="/ceph/project/tsslab/mweinber/2023_datasets/datasets_scATACseq/Wang_et_al_10.1016j.celrep.2020.108472/GSE153479_fragments.tsv.gz"

# supply path to a metadata .tsv file containing cell barcodes as row names, a "sample_id" column of sample identifiers (for sample integration), and optional metadata columns.
# Does not need to be supplied if scATAC-seq mapping (above) has been run, in which case the generated mapping outputs will be used.
metadata_ATAC="/ceph/project/tsslab/mweinber/2023_datasets/datasets_scATACseq/Wang_et_al_10.1016j.celrep.2020.108472/signac/Wang_metadata_ATAC_barcodes.tsv"

# indicate the Leiden clustering resolution to be used
leiden_res_ATAC=0.4

# (optional) supply filepath to Seurat (file ending ".rds") or Scanpy (file ending ".h5ad") scRNA-seq object to be used for cell type annotation of scATAC-seq data, metadata needs to contain a column "cell_type" with cell type annotation.
# If you would like to use the output of the scRNA-seq analysis above for cell type label transfer,
# set scRNA_path to (if scRNA_analysis parameter is "scanpy"): ${out_dir}/scanpy/{project}_scRNAseq_no_doublets_annotated.h5ad 
# (no quotation marks)
# or set scRNA_path to (if scRNA_analysis parameter is "seurat"): ${out_dir}/seurat/{project}_scRNAseq_no_doublets_annotated.rds 
# (no quotation marks)
scRNA_path="/ceph/project/tsslab/mweinber/2023_datasets/datasets_scRNAseq/Wang_et_al_10.1016j.celrep.2020.108472/Scanpy/Wang_analysed_final.h5ad"


#################################################################################################################################################





### Analysis

## load modules
module load R-cbrg
module load python-cbrg



## create output directory if it does not exist
[ ! -d "$out_dir" ] && mkdir -p "$out_dir"



## prepare genome files
# specify additional genome information
if [ "$species" == "human" ] ; then
	species_latin="homo_sapiens"
	species_latin_2="Homo_sapiens"
	species_latin_short="hsapiens"
	genome="GRCh38"
	genome_ucsc="hg38"
	chrom_number=22 #X,Y chromosomes excluded
elif [ "$species" == "mouse" ] ; then
	species_latin="mus_musculus"
	species_latin_2="Mus_musculus"
	species_latin_short="mmusculus"
	genome="GRCm39"
	genome_ucsc="mm39"
	chrom_number=19 #X,Y chromosomes excluded
elif [ "$species" == "zebrafish" ] ; then
	species_latin="danio_rerio"
	species_latin_2="Danio_rerio"
	species_latin_short="drerio"
	genome="GRCz11"
	genome_ucsc="danRer11"
	chrom_number=25
fi



## download genome files
# check if directory for genome files exists, if not create it using mkdir
genome_dir=${out_dir}/genomes
[ ! -d "$genome_dir" ] && mkdir -p "$genome_dir"

source ${script_dir}/2_genome_files_SH.sh



## run cellranger count + aggr
cellranger_dir=${out_dir}/cellranger
[ ! -d "$cellranger_dir" ] && mkdir -p "$cellranger_dir"

if [ "$scRNA_mapping" == "Yes" ] ; then

	prefix_aggr=$project
	source ${script_dir}/3_cellranger_count_aggr.sh
fi



## construct Anndata or Seurat object, remove doublets and cluster
if [ "$scRNA_clustering" = "Yes" ] ; then

	if [ ! -d "$cellranger_out_dir" ] ; then
		cellranger_out_dir=${out_dir}/cellranger/${project}/outs/count/filtered_feature_bc_matrix
	fi

	if [ ! -f "$metadata" ] ; then	
		metadata=${out_dir}/cellranger/cellranger_aggr_cell_metadata.tsv
	fi

	if [ "$scRNA_analysis" == "scanpy" ] ; then
		# run Scanpy clustering 
		python ${script_dir}/4_scanpy.py -i $cellranger_out_dir -m $metadata -o ${out_dir}/scanpy -n $project -r $leiden_res
	elif [ "$scRNA_analysis" == "seurat" ] ; then
		# run Seurat clustering
		R -f ${script_dir}/4_seurat.R --args data_dir=$cellranger_out_dir metadata_file=$metadata out_dir=${out_dir}/seurat project_name=$project leiden_res=$leiden_res
	fi
fi



## add cell type annotation and generate UMAP plots
if [ "$scRNA_annotation" == "Yes" ] ; then

	# define input scRNA-seq object
	if [ ! -f "$scRNA_annotation_input" ] ; then
		if [ "$scRNA_analysis" == "scanpy" ] ; then 
			scRNA_annotation_input=${out_dir}/scanpy/${project}_scRNAseq_analysed_no_doublets.h5ad
		elif [ "$scRNA_analysis" == "seurat" ] ; then
			scRNA_annotation_input=${out_dir}/seurat/${project}_scRNAseq_analysed_no_doublet.rds
		fi
	fi

	# run annotation script
	if [ "${scRNA_annotation_input##*.}" == "h5ad" ] ; then 
		python ${script_dir}/5_scanpy_plotting.py -i $scRNA_annotation_input -a $cluster_anno -o ${out_dir}/scanpy -n $project -r $leiden_res
	elif [ "${scRNA_annotation_input##*.}" == "rds" ] ; then
		R -f ${script_dir}/5_seurat_plotting.R --args in_file=$scRNA_annotation_input annotation=$cluster_anno out_dir=${out_dir}/seurat project_name=$project leiden_res=$leiden_res
	fi
fi



## perform mRNA velocity analysis
if [ "$scRNA_velocity" == "Yes" ] ; then

	# create output directory
	scvelo_dir=${out_dir}/scvelo
	if [ ! -d "$scvelo_dir" ] ; then
  		mkdir -p "$scvelo_dir"
	fi

	# define locations of input BAM files
	if [ ! -d "$scRNA_velocity_cellranger_dir" ] ; then	
		$scRNA_velocity_cellranger_dir=$cellranger_dir
	fi

	# define transcriptome data to use
	if [ $species = "human" ] ; then
		transcriptome=${genome_dir}/refdata-cellranger-GRCh38-3.0.0
	elif [ $species = "mouse" ] ; then
		transcriptome=${genome_dir}/refdata-gex-mm10-2020-A
	elif [ $species = "zebrafish" ] ; then
		transcriptome=${genome_dir}/refdata-cellranger-${genome}
	fi

	#repeats=${genome_dir}/${genome}_rmsk.txt

	# generate loom files of spliced and unspliced reads
	source ${script_dir}/6_velocity_loom_files_SH.sh


	# run python script to merge loom files
	python ${script_dir}/7_velocity_loom_merge.py -i ${scvelo_dir}/loom_files -m $cellranger_dir


	# define input scRNA-seq object for scvelo analysis
	if [ ! -f "$scvelo_input" ] ; then
		if [ "$scRNA_analysis" == "scanpy" ] ; then
			scvelo_input=${out_dir}/scanpy/${project}_analysed_final.h5ad
		elif [ "$scRNA_analysis" == "seurat" ] ; then
			scvelo_input=${out_dir}/seurat/${project}_scRNAseq_seurat_no_doublets_annotated.rds
		fi
	fi


	# if input is Seurat object, export data necessary to construct Anndata object for scvelo analysis
	if [ "${scvelo_input##*.}" == "rds" ] ; then
		R -f ${script_dir}/8_export_seurat_data_for_anndata.R --args in_file=$scvelo_input out_dir=${scvelo_dir}/seurat_data project_name=$project
	fi


	# perform scvelo analysis
	python ${script_dir}/9_scvelo.py -i $scvelo_input -o ${scvelo_dir}/${project} -mtx ${scvelo_dir}/seurat_data/velocity_counts_${project}.mtx \
		-md ${scvelo_dir}/seurat_data/velocity_metadata_${project}.csv -g ${scvelo_dir}/seurat_data/velocity_genes_${project}.csv \
		-p ${scvelo_dir}/seurat_data/velocity_pca_${project}.csv -l ${scvelo_dir}/loom_files/merged.loom -n $project
fi




## perform scATAC-seq mapping
# download genome files
source ${script_dir}/10_genome_files_ATAC_SH.sh

# perform mapping
if [ "$scATAC_mapping" == "Yes" ] ; then
	prefix_aggr=$project
	source ${script_dir}/11_cellranger_count_aggr_ATAC.sh
fi



## perform scATAC-seq clustering + annotation
if [ "$scATAC_clustering" == "Yes" ] ; then

	if [ ! -f "$cellranger_out_ATAC" ] ; then
		cellranger_out_ATAC=${out_dir}/cellranger_ATAC/${project}/outs/fragments.tsv.gz
	fi

	if [ ! -f "$metadata_ATAC" ] ; then	
		metadata_ATAC=${out_dir}/cellranger_ATAC/cellranger_aggr_cell_metadata.tsv
	fi

	R -f ${script_dir}/12_signac.R --args species=$species in_file=$cellranger_out_ATAC metadata=$metadata_ATAC \
		out_dir=${out_dir}/signac_scATAC project_name=$project scRNA_path=$scRNA_path leiden_res=$leiden_res_ATAC
fi





echo All done! 

