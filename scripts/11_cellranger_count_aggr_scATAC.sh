#!/bin/bash

#SBATCH --mem=60G
#SBATCH --cpus-per-task=20


#########################################################
### This script:
# creates a gene expression count table using Cellranger
# input: 10x scATAC-seq fastq files
#########################################################


#####################################################################################################################################################
################################################################ USER-DEFINED VARIABLES #############################################################


## if you want files to be md5sum checked, specify containing directories
# copy md5sum files into directories that should be md5sum checked
md5sum_dirs_ATAC=(
)


####################################################################################################################################################



module load cellranger-atac



## check MD5sums
# loop through list of md5sum directories
len=${#md5sum_dirs_ATAC[@]}

for (( k=0; k<$len; k++ ))
do
	cd ${md5sum_dirs_ATAC[k]}

	# check MD5 sums are correct
	sum_file=$(find . -maxdepth 1 -regex '.*md5sum\.txt')
	md5sum -c $sum_file

	wait
done



## run the cellranger count command
cellranger_dir_ATAC=${out_dir}/cellranger_ATAC
[ ! -d "$cellranger_dir_ATAC" ] && mkdir -p "$cellranger_dir_ATAC"

cd $cellranger_dir_ATAC


# choose the correct genome
if [ $species = "human" ] ; then
  	genome_ref="${genome_dir}/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
elif [ $species = "mouse" ] ; then
  	genome_ref="${genome_dir}/refdata-cellranger-arc-mm10-2020-A-2.0.0"
elif [ $species = "zebrafish" ] ; then
  	genome_ref="${genome_dir}/refdata-cellranger-atac-${genome}"
fi


for j in $(sed '1d' ${sample_info_ATAC} | awk -F\\t '{print $1}')
do
	# extract sample id
	sample=$(grep $j $sample_info_ATAC | awk -F\\t '{print $1}' | tr , _)

	if [ ! -f "${cellranger_dir_ATAC}/${sample}/outs/fragments.tsv.gz" ]; then

		# run cellranger count command
		echo Cellranger is counting ${j}

		fastq_dir=$(grep $j $sample_info_ATAC | awk -F\\t '{print $2}')	

		cellranger-atac count --id=${sample} \
		--fastqs=${fastq_dir} \
		--reference=${genome_ref} \
		--sample=${j}
	fi
done



## create sample data .csv file (required for aggr command)
# create header for aggregator input file
sed 1q ${sample_info_ATAC} | awk -F\\t -vOFS=, '{$1=$2=""; print $0}' | sed 's/^,//' > ${cellranger_dir_ATAC}/cellranger_aggr_input.csv
sed -i "s/^/library_id,fragments,cells/" ${cellranger_dir_ATAC}/cellranger_aggr_input.csv

# create sample entries
for j in $(sed '1d' ${sample_info_ATAC} | awk -F\\t '{print $1}')
do
	# extract sample id
	sample=$(grep ${j} ${sample_info_ATAC} | awk -F\\t '{print $1}' | tr , _)

	# create fragments.tsv.gz path
	fragments_path=${cellranger_dir_ATAC}/${sample}/outs/fragments.tsv.gz

	# create singlecell.csv path
	singlecell_path=${cellranger_dir_ATAC}/${sample}/outs/singlecell.csv

	# extract additional column entries and paste together with comma as separator
	categories=$(grep ${j} ${sample_info_ATAC} | awk -F\\t -vOFS=, '{$1=$2=""; print $0}' | sed 's/^,//')

	# append to .csv file
	echo ${sample},${fragments_path},${singlecell_path}${categories} >> ${cellranger_dir_ATAC}/cellranger_aggr_input.csv
done



## run cellranger aggr command to combine all sample counts
echo Cellranger is aggregating counts of all samples

cellranger-atac aggr --id=${prefix_aggr} \
	--csv=${cellranger_dir_ATAC}/cellranger_aggr_input.csv \
	--reference=${genome_ref}

wait



## generate a tsv file containing cell barcode metadata
# create header for cell metadata file
sed 1q ${cellranger_dir_ATAC}/cellranger_aggr_input.csv | awk -F, '{for(i=4;i<=NF;i++) printf $i"\t"; print ""}' > ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata.tsv
sed -i "s/^/barcode\tsample_id\t/" ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata.tsv

# populate metadata file
sample_count=1

for j in $(sed '1d' ${cellranger_dir_ATAC}/cellranger_aggr_input.csv | awk -F, '{print $1}')
do
	# extract cell barcodes
	sed '1d' ${j}/outs/singlecell.csv | awk -F, '{print $1}' > ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata_${j}.tsv

	# rename barcodes by replacing last letter "1" with the sample count (similar to what cellranger aggr does)
	sed "s/\(.\)$/${sample_count}/" ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata_${j}.tsv > ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata_${j}_2.tsv

	# add sample ID
	sample_id=$(grep $j ${cellranger_dir_ATAC}/cellranger_aggr_input.csv | awk -F, '{print $1}')
	sed -i "s/$/\t$sample_id/" ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata_${j}_2.tsv

	# add sample metadata
	# get metadata from aggregator input sheet, tab-delimited, trailing tab deleted
	categories=$(grep $j ${cellranger_dir_ATAC}/cellranger_aggr_input.csv | awk -F, '{for(i=4;i<=NF;i++) printf $i"\t"; print ""}' | sed 's/.$//')
	sed -i "s/$/\t$categories/" ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata_${j}_2.tsv

	# append everything to combined output file
	cat ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata_${j}_2.tsv >> ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata.tsv

	rm ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata_${j}.tsv
	rm ${cellranger_dir_ATAC}/cellranger_aggr_cell_metadata_${j}_2.tsv

	sample_count=`expr $sample_count + 1`
done



echo All done!