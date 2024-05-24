#!/bin/bash

#SBATCH --mem=60G
#SBATCH --cpus-per-task=20


#########################################################
### This script:
# creates a gene expression count table using Cellranger
# input: 10x scRNAseq fastq files
#########################################################


#####################################################################################################################################################
################################################################ USER-DEFINED VARIABLES #############################################################


## if you want files to be md5sum checked, specify containing directories
# copy md5sum files into directories that should be md5sum checked
md5sum_dirs=(
)


####################################################################################################################################################



module load cellranger/7.2.0



## check MD5sums
# loop through list of md5sum directories
len=${#md5sum_dirs[@]}

for (( k=0; k<$len; k++ ))
do
	cd ${md5sum_dirs[k]}

	# check MD5 sums are correct
	sum_file=$(find . -maxdepth 1 -regex '.*md5sum\.txt')
	md5sum -c $sum_file

	wait
done



## run the cellranger count command
cd $cellranger_dir


# choose the correct transcriptome
if [ $species = "human" ] ; then
  	transcriptome_ref="${genome_dir}/refdata-gex-GRCh38-3.0.0"
elif [ $species = "mouse" ] ; then
  	transcriptome_ref="${genome_dir}/refdata-gex-mm10-2020-A"
elif [ $species = "zebrafish" ] ; then
  	transcriptome_ref="${genome_dir}/refdata-cellranger-${genome}"
fi


for j in $(sed '1d' ${sample_info} | awk -F\\t '{print $1}')
do
	# extract sample id
	sample=$(grep $j $sample_info | awk -F\\t '{print $1}' | tr , _)

	if [ ! -f "${cellranger_dir}/${sample}/outs/molecule_info.h5" ]; then
		# extract additional columns
		categories=$(grep $j $sample_info | awk -F\\t '{for(i=3;i<=NF;i++) printf $i","; print ""}')

		# populate sample data file
		echo ${sample},${cellranger_dir}/${sample}/outs/molecule_info.h5,${categories} >> cellranger_aggr_input.csv
	
		# run cellranger count command
		echo Cellranger is counting ${j}

		fastq_dir=$(grep $j $sample_info | awk -F\\t '{print $2}')	

		cellranger count --id=${sample} \
		--fastqs=${fastq_dir} \
		--transcriptome=${transcriptome_ref} \
		--sample=${j}
	fi
done



## create sample data .csv file (required for aggr command)
# create header for aggregator input file
sed 1q ${sample_info} | awk -F\\t -vOFS=, '{$1=$2=""; print $0}' | sed 's/^,//' > ${cellranger_dir}/cellranger_aggr_input.csv
sed -i "s/^/sample_id,molecule_h5/" ${cellranger_dir}/cellranger_aggr_input.csv

# create sample entries
for j in $(sed '1d' ${sample_info} | awk -F\\t '{print $1}')
do
	# extract sample id
	sample=$(grep ${j} ${sample_info} | awk -F\\t '{print $1}' | tr , _)

	# create molecule_info h5 path name
	h5_path=${cellranger_dir}/${sample}/outs/molecule_info.h5

	# extract additional column entries and paste together with comma as separator
	categories=$(grep ${j} ${sample_info} | awk -F\\t -vOFS=, '{$1=$2=""; print $0}' | sed 's/^,//')

	# append to .csv file
	echo ${sample},${h5_path}${categories} >> ${cellranger_dir}/cellranger_aggr_input.csv
done



## run cellranger aggr command to combine all sample counts
echo Cellranger is aggregating counts of all samples

cellranger aggr --id=${prefix_aggr} \
--csv=${cellranger_dir}/cellranger_aggr_input.csv

wait



## generate a tsv file containing cell barcode metadata
# create header for cell metadata file
sed 1q ${cellranger_dir}/cellranger_aggr_input.csv | awk -F, '{for(i=3;i<=NF;i++) printf $i"\t"; print ""}' > ${cellranger_dir}/cellranger_aggr_cell_metadata.tsv
sed -i "s/^/barcode\tsample_id\t/" ${cellranger_dir}/cellranger_aggr_cell_metadata.tsv

# populate metadata file
sample_count=1

for j in $(sed '1d' ${cellranger_dir}/cellranger_aggr_input.csv | awk -F, '{print $1}')
do
	# extract cell barcodes
	gunzip ${j}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
	cat ${j}/outs/filtered_feature_bc_matrix/barcodes.tsv | awk -F\\t '{print $(NF-0)}' > ${cellranger_dir}/cellranger_aggr_cell_metadata_${j}.tsv

	# rename barcodes by replacing last letter "1" with the sample count (similar to what cellranger aggr does)
	sed "s/\(.\)$/${sample_count}/" ${cellranger_dir}/cellranger_aggr_cell_metadata_${j}.tsv > ${cellranger_dir}/cellranger_aggr_cell_metadata_${j}_2.tsv

	# add sample ID
	sample_id=$(grep $j ${cellranger_dir}/cellranger_aggr_input.csv | awk -F, '{print $1}')
	sed -i "s/$/\t$sample_id/" ${cellranger_dir}/cellranger_aggr_cell_metadata_${j}_2.tsv

	# add sample metadata
	# get metadata from aggregator input sheet, tab-delimited, trailing tab deleted
	categories=$(grep $j ${cellranger_dir}/cellranger_aggr_input.csv | awk -F, '{for(i=3;i<=NF;i++) printf $i"\t"; print ""}' | sed 's/.$//')
	sed -i "s/$/\t$categories/" ${cellranger_dir}/cellranger_aggr_cell_metadata_${j}_2.tsv

	# append everything to combined output file
	cat ${cellranger_dir}/cellranger_aggr_cell_metadata_${j}_2.tsv >> ${cellranger_dir}/cellranger_aggr_cell_metadata.tsv

	rm ${cellranger_dir}/cellranger_aggr_cell_metadata_${j}.tsv
	rm ${cellranger_dir}/cellranger_aggr_cell_metadata_${j}_2.tsv
	gzip ${j}/outs/filtered_feature_bc_matrix/barcodes.tsv

	sample_count=`expr $sample_count + 1`
done



echo All done!