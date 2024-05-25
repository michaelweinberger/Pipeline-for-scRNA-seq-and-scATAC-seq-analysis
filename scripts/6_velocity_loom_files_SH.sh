#!/bin/bash

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=2-00:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=12
#SBATCH --partition=long



###############################################
### This script:
# creates spliced and unspliced read loom files
###############################################



################################################################ USER-DEFINED VARIABLES #############################################################


####################################################################################################################################################



module load python-cbrg
module load samtools/1.17



### generate output directory if not existing
if [ ! -d "$scvelo_dir" ] ; then
  mkdir -p "$scvelo_dir"
fi

cd $scvelo_dir



### loop over all subdirectories of cellranger directory
for j in $(find $scRNA_velocity_cellranger_dir -maxdepth 3 -type d)
do

	# check if directory structure is that of a cellranger count output
	if [ -f "${j}/outs/possorted_genome_bam.bam" ] ; then

		if [ ! -f "${scvelo_dir}/loom_files/${j##*/}.loom" ]; then

			echo Generating loom file for ${j}

			velocyto run10x $j \
                		 	${transcriptome}/genes/genes.gtf

			wait
		fi
	fi
done



### collect all velocyto loom files
# generate loom file directory if not existing
if [ ! -d "${scvelo_dir}/loom_files" ] ; then
  	mkdir -p "${scvelo_dir}/loom_files"
fi

# collect loom files
for loom in $(find $scRNA_velocity_cellranger_dir -name \*.loom)
do
	mv $loom ${scvelo_dir}/loom_files
done



echo All done!