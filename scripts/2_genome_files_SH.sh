#!/bin/bash

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=0-02:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=12




#############################################################
####                                                    #####
####              User defined variables                #####
####                                                    #####
#############################################################



############################################################



module load samtools/1.17
module load cellranger/7.2.0



### generate genome files
### download human or mouse cellranger reference data
if [ $species = "human" ] ; then
	if [ ! -d ${genome_dir}/refdata-cellranger-GRCh38-3.0.0 ] ; then
		cd $genome_dir
		wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
		tar -zxvf refdata-gex-GRCh38-3.0.0.tar.gz
		rm refdata-gex-GRCh38-3.0.0.tar.gz
		cd $script_dir
	fi
elif [ $species = "mouse" ] ; then
	if [ ! -d ${genome_dir}/refdata-gex-mm10-2020-A ] ; then
		cd $genome_dir
		wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
		tar -zxvf refdata-gex-mm10-2020-A.tar.gz
		rm refdata-gex-mm10-2020-A.tar.gz
		cd $script_dir
	fi
fi



### download zebrafish genome files and prepare cellranger reference data
if [ $species = "zebrafish" ] ; then
	# download fasta file if not present
	if [ ! -f "${genome_dir}/${genome}.fa" ]; then
		rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/${species_latin}/dna/${species_latin_2}.${genome}.dna.primary_assembly.fa.gz $genome_dir
		mv ${genome_dir}/${species_latin_2}.${genome}.dna.primary_assembly.fa.gz ${genome_dir}/${genome}.fa.gz
		gunzip ${genome_dir}/${genome}.fa.gz
	fi


	# generate chromosome sizes file if not present
	#if [ ! -f "${genome_dir}/${genome}.chrom.sizes" ]; then
	#	awk -F '\t' '{print $1, $2}' ${genome_dir}/${genome}.fa.fai > ${genome_dir}/${genome}.chrom.sizes
	#fi


	# download gtf file if not present
	if [ ! -f "${genome_dir}/${genome}.gtf" ]; then
		rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/gtf/${species_latin}/${species_latin_2}.${genome}.110.gtf.gz $genome_dir
		mv ${genome_dir}/${species_latin_2}.${genome}.110.gtf.gz ${genome_dir}/${genome}.110.gtf.gz
		gunzip ${genome_dir}/${genome}.110.gtf.gz
	fi


	# prepare cellranger reference data
	cellranger mkgtf \
    		${genome_dir}/${genome}.110.gtf \
    		${genome_dir}/${genome}.110.filtered.gtf \
    		--attribute=gene_biotype:protein_coding

	cellranger mkref \
    		--genome=refdata-cellranger-${genome} \
    		--fasta=${genome_dir}/${genome}.fa \
    		--genes=${genome_dir}/${genome}.110.filtered.gtf \
		--output-dir=$genome_dir
fi



### download file containing positions of repetitive elements in genome
if [ ! -f "${genome_dir}/${genome}_rmsk.txt" ]; then
	rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/${genome_ucsc}/database/rmsk.txt.gz $genome_dir
	wait
	gunzip ${genome_dir}/rmsk.txt.gz
	mv ${genome_dir}/rmsk.txt ${genome_dir}/${genome}_rmsk.txt
fi



echo All done!
