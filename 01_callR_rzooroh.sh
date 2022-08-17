#!/bin/bash
# header ...

module load R/4.0.0
module load bcftools/1.14
module load parallel/20200522

# Starting with, for each population (4x), for each chromosome (18x), a file .gen.gz and .samples (72x of each)

# Loop to remove 2 extra lines at the top of the .samples files
FILE=($(for i in *samples*
do
	echo $(basename ${i})
done))
for INDS in ${FILE[@]}
do
	tail -n +3 ${INDS} > ${INDS}_temp
	mv ${INDS}_temp ${INDS}
done

# Function to call the R script, with the varible -c being $1
	# e.g. used: run_r_script chromosome_name
function run_r_script {
	Rscript 02_Rscript_rzooroh.R -c $1
}
export -f run_r_script

# Function to get a list of all chromosomes found in the input vcf
	# e.g. used: vcf_chromosomes input.vcf
function vcf_chromosomes {
	bcftools query -f '%CHROM\n' $1 | uniq
}
export -f vcf_chromosomes

# Make an object which is the list of chromosomes
chrom_set=`vcf_chromosomes mLutLut_dp_q.vcf.gz`
# Apply the function to call the R script, in parallel across the list of chromosomes
parallel --verbose -j 18 run_r_script ::: ${chrom_set}
