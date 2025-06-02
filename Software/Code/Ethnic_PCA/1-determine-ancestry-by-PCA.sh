#!/bin/bash
#SBATCH --job-name=PCA

# Determining ancestry from SNP data 
# Author: Laura Budurlean

# Input: VCF files with SNV's called (using hg38) from your own samples
# Output: The final output will be a plot with your samples overlapped with super-populations from 1000 genomes data (https://www.internationalgenome.org/home) 

# You will need plink2, bcftools, vcftools

# Variables
name=$1
PIPELINEDIR="MY_PIPELINE_DIR" # Replace with your actual pipeline directory path
All_chr_vcf="${PIPELINEDIR}/Data/${name}/VCF/post_imputation/${name}_chr_all.vcf" #post-imputation files
#All_chr_vcf="/storage/home/x_tomisinecvandoren/Ksein/FinalPipelineV3/Data/${name}/VCF/${name}_no_monomorph.vcf" #préimputation  files- useful for testing with small computing power
Analysis_folder="${PIPELINEDIR}/Data/${name}"
PCA_folder="${PIPELINEDIR}/Data/${name}/PCA"
Code_path="${PIPELINEDIR}/Software/Code"


mkdir -p $PCA_folder

PLINK_PREFIX=$PCA_folder"/All"
PLINK_GENOME_PREFIX=$PCA_folder"/All.geno"
SITE2_BED=$PCA_folder"/site2.bed"
THREADS=16
#####################################################
# Start time
START_TIME=$(date +%s)
#####################################################
: '
' #FIN commentairebas
# 1) bgzip and index sample VCF file.
echo 'PARTIE1'

bgzip $All_chr_vcf #--output $All_chr_vcf".gz"
bcftools sort $All_chr_vcf".gz" -Oz -o $All_chr_vcf".gz"
bcftools index $All_chr_vcf".gz" -o $All_chr_vcf".gz.csi"
  
# 2) remove all duplicate mutation sites from Start_file.vcf
echo 'PARTIE2'
echo 'remove all duplicate mutation sites from Start_file.vcf'
bcftools norm --rm-dup both -O z -o $PCA_folder/$name"_norm.vcf.gz" "${All_chr_vcf}.gz"
  
# 3) Convert Start_file_rm_dups.vcf.gz to plink format and also do some more VCF quality filtering. Output files with prefix "All" or another identifier of your choice. Depending on the plink version you are using if you receive a duplicate error  try increasing the default set-missing-var-id character limit from its default 23, using this --new-id-max-allele-len and set it to 100 or even higher until you dont get dups. This is a known issue in plink1.9 and plink2.0 solves it. Read more about it here: https://www.cog-genomics.org/plink/1.9/data#set_missing_var_ids
echo 'PARTIE3'
echo 'Convert Start_file_rm_dups.vcf.gz to plink format and also do some more VCF quality filtering'
plink2 --vcf $PCA_folder/$name"_norm.vcf.gz"\
    --set-missing-var-ids @:#[b19]\$1\$2 \
    --new-id-max-allele-len 250 \
    --allow-extra-chr \
    --chr 1-22 \
    --double-id \
    --max-alleles 2 \
    --vcf-filter \
    --vcf-min-gq 15 \
    --make-bed \
    --out $PLINK_PREFIX

# 4) Remove rare variants from "All" prefix files and create output files with prefix "All.geno" or another identifier of your choice. 
 echo 'PARTIE4'
 plink2 --bfile $PLINK_PREFIX --maf 0.1 --make-bed --out $PLINK_GENOME_PREFIX

    # Troubleshooting Notes: 
    # Check if the All.geno.bim file has "chr" prefixes, if it does remove them.
    # You only need to do this if you see that you have "chr" prefixes in the All.geno.bim file.
    # Create a file with old names and new names and run plink2 with updated names list.
    # Otherwise skip this step and run plink normally without -update-name option.
    # You can try this to update the names if you have "chr" prefixes in your bim file.
        #	awk '{print $2}' All.geno.bim > tmp1
        #	sed 's/chr//g' tmp1 > tmp2
        #	paste tmp1 tmp2 > updatenames.list
        #	plink2 --bfile All.geno --update-name updatenames.list --snps-only --make-bed --out All.geno.snps

# 5) Retain SNPs only, bgzip, and index. This file is what contains your samples variants that will be input into PCA later.
 echo 'PARTIE5'
 plink2 --bfile $PLINK_GENOME_PREFIX --snps-only --make-bed --recode vcf --out "${PLINK_GENOME_PREFIX}.snps"
 bgzip "${PLINK_GENOME_PREFIX}.snps.vcf"
 bcftools index "${PLINK_GENOME_PREFIX}.snps.vcf.gz"

# 6) Create "site2.bed" file which is SNP positions present in your samples, that you will then use to filter out variants in the 1000 genomes samples. You need to make 2 files (1 for your samples and 1 for 1000 genomes) that only have the SNPs present in both files.
 echo 'PARTIE6'
 awk -F '\t' '{print $1"\t"$4-1"\t"$4}' "${PLINK_GENOME_PREFIX}.snps.bim" > "${SITE2_BED}"

# 7) Filter your sites against the 1000 genome sites. You will need phased shapeit2 files from 1000genomes "shapeit2_integrated..." (http://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/)

#Set up an array as a bash script for every chromosome.

echo 'PARTIE7'
#séparation en 22 chromosomes
for number in $(seq 1 22); do
    bcftools filter --threads $THREADS -T $SITE2_BED -Oz -o $PCA_folder"/chr${number}.1kg.sites.vcf.gz" "/storage/home/x_tomisinecvandoren/1K_Genome_hg19/ALL.chr${number}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
done 
#Check threads


#remove
#number=13
#bcftools filter --threads $THREADS -T $SITE2_BED -Oz -o $PCA_folder"/chr${number}.1kg.sites.vcf.gz" "/storage/home/x_tomisinecvandoren/1K_Genome_hg19/ALL.chr${number}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

# 8) Combine all the chromosome files using vcftools, bgzip and index.
 echo 'PARTIE8'
 vcf-concat $PCA_folder"/"chr*.1kg.sites.vcf.gz > $PCA_folder"/merged.1kg.sites.vcf"
 bgzip $PCA_folder"/merged.1kg.sites.vcf"
 bcftools index "${PCA_folder}/merged.1kg.sites.vcf.gz"

# 9) Intersect the All.geno.snps.vcf.gz and merged.1kg.sites.vcf.gz. Set the # of threads as needed
 echo 'PARTIE9'
 bcftools isec --threads ${THREADS} -p "${PCA_folder}/isec" -n =2 "${PLINK_GENOME_PREFIX}.snps.vcf.gz" "${PCA_folder}/merged.1kg.sites.vcf.gz" #-o $PCA_folder
 
 bgzip  "${PCA_folder}/isec/0000.vcf"
 bgzip  "${PCA_folder}/isec/0001.vcf"
 bcftools index "${PCA_folder}/isec/0000.vcf.gz"
 bcftools index "${PCA_folder}/isec/0001.vcf.gz"
 bcftools merge -m none "${PCA_folder}/isec/0000.vcf.gz" "${PCA_folder}/isec/0001.vcf.gz" -Oz -o "${PCA_folder}/merged_ALL_1kg_samples.vcf.gz"

# 10) Set missing genotypes to 0|0
 #bcftools +setGT merged_ALL_1kg_samples.vcf.gz -- -t . -n 0p\n > merged_ALL_1kg_samples_setGT.vcf.gz
 echo 'PARTIE10'
 plink2 --vcf $PCA_folder/merged_ALL_1kg_samples.vcf.gz --set-missing-var-ids @:#[b19]\$r\$a --make-bed --out $PCA_folder/unpruned1kg_setGT
 plink2 --bfile $PCA_folder/unpruned1kg_setGT --indep-pairwise 100kb 1 .1 --out $PCA_folder
 plink2 --bfile $PCA_folder/unpruned1kg_setGT --extract $PCA_folder/plink2.prune.in --maf 0.1 --pca biallelic-var-wts --make-bed --out $PCA_folder/All.geno.prune.setGT
 
# 11) PCA scoring. You may score as many principal components as needed. I recommend at least the first three for samples that map ambiguously when graphing PC1 and PC2. You may wish to graph PC1 and PC3 as well. Keeping in mind that the first two principal components will explain the most amount of variation in the data.
 echo 'PARTIE11'
 plink2 --bfile $PCA_folder/unpruned1kg_setGT --score $PCA_folder/All.geno.prune.setGT.eigenvec.var 2 3 5 --out $PCA_folder/1kg-PC1-score
 plink2 --bfile $PCA_folder/unpruned1kg_setGT --score $PCA_folder/All.geno.prune.setGT.eigenvec.var 2 3 6 --out $PCA_folder/1kg-PC2-score
 plink2 --bfile $PCA_folder/unpruned1kg_setGT --score $PCA_folder/All.geno.prune.setGT.eigenvec.var 2 3 7 --out $PCA_folder/1kg-PC3-score

# Continue in R to plot samples with 1000 genomes super-population data. See "plot.R" script.
echo 'PARTIE12'
Rscript ${Code_path}/Ethnic_PCA/2-plot.R $PCA_folder $Code_path/Ethnic_PCA/sample_info_file-ADD_SAMPLES.csv $Analysis_folder"/VCF/${name}_no_monomorph.vcf"

# End time
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
# Print runtime
echo "Total runtime: $RUNTIME seconds"
