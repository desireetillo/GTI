#!/bin/bash

#SBATCH --partition="norm,ccr"
#SBATCH --mem=3g
#SBATCH --mail-type="ALL"
#SBATCH --gres=lscratch:20

module load R
module load bedtools

# EDIT ME
src_directory=/data/GAU/projects/Lebensohn_2021/GT_pipeline/src

ctrl_bed=$1
treat_bed=$2
output_dir=$3

prefix=$(basename $treat_bed .bed)
echo $prefix

mkdir -p $output_dir/Results
perl $src_directory/compute_gene_insertion_enrichment.v3.pl $ctrl_bed $treat_bed $output_dir/Results/${prefix}.gene_level
