# Scripts for processing gene trap insertion screens

Built for running on NIH's [biowulf HPC](https://hpc.nih.gov)

## Overview


## Dependencies
   
  * trimmomatic/0.39
  * bowtie/1.2.2
  * samtools/1.13 
  * deeptools/3.5.1
  * ucsc/442
  * bedtools/2.30.0
  * perl
  * R>=3.6
  
     
Additional files:

* Annotation files (bowtie index for hg38, hg38 genome sequence fasta).
* Gene (exon, intron) annotations (bed format) in the annotations/ directory 
* A set of utility scripts in the src/ directory


## Setup instructions

* If setting up for the first time do the following:  
	* Make perl utils executable:
	  	`cd src;
  		chmod+x *pl`
	* Edit path to `src/` directory in the following scripts:
		- line 10 of `Step1_GT_trim_align.sh`
		- line 12 of `Step2_enrichment.sh`
		- line 10 of `src/compute_gene_insertion_enrichment.v3.pl`
	* unzip exons: `cd annotations; gunzip ucsc_refseq_hg38.exons.bed.gz`
	* Edit path to genome fasta (line 14) and bowtie index (line 15) of `Step1_GT_trim_align.sh`

## Steps

**1. Processing of fastqs**

Steps:

* Adapter trims reads using trimmomatic, retains only reads that are >25 bp in length
* Aligns to hg38 using bowtie
* Produces output files for enrichment analysis and viewing on IGV/genome browser (.bed/.bb/.bigwig)
	
Usage:
	`sbatch Step1_align.sh <fastq_file> <output_dir>`

`<fastq_file>`: input fastq file.  
`<output_dir>`: output directory name (best to use same output directory for all samples in a related screen so that all output files are in the same place.

To run processing on all fastqs in a given directory:

```
for i in <fastq_dir>/*gz; do  
	sbatch Step1_GT_trim_align.sh $i <output_dir>
done
```

_List of outputs:_

Output directory | Contents
-----|-----
\<output\_dir>/trim | Trimmed fastqs for each sample
\<output\_dir>/bam/| Alignment bam files
\<output\_dir>/bed/| Genomic locations in bed and bigBed (.bed and .bb) format for each unique insertion
\<output\_dir>/bigwigs/ |  Signal tracks in bigwig (.bw) format that contain the coverage in 50bp bins for all unique insertions for each sample, in which the orientation of the insertion is relative to the human genome
\<output\_dir>/*alignment_stats | alignment statistics

After all fastqs are processed, Alignment statistics can be compiled by running

`bash Step1b_compile_all_stats.sh <output_dir>`
     
where `<output_dir>` is the output directory specified in Step1
     
This will write alignment statistics to a single `<output_dir>.csv` file
     

 
**2. Enrichment analysis**

This step performs enrichment analyses for each GT screen given a control bed file and a sorted bed file.
 
Computes counts for insertion types in control (unsorted) and sorted samples, enrichment FDR corrected pvalues (v1 and v2 p-values), and Intronic Gene Trap Insertion Orientation Bias (IGTIOB) score for the sorted sample.

Score | Description 
-----|-----
v1 p-value | FDR-corrected one-sided Fisher exact test  p-value computed from comparing the relative frequency at which that gene harbored inactivating  GT insertions (sense, antisense exonic + sense intronic) in the sorted cells compared to the relative frequency at which the gene carried any GT insertion in the control dataset.  
v2 p-value | FDR-corrected one-sided Fisher exact test  p-value	 computed from comparing all insertions in sorted cells vs. all insertions in the control population
Intronic Gene Trap Insertion Orientation Bias  (IGTIOB) score | IGTIOB = log2(S/A) x ln(S x A), where ‘S’ and ‘A’ equal one plus the number of unique sense or antisense intronic GT insertions, respectively, in a given gene. High positive IGTIOB scores generally indicate genes whose disruption promotes the phenotype enriched for during the screen.


Output will be writen to tab-delimited text files `<treatment_prefix>`.tab

Script also generates "volcano"-like plots (IGTIOB vs -log10(FDR corrected v2 p-value))

Usage:
      `sbatch Step2_enrichment.sh  <control_bed_file> <sorted_bed_file> <output_directory>`

All outputs will be in `output_directory/Results`


**3. Compile enrichment results**

   Usage:
 
   `python src/files2xlsx.py <output_dir>/Results <output_prefix>`

   Compiles all enrichment results from Step 2 into a single excel file `<output_prefix>.xlsx`
