# Scripts for processing gene trap insertion screens

Built for running on NIH's biowulf HPC

## Dependencies
   
  * trimmomatic/0.39
  * bowtie/1.2.2
  * samtools/1.13 
  * deeptools/3.5.1
  * ucsc/418
  * bedtools/2.30.0
  * python/2.7
  * perl/5.24
  * R>=3.6
  
     
Additional files:

* annotation files (bowtie index for hg38, hg38 genome sequence fasta).
* gene (exon, intron) annotations (bed format) in the annotations/ directory 
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

## Steps

**1. Processing of fastqs**

Steps:

* adapter trims reads using trimmomatic, retains only reads that are >25 bp in length
* aligns to hg38 using bowtie
* alignments are converted into genome track formats:
* locations of all unique insertions in bed files (.bed and .bb) format 
* bigwigs that contain the coverage in 50bp bins of all unique insertions for each sample, in which the orientation of the insertion is relative to the human genome.
	
Usage:
	`sbatch Step1_align.sh <fastq_file> <output_dir>`

`<fastq_file>`: input fastq file.  
`<output_dir>`: output directory name (best to use same output directory for all samples in a related screen so that all output files are in the same place.
	
After all fastqs are processed, alignment statistics can be compiled by running

`bash Step1b.sh <output_dir>`
     
where `<output_dir>`: output_directory specified in Step1
     
Alignment statistics into a single `<output_dir>.csv` file
     

 
**2. Enrichment analysis**

Performs enrichment analyses for each GT screen given a control bed file and a treatment bed file
 
- computes counts for insertion types in control and treatment samples, enrichment FDR corrected pvalues (v1 and v2 p-values), and IGTIOB scores
- statistics in `<treatment_prefix>`.tab files
- also generates "volcano" plots (IGTIOB vs -log10(FDR corrected v2 p-value))

Usage:
      `sbatch Step2_enrichment.sh  <control_bed_file> <treatment_bed_file> <output_directory>`

All outputs will be in `output_directory/Results`

**3. Compile outputs.**

   Usage:
 
   `python src/files2xlsx.py <output_dir>/Results <output_prefix>`

   Compiles all enrichment results from step 2 into a single excel file `<output_prefix>.xlsx`
