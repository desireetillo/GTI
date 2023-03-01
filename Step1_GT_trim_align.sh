#!/bin/bash

#SBATCH --partition="ccr"
#SBATCH --mem=24g
#SBATCH --mail-type="ALL"
#SBATCH --cpus-per-task=16
#SBATCH --time=8:00:00

# src dir (EDIT THIS)
pipeline_path="/data/GAU/projects/Lebensohn_2021/Gene_Trap_Insertions/src"


# required files:
genome_seq="/fdb/igenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
bowtie_index_path="/data/GAU/apps/baims_pipeline/baims_pipeline_asset/resources/GRCh38_index/GCA_000001405.15_GRCh38_no_alt_analysis_set"
chromsizes="annotations/GRCh38.chrom.sizes"

# inputs
fastq=$1  # full path to fastq file
outdirname=$2
prefix=$(basename $fastq _R1_001.fastq.gz)

# outputs
workdir=`pwd`
outdir="$workdir/$outdirname"
mkdir -p $outdir;
mkdir -p $outdir/bam
mkdir -p $outdir/bed
mkdir -p $outdir/trim
mkdir -p $outdir/bigwigs

# output files
trimmed_fastq="$outdir/trim/trim_plus25.$prefix.fastq.gz"
samfile="$outdir/bam/$prefix.sam"
bedOutput="$outdir/bed/$prefix.bed"
bbOutput="$outdir/bed/$prefix.bb"
bamfile="$outdir/bam/$prefix.bam"
statsfile="$outdir/${prefix}_alignment_stats"
insfile="$outdir/${prefix}.insert_coverage.tab"
out_bw_plus="$outdir/bigwigs/$prefix.plus.bw"
out_bw_minus="$outdir/bigwigs/$prefix.minus.bw"

# trim
module load trimmomatic/0.39;
java -jar $TRIMMOJAR SE -threads $SLURM_CPUS_PER_TASK $fastq $trimmed_fastq -trimlog $trimmed_fastq.log ILLUMINACLIP:/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 MINLEN:25

# align
module load bowtie/1.2.2;
bowtie_path="/usr/local/apps/bowtie/1.2.2/bin"
echo "Path to the Bowtie folder is $bowtie_path"
echo "Beginning bowtie alignment for $trimmed_fastq"
${bowtie_path}/bowtie -S --sam-nohead -p $SLURM_CPUS_PER_TASK -v 3 -k 5 --best ${bowtie_index_path} $trimmed_fastq $samfile

# get alignment stats
echo "Computing alignment stats"
#module load python/2.7
python $pipeline_path/initial_sam_analysis.py -i $samfile -o $statsfile                           
python $pipeline_path/samToBED_print_all.py -i $samfile | sort -k1,1 -k2,2n  |  cut -f 1,2 | awk '{ cnts[$0] += 1 } END { for (v in cnts) print v, cnts[v] }' | tr '\t' ':' |  ${pipeline_path}/space2tab.pl | cut -f 1,2 > $insfile

# create bed and bb files
echo "Creating bed files"

python $pipeline_path/samToBED.py -i $samfile -o tmp.$prefix.bed
sort -k1,1 -k2,2n tmp.$prefix.bed >$bedOutput

module load ucsc/442
cat $bedOutput | ${pipeline_path}/add_column.pl -s "$prefix" | ${pipeline_path}/add_column.pl -s "1000" | ${pipeline_path}/uniquify.pl -f -c 6  | ${pipeline_path}/${pipeline_path}/cut.pl -f 1-3,7,8,6 >tmp.$prefix.bed
bedToBigBed tmp.$prefix.bed $chromsizes $bbOutput
rm -f tmp.$prefix.bed

# convert alignment to bam file and generate stranded bigwigs
## rehead samfile, generate bam file
echo "Creating bams and  stranded bigwigs"
module load samtools/1.13
echo "samtools view -h  -T $genome_seq $samfile |  samtools view -S -b -F 4 - >tmp.$prefix.bam"
samtools view -h  -T $genome_seq $samfile |  samtools view -S -b -F 4 - >tmp.$prefix.bam

## sort
samtools sort -@$SLURM_CPUS_PER_TASK tmp.$prefix.bam  >$bamfile 
samtools index $bamfile
rm -f tmp.$prefix.bam $samfile

# generate stranded bigwigs
module load deeptools/3.5.1
bamCoverage --bam $bamfile -o $out_bw_plus  --samFlagExclude 16 --ignoreDuplicates -p max --scaleFactor 1
bamCoverage --bam $bamfile -o $out_bw_minus --samFlagInclude 16 --ignoreDuplicates -p max --scaleFactor -1
