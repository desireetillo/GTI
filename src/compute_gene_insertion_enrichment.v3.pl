#!/usr/bin/perl
use strict;

# GT insertion analysis script
# 08.09.201908.09.2019, updated 08.29.19
# takes in sorted read bed files from BAIMS script and computes insertions
# then feeds to Rscript to compute p-vals and generate figures

# EDIT ME
my $pipeline_path="/data/GAU/projects/Lebensohn_2021/Gene_Trap_Insertions";
$ENV{PERL_HOME} = "$pipeline_path/src/perl_lib"

my $ctrl_bed = $ARGV[0];
my $sort_bed = $ARGV[1];
my $out_file = $ARGV[2];

# annotation locations
my $exonAnnotations = "$pipeline_path/annotations/ucsc_refseq_hg38.exons.bed";
my $intronAnnotations = "$pipeline_path/annotations/ucsc_refseq_hg38.introns.bed";
my $geneAnnotations = "$pipeline_path/annotations/ucsc_refseq_hg38.genes.bed";

# rscript location
my $rscript="$pipeline_path/src/analyze_GT_data.v3.R";

if(@ARGV == 3){
    my %stats = ();
    
# initialize gene list
   
    open(F, "$geneAnnotations") || die;
    while(my $line = <F>){
	chomp $line;
	my($c,$s,$e,$id,@etc)=split(/\t/, $line);
	$stats{$id}{"Unsorted_Total"}=0;
	$stats{$id}{"Sense_Intron"}=0;
	$stats{$id}{"Antisense_Intron"}=0;
	$stats{$id}{"Sense_Exon"}=0;
	$stats{$id}{"Antisense_Exon"}=0;
	$stats{$id}{"Seen"}=0;    
    }
    
#Unsorted cellsSorted cellsDerived from columns B and CDerived from columns F and G
    
    
#my $ctrl_total = `wc -l $ctrl_bed`;
#my $sort_total = `wc -l $sort_bed`;
    
# FOR Unsorted cells
# All GT insertions in genes
    
    
    my @counts= `intersectBed -a $ctrl_bed -b $geneAnnotations -wb  |  awk  \'\$2 >= \$8\' | $pipeline_path/src/cut.pl -f 10 | sort | uniq -c | $pipeline_path/src/space2tab.pl | $pipeline_path/src/cut.pl -f 3,2`;
    foreach my $x(@counts){
	chomp $x;
	my($i,$c) = split(/\t/, $x);
	$stats{$i}{"Unsorted_Total"}=$c;
	$stats{$i}{"Seen"}++;
    }
    
# sense introns
    my @counts= `intersectBed -a $sort_bed -b $intronAnnotations -wb -s |  awk  \'\$2 >= \$8\' | $pipeline_path/src/cut.pl -f 10 | sort | uniq -c | $pipeline_path/src/space2tab.pl | $pipeline_path/src/cut.pl -f 3,2`;
    foreach my $x(@counts){
	chomp $x;
	my($i,$c) = split(/\t/, $x);
	$stats{$i}{"Sense_Intron"}=$c;
	$stats{$i}{"Seen"}++;
    }
    
# antisense introns
    my @counts= `intersectBed -a $sort_bed -b $intronAnnotations -wb  -S |  awk  \'\$2 >= \$8\' | $pipeline_path/src/cut.pl -f 10 | sort | uniq -c | $pipeline_path/src/space2tab.pl | $pipeline_path/src/cut.pl -f 3,2`;
    foreach my $x(@counts){
	chomp $x;
	my($i,$c) = split(/\t/, $x);
	$stats{$i}{"Antisense_Intron"}=$c;
	$stats{$i}{"Seen"}++;
	
    }
# sense exons
    my @counts= `intersectBed -a $sort_bed -b $exonAnnotations -wb -s |  awk  \'\$2 >= \$8\' | $pipeline_path/src/cut.pl -f 10 | sort | uniq -c | $pipeline_path/src/space2tab.pl | $pipeline_path/src/cut.pl -f 3,2`;
    foreach my $x(@counts){
	chomp $x;
	my($i,$c) = split(/\t/, $x);
	$stats{$i}{"Sense_Exon"}=$c;
	$stats{$i}{"Seen"}++;
	
    }
# antisense exons
    my @counts= `intersectBed -a $sort_bed -b $exonAnnotations -wb  -S |  awk  \'\$2 >= \$8\' | $pipeline_path/src/cut.pl -f 10 | sort | uniq -c | $pipeline_path/src/space2tab.pl | $pipeline_path/src/cut.pl -f 3,2`;
    foreach my $x(@counts){
	chomp $x;
	my($i,$c) = split(/\t/, $x);
	$stats{$i}{"Antisense_Exon"}=$c;
	$stats{$i}{"Seen"}++;
	
    }
    
    
# print out non redundant data
    
    my %nrStats = ();
    my %tmp =();
    
    foreach my $x(sort keys %stats){
	if ($stats{$x}{"Seen"} > 0){ # only print out gene data for which we have non-zero counts	
	    my ($txID, $altName,$coords)  = split(/\;/,$x);
	    my $total_all_insertions = $stats{$x}{"Sense_Exon"} + $stats{$x}{"Antisense_Exon"} + $stats{$x}{"Sense_Intron"} + $stats{$x}{"Antisense_Intron"};
	    my $total_inactivating = $stats{$x}{"Sense_Exon"} + $stats{$x}{"Antisense_Exon"} + $stats{$x}{"Sense_Intron"};
	    my $i_score = compute_igtiob( $stats{$x}{"Sense_Intron"}+1,$stats{$x}{"Antisense_Intron"}+1);
	    my $url = "https:\/\/www.genecards.org\/cgi-bin\/carddisp.pl\?gene=$altName";
	    my $vals = $stats{$x}{"Unsorted_Total"} . "\t" . $total_all_insertions . "\t"  . $total_inactivating . "\t" . $stats{$x}{"Sense_Exon"}  . "\t" .  $stats{$x}{"Antisense_Exon"} . "\t" . $stats{$x}{"Sense_Intron"} . "\t" . $stats{$x}{"Antisense_Intron"} . "\t" . $i_score;
	    my $cmpstring= $altName . "\t" . $url . "\t" .  $vals;
	    push(@{$nrStats{$cmpstring}},$txID);  
	}
    }
    
    my $tmp_file = "$out_file.temp";
    print "Writing temp data to file $tmp_file\n";
    open(TMP, ">$tmp_file") || die  "Couldnt open $tmp_file for writing: $!\n";
    print TMP "Transcripts\tGene\tGeneCardsURL\tInsertionsGenes_Unsorted\tTotalAllInsertions\tTotalInactivating\tSenseExon\tAntisenseExon\tSenseIntron\tAntisenseIntron\tIGTIOB\n";
    foreach my $x(keys %nrStats){
	print TMP  join(";",@{$nrStats{$x}}), "\t", $x , "\n";
    }
    close TMP;

    print "computing P-vals\n";
    print "Rscript --vanilla $rscript $tmp_file $out_file\n";
    system("Rscript --vanilla $rscript $tmp_file $out_file");
    print "removing temp file $tmp_file\n";
    system("rm -f $tmp_file");
    print "done\n";
}
else{
    print "compute_gene_insertion_enrichment.v3.pl ctrl_bed sort_bed out_file\n";
    
}

sub compute_igtiob{
	my($s,$a) = @_;
	my ($score);
	my $x = log2($s/$a);
	my $y= log($s * $a);
	$score = $x * $y;
	return $score;
}

sub log2{
    my($n) = shift;
	return log($n)/log(2);
}

#-s Force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
# -S Require different strandedness. That is, only report hits in B that overlap A on the _opposite_ strand. By default, overlaps are reported without respect to strand.
# FOR SORTED CELLS#

# IGTIOB score=log2(S/A) x ln(SxA) S=1+senseIntronicsertion A=1+antisenseIntronicInsertion
# =LOG2((S+1)/(A+1)) *LN((S+1)*(1+A))

