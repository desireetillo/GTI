#!/bin/bash

aln_out=$1
op=${aln_out/\//}
outfile=$op.stats.csv

echo "Sample,# of reads,# uniquely mapped reads,# unique insertions,coverage/insertion" >$outfile

for stats_file in $aln_out/*alignment_stats; do
    prefix=$(basename $stats_file _alignment_stats)
    bed_file=$aln_out/bed/$prefix.bed
    nreads="$(grep "Number of reads:" $stats_file | cut -f 2 -d ":" | sed 's/^[ \t]*//')"
    areads="$(grep "Number of reads aligned uniquely:"  $stats_file | cut -f 2 -d ":" | sed 's/^[ \t]*//')"
    uins="$(wc -l $bed_file | cut -f 1 -d " ")"
    echo -n "$prefix,$nreads,$areads,$uins," >>$outfile
    echo "scale=2 ; $areads / $uins" | bc >>$outfile
done

echo "wrote stats to $outfile"
