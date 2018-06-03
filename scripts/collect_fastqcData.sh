#!/bin/bash
set -euo pipefail

file=${1##*/}
file=${file%%_fastqc}
stat_f=$1"/fastqc_data.txt"
# using awk is much more efficient than grep and then awk
n_reads=$( awk '/Total Sequences/ {print $NF}' $stat_f )
n_dedup=$( awk '/#Total Deduplicated Percentage/ {print $NF}' $stat_f | xargs printf "%.*f\n" 2 )
gc=$( awk '/%GC/ {print $NF}' $stat_f )
printf "%s\t%s\t%s\t%s\n" $file  $n_reads $n_dedup $gc
printf "%s\t%s\t%s\t%s\n" $file  $n_reads $n_dedup $gc > "./QC/"$file"_stats.txt"
mv $stat_f "./QC/"$file"_QC.txt" 
rm -drf $1
