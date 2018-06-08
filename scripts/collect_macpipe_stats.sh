#!/bin/bash
set -exuo pipefail

# basic trinity stats
N50_f="assembly_stats/"$1"_N50.txt"
ExN50_f="assembly_stats/"$1"_ExN50.txt"
t_genes=$( head -n 6  $N50_f | tail -n 1 | awk '{print $NF}')
t_trans=$( head -n 7  $N50_f | tail -n 1 | awk '{print $NF}')
GC_con=$( awk '/Percent GC:/ {print $NF}' $N50_f )
N50=$( head -n 18  $N50_f | tail -n 1 | awk '{print $NF}')
Median=$( head -n 20  $N50_f | tail -n 1 | awk '{print $NF}')
Ex50N50=$( grep -E "^50" $ExN50_f | awk '{print $2}' )
Ex90N50=$( grep -E "^90" $ExN50_f | awk '{print $2}' )

basic_header="N_genes\tN_transcripts\tGC_contigs\tN50\tMedian_contig_length\tEx50N50\tEx90N50"
basic_stats="$t_genes\t$t_trans\t$GC_con\t$N50\t$Median\t$Ex50N50\t$Ex90N50"

# gather busco results
busco_f="busco/summary/short_summary_"$1".busco.txt"
busco_f_trans="busco/summary/short_summary_"$1".good.busco.txt"
# we need to escape the tabs in all the previous commands and then only replace them at the very end
# otherwise we get a mix of whitespace and proper tabs
busco_stat=$( head -n 8 $busco_f | tail -n 1 |  sed -E 's/C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%,n:([0-9.]+)/\1\\t\2\\t\3\\t\4\\t\5\\t\6/' )



busco_stats_trans=$( head -n 8 $busco_f_trans | tail -n 1 |  sed -E 's/C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\],F:([0-9.]+)%,M:([0-9.]+)%,n:([0-9.]+)/\1\\t\2\\t\3\\t\4\\t\5\\t\6/')

busco_header="busco_complete\tbusco_complete_s\tbusco_complete_d\tbusco_fragmented\tbusco_missing\tbusco_genes\tbusco_good_complete\tbusco_good_complete_s\tbusco_good_complete_d\tbusco_good_fragmented\tbusco_good_missing\tbusco_good_genes"
busco_stats=$( echo ${busco_stat}${busco_stats_trans} )

# transrate
transrate_f="transrate/summary/"$1".assemblies.csv"
transrate_stats=$( tail -n 1 $transrate_f | awk -F "," '{$1=""; print $0}' | sed 's/ /\\t/g' )
transrate_header=$( head -n 1 $transrate_f | awk -F "," '{$1=""; print $0}' | sed 's/ /\\t/g')


# transdecoder

transdecoder_fp="transdecoder/"$1"_longest_orfs.pep.length"

transdecoder_l=$( awk '{print $NF}' $transdecoder_fp | sort -n | awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1< min) {min=$1}; total+=$1; count+=1} END {print total/count, min, max}' | sed 's/ /\\t/g' )
transdecoder_c=$( wc -l "transdecoder/"$1"_longest_orfs.pep" | awk '{print $1}' | sed 's/ /\\t/g' )

transdecoder_header="transdecoder_genes\ttransdecoder_mean_length\ttransdecoder_min_length\ttransdecoder_max_length"
transdecoder_stats=$transdecoder_c"\t"$transdecoder_l


# mapping
mapping_f="mapping/"$1"_snap_to_ref.txt"
mapping_header="snap_n_reads\tsnap_map_high\tsnap_map_low\tsnap_no_map\tsnap_too_short"
mapping_stats=$( tail -n 1 $mapping_f | awk '{print $1"\\t"$2"\\t"$4"\\t"$6"\\t"$8}' | tr -d ',' )

echo -e "assembly\t"$basic_header"\t"$busco_header"\t"$transrate_header"\t"$transdecoder_header"\t"$mapping_header
echo -e $1"\t"$basic_stats"\t"$busco_stats"\t"$transrate_stats"\t"$transdecoder_stats"\t"$mapping_stats


