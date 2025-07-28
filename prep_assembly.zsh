#!/bin/zsh

# Create a FAI index, calculate summary statistics, gzip the assembly and extract Minia3 coverage statistics

s=(-samtools '/opt/Bio/samtools/1.9/bin/samtools')
n=(-n50 '/filer-dg/agruppen/seq_shared/mascher/code_repositories/genome_assembly_pipeline/shell/n50')

zparseopts -D -K -- -samtools:=s -n50:=n

samtools=$s[2]
n50=$n[2]

$samtools faidx $fa && $n50 $fa.fai > $fa.n50

grep '^>' $fa | tr -d '>' | cut -d ' ' -f 1,5- | awk 'NF > 1' \
 | awk '{for(i=2;i<=NF;i++) {split($i, a, ":"); print $1"\t"a[3]}}' \
 > ${fa:r}_links.tsv &
grep '^>' $fa | tr -d '>' | tr : ' ' \
 | awk '{print $1"\t"$4"\t"$7"\t"$10}' > ${fa:r}_cov.tsv &

wait

gzip $fa
