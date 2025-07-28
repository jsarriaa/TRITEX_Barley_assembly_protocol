#!/bin/zsh
# Align Hi-C reads to a genome assembly and assign read pairs to restriction fragments

###########
# OLD ORIGINAL PATHS (for reference, not used)
# mini=(-minimap ...)
# cut=(-cutadapt ...)
# sam=(-samtools ...)
# bedt=(-bedtools ...)
# novo=(-novosort ...)
# bgz=(-bgzip ...)
# indexsize='50G'
# batchmem='5G'
###########

add_dist='
BEGIN{
 OFS="\t"
}
$9 == "-" { d1=$2-$12 }
$10 == "-" { d2=$5-$15 }
$9 == "+" { d1=$13-$3 }
$10 == "+" { d2=$16-$6 }
d1 < 0 { d1 = 0 }
d2 < 0 { d1 = 0 }
$1 == $4 {
 if($12 == $15)
  rd = 0
 else if($13 < $15)
  rd = $15-2-($13+3)
 else
  rd = $12-2-($16+3)
}
$1 != $4 { rd=-1 }
{
 l1=d1+$3-$2
 l2=d2+$6-$5
 print $0"\t"d1,d2,l1,l2,l1+l2,rd
}'

dist_stat='
BEGIN{
 OFS="\t"
}
$22 == 0 { same++ }
$22 == -1 { diffcl++ }
#$22 == 1 { adj++ }
$22 > 1 { samecl++ }
END {
 print NR,0+same,0+diffcl,0+samecl
}'

# Parse options
zparseopts -D -K -- -i:=i -indir:=ind -linker:=l -minlen_read:=m -threads:=t \
-mem:=e -bed:=b -ref:=r -tmp:=p -minq:=q -maxlen:=x -onlyfrag:=o \
-cutadapt:=cut -minimap:=mini -samtools:=sam \
-bedtools:=bedt -novosort:=novo -bgzip:=bgz

# Remove all parsed options from $@, leaving only positional arguments
while [[ "$1" == -* ]]; do
  shift
done

# Assign input_dir from -i, -indir, or positional
if [[ -n "${i[2]}" ]]; then
  input_dir="${i[2]}"
elif [[ -n "${ind[2]}" ]]; then
  input_dir="${ind[2]}"
elif [[ -n "$1" ]]; then
  input_dir="$1"
else
  input_dir=""
fi

# Assign all other variables from parsed options or defaults
minimap=${mini[2]:-'/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/minimap2'}
cutadapt=${cut[2]:-'/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/cutadapt'}
samtools=${sam[2]:-'/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/samtools'}
bedtools=${bedt[2]:-'/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/bedtools'}
novosort=${novo[2]:-'none'}
bgzip=${bgz[2]:-'/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/bgzip'}

linker=${l[2]:-"AAGCTAGCTT"}  # Use the proper linker for your restriction enzyme, e.g., GATCGATC for DpnII
minlen=${m[2]:-30}
threads=${t[2]:-1}
mem=${e[2]:-"5G"}
ref=${r[2]:-""}
tmp=${p[2]:-"/tmp"}
minq=${q[2]:-10}
maxlen=${x[2]:-500}
bed=${b[2]:-""}
onlyfrag=${o[2]:-0}

echo "Parsed options:"
echo "Linker: $linker"
echo "Min length: $minlen"
echo "Threads: $threads"
echo "Memory: $mem"
echo "Reference: $ref"
echo "Temporary directory: $tmp"
echo "Min quality: $minq"
echo "Max length: $maxlen"
echo "Only fragments: $onlyfrag"
echo "Input directory: $input_dir"

# Validate arguments
if [ -z "$input_dir" ]; then
  echo "ERROR: Input directory is not set."
  exit 1
fi

if [ -z "$threads" ] || [ "$threads" -eq 1 ]; then
  echo "ERROR: Number of threads (-threads) is not set or is 1."
  exit 1
fi

if [ -z "$ref" ]; then
  echo "ERROR: Reference file (-ref) is not set."
  exit 1
fi

if [ -z "$bed" ]; then
  echo "ERROR: BED file (-bed) is not set."
  exit 1
fi

if [[ $threads -gt 16 ]]; then
  novosort_threads=16
else
  novosort_threads=$threads
fi

base=$input_dir/${input_dir:t}
bam=${base}.bam
minimaperr=${base}_minimap.err
samtoolserr=${base}_samtools.err
sorterr1=${base}_samsort1.err
sorterr2=${base}_samsort2.err
cutadapterr=${base}_cutadapt.err
sorterr0=${base}_samtools0.err

b=$base:t
rgentry="@RG\tID:$b\tPL:ILLUMINA\tPU:$b\tSM:$b"

# Note the find commands. If sample names contain R1 and R2, this will not work.
if [[ $onlyfrag -eq 0 ]]; then
  file=$(find $input_dir | grep _R1 | grep 'f.*q.gz$')
  echo -e "poook$input_dir,....$file"
  $cutadapt -f fastq --interleaved -a $linker -A $linker -O 1 -m $minlen 2> $cutadapterr \
    <(find $input_dir | grep _R1 | grep 'f.*q.gz$' | sort | xargs zcat) \
    <(find $input_dir | grep _R2 | grep 'f.*q.gz$' | sort | xargs zcat) \
    | $minimap -ax sr -R $rgentry -t $threads -2 -I $indexsize -K$batchmem $ref /dev/stdin 2> $minimaperr \
    | $samtools fixmate -m - - 2> samtools_fixmate.err \
    | $samtools sort -m 6G -@ 8 2> samtools_sort_pos.err \
    | $samtools markdup -r - - 2> samtools_markdup.err \
    | $samtools sort -@ $novosort_threads -m $mem -n -o $bam 2> samtools_sort_name.err

  echo $pipestatus > ${base}_pipestatus.txt
fi

f="${base}_reads_to_fragments.bed.gz"
err="${base}_reads_to_fragments.err"

{
  {
    $samtools view -H $bam
    $samtools view -q 10 -F 3332 $bam \
      | tee >(cut -f 1 | uniq -c | awk '$1 == 2' | wc -l > ${base}_both_mapped_q10.len) \
      | cut -f -16 \
      | awk 'old != $1 {old=$1; printf "\n"$0; next} {printf "\t"$0}' \
      | awk -F '\t' 'NF == 32' \
      | awk '{for(i=1; i<=15; i++) printf $i"\t"; printf $16"\n"$17; for(i=18; i<=32; i++) printf "\t"$i; print ""}'
  } | $samtools view -Su - \
    | $bedtools pairtobed -f 1 -type both -abam - -bedpe -b $bed \
    | awk 'old != $7 {old=$7; printf "\n"$0; next} {printf "\t"$0}' \
    | awk 'NF == 26'  | cut -f 1-13,24-26 | awk "$add_dist" \
    | tee >($bgzip -c > $f) \
          >(wc -l > $f:r.len) \
          >(awk -v maxlen=$maxlen '$21 > maxlen && $1 == $4' \
             | cut -f 1,2,3,5,6 | awk '{print $5 - $2}' \
             | awk -v maxlen=$maxlen '$1 <= maxlen' | sort | uniq -c > ${base}_length_dist_PE.txt ) \
          >(awk -v maxlen=$maxlen '$21 > maxlen' | wc -l > ${base}_pe_count.txt) \
          >(awk -v maxlen=$maxlen '$21 <= maxlen' | cut -f 21 | sort | uniq -c > ${base}_length_dist.txt) \
          >(awk -v maxlen=$maxlen '$21 <= maxlen' | awk "$dist_stat" > ${base}_frag_stat.txt) \
    | awk -v maxlen=$maxlen '$21 <= maxlen' \
    | cut -f 11,12,14,15 | awk '$1 != $3 || $2 != $4' \
    | awk '{print $0"\n"$3"\t"$4"\t"$1"\t"$2}' \
    | $bgzip > ${base}_fragment_pairs.tsv.gz
} 2> $err







############ THISIS A MESS REWRITE THE FCKING SCRIPT
## REMAKING BUT IT BASH