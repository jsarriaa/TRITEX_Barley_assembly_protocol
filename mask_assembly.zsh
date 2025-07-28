#!/bin/zsh


# Original paths for tools
###########
# samtools='/opt/Bio/samtools/1.9/bin/samtools'
# bbduk='/filer-dg/agruppen/DG/mascher/source/BBMap_37.93/bbduk.sh'
# bit='/filer-dg/agruppen/DG/mascher/source/faToTwoBit'
# tbinfo='/filer-dg/agruppen/DG/mascher/source/twoBitInfo'
# bedtools='/opt/Bio/bedtools/2.27.1/bin/bedtools'
# kmer='/filer-dg/agruppen/DG/mascher/source/BBMap_37.93/kmercountexact.sh'
###########

debug_verbose=true

samtools="/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/samtools"
bedtools="/scratch/software-phgv2/miniconda3/envs/IPK-tritex/bin/bedtools"
bbduk="/scratch/software-phgv2/miniconda3/envs/IPK-tritex-2/bin/bbduk.sh"
bit="/scratch/software-phgv2/miniconda3/envs/IPK-tritex-2/bin/faToTwoBit"
tbinfo="/scratch/software-phgv2/miniconda3/envs/IPK-tritex-2/bin/twoBitInfo"
kmer="/scratch/software-phgv2/miniconda3/envs/IPK-tritex-2/bin/kmercountexact.sh"

# Print everything if debug_verbose:
if [ "$debug_verbose" = true ]; then
  echo "Debug mode is ON"
  echo "samtools: $samtools"
  echo "bedtools: $bedtools"
  echo "bbduk: $bbduk"
  echo "bit: $bit"
  echo "tbinfo: $tbinfo"
  echo "kmer: $kmer"
else
  echo "Debug mode is OFF"
fi



s=(-size 31)
m=(-mincount 2)
e=(-mem '500G')

if [ "$debug_verbose" = true ]; then
  echo "Default parameters:"
  echo "Size: ${s[2]}"
  echo "Min count: ${m[2]}"
  echo "Memory: ${e[2]}"
fi

zparseopts -D -K -- -size:=s -fasta:=f -mincount:=m -out:=d -mem:=e

if [ "$debug_verbose" = true ]; then
  echo "Parsed options:"
  echo "Size: $s"
  echo "Fasta: $f"
  echo "Min count: $m"
  echo "Output directory: $d"
  echo "Memory: $e"
fi

fa=$f[2]
size=$s[2]
mincount=$m[2]
out=$d[2]
mem=$e[2]

if [ "$debug_verbose" = true ]; then
  echo "Final parameters:"
  echo "FASTA file: $fa"
  echo "K-mer size: $size"
  echo "Min count: $mincount"
  echo "Output directory: $out"
  echo "Memory: $mem"
fi

mkdir -p $out
base=$out/${fa:t:r}

if [ "$debug_verbose" = true ]; then
  echo "Base name for output files: $base"
  echo "Creating output directory: $out"
fi

# stop if the input file does not exist
if [ ! -f "$fa" ]; then
  echo "ERROR: Input FASTA file does not exist: $fa"
  exit 1
fi

# stop if there is not base name or the out dir has not bein created
if [ -z "$base" ] || [ ! -d "$out" ]; then
  echo "ERROR: Base name or output directory is not set correctly."
  exit 1
fi

if [ "$debug_verbose" = true ]; then
  echo "Running command: $kmer -Xmx$mem mincount=$mincount in=$fa out=${base}_kmercount.txt k=$size"
fi

$kmer -Xmx$mem mincount=$mincount in=$fa \
 out=${base}_kmercount.txt k=$size \
  > ${base}_kmercount.out 2> ${base}_kmercount.err 

# stop if the kmer count file was not created
if [ ! -f "${base}_kmercount.txt" ]; then
  echo "ERROR: K-mer count file was not created: ${base}_kmercount.txt"
  exit 1
fi

if [ "$debug_verbose" = true ]; then
  echo "Running command: $bbduk kmask='N' k=$size -Xmx$mem in=$fa ref=${base}_kmercount.txt out=${base}_masked.fasta"
fi

$bbduk kmask='N' k=$size -Xmx$mem in=$fa ref=${base}_kmercount.txt \
 out=${base}_masked.fasta > ${base}_bbduk.out 2> ${base}_bbduk.err

# stop if the masked fasta file was not created
if [ ! -f "${base}_masked.fasta" ]; then
  echo "ERROR: Masked FASTA file was not created: ${base}_masked.fasta"
  exit 1
fi

if [ "$debug_verbose" = true ]; then
  echo "Running command: $samtools faidx ${base}_masked.fasta"
fi

$samtools faidx ${base}_masked.fasta 

# stop if the FASTA index file was not created
if [ ! -f "${base}_masked.fasta.fai" ]; then
  echo "ERROR: FASTA index file was not created: ${base}_masked.fasta.fai"
  exit 1
fi


fa="${base}_masked.fasta"

if [ "$debug_verbose" = true ]; then
  echo "Running command: $bbduk -Xmx$mem in=$fa out=${fa:r}_singlecopy.fasta k=$size mincount=$mincount overwrite=true"
fi

$bit -long $fa ${fa:r}.bit

# stop if the twoBit file was not created
if [ ! -f "${fa:r}.bit" ]; then
  echo "ERROR: TwoBit file was not created: ${fa:r}.bit"
  exit 1
fi

if [ "$debug_verbose" = true ]; then
  echo "Running command: $tbinfo -nBed ${fa:r}.bit ${fa:r}_gaps.bed"
fi

$tbinfo -nBed ${fa:r}.bit ${fa:r}_gaps.bed

# stop if the gaps BED file was not created
if [ ! -f "${fa:r}_gaps.bed" ]; then
  echo "ERROR: Gaps BED file was not created: ${fa:r}_gaps.bed"
  exit 1
fi

if [ "$debug_verbose" = true ]; then
  echo "Running command: $bedtools complement -g =(cut -f -2 $fa.fai) -i ${fa:r}_gaps.bed > ${fa:r}_noGaps.bed"
fi

$bedtools complement -g =(cut -f -2 $fa.fai) -i ${fa:r}_gaps.bed > ${fa:r}_noGaps.bed

# stop if the noGaps BED file was not created
if [ ! -f "${fa:r}_noGaps.bed" ]; then
  echo "ERROR: NoGaps BED file was not created: ${fa:r}_noGaps.bed"
  exit 1
fi

# Print a message indicating successful completion
if [ "$debug_verbose" = true ]; then
  echo "Masking and preparation of the reference genome completed successfully."
fi