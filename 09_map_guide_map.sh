# Maping marker sequences to contig assembly
# Output PAF file will be input into the next R step

#!/bin/bash

debug_verbose=true

projectdir='/scratch/GDB136/IPK'
projectname='GDB136'
refname='MorexV3' # Reference genome name
minimap2=$(which minimap2)

if [ "$debug_verbose" = true ]; then
  echo "Debug mode is ON"
  echo "Project directory: $projectdir"
  echo "Project name: $projectname"
  echo "Reference genome name: $refname"
else
  echo "Debug mode is OFF"
fi

ref="$projectdir/assembly/$projectname/$projectname.p_ctg.fa"
qry="${refname}_singlecopy_100bp.fasta" # this is the guide map sequence filesize='20G'
mem='5G'
threads=30
prefix=$(basename "$ref" .fa)_"$projectname" 

if [ "$debug_verbose" = true ]; then
  echo "Reference FASTA file: $ref"
  echo "Query FASTA file: $qry"
  echo "Size: $size"
  echo "Memory: $mem"
  echo "Threads: $threads"
  echo "Output prefix: $prefix"
fi

if [ "$debug_verbose" = true ]; then
  echo "Running the command: $minimap2 -t $threads -2 -I $size -K$mem -x asm5 $ref $qry | bgzip > $prefix.paf.gz"
fi

{
$minimap2 -t $threads -2 -I $size -K$mem -x asm5 $ref $qry | bgzip > $prefix.paf.gz
} 2> $prefix.err 