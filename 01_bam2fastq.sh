#  Script to convert BAM files to FASTQ format using samtools (paired-end reads)
# Usage: ./01_bam2fastq.sh

projectdir='/scratch/GDB136/IPK'
prefix='GDB136'

cd "$projectdir/Hi-Fi" || { echo "ERROR: no existe $projectdir/Hi-Fi"; return 1; }

# Capture what it ends in ".bam"
bam_file=$(ls *.bam 2>/dev/null)
if [ -z "$bam_file" ]; then
  echo "ERROR: No .bam files found in $projectdir/Hi-Fi"
  return 1
fi

echo "Converting BAM file to FASTQ format: $bam_file"
echo "command: samtools fastq ${bam_file} > ${prefix}_HiFi.fastq"

samtools fastq "${bam_file}" > "${prefix}_HiFi.fastq"

# Check if the fastq file was created successfully
if [ ! -s "${prefix}_HiFi.fastq" ]; then
  echo "ERROR: Failed to create FASTQ file from BAM"
  return 1
fi

#come back to the project directory
cd "$projectdir" || { echo "ERROR: no existe $projectdir"; return 1; }
