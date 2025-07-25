# Script to split a FASTQ into different files
# Usage: ./01.5_splitFASTQ.sh <input.fastq>

input_fastq="$1"
input_folder ="$(dirname "$input_fastq")"

if [ -z "$input_fastq" ]; then
  echo "Usage: $0 <input.fastq>"
  return 1
fi

parts=20

seqkit split -p $parts "${input_fastq}"

# Check if in the folder there are $parts files
# Folder = input_fastq.split

split_dir="${input_fastq}.split"
if [ ! -d "$split_dir" ]; then
  echo "ERROR: Split directory $split_dir does not exist."
  return 1
fi  

split_files=$(ls "$split_dir" | wc -l)
if [ "$split_files" -ne "$parts" ]; then
  echo "ERROR: Expected $parts files in $split_dir, but found $split_files files."
  return 1
fi  

echo "Successfully split $input_fastq into $parts"

# Move files to original folder

for file in "$split_dir"/*; do
  mv "$file" "$input_folder"
done    

# Remove the split directory
rm -r "$split_dir"

# Mask original file
mv "$input_fastq" "${input_fastq}.raw"
echo "Original file $input_fastq has been masked as ${input_fastq}.raw"
