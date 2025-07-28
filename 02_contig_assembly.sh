projectdir='/scratch/GDB136/IPK'

cd "$projectdir/assembly" || { echo "ERROR: no existe $projectdir/assembly"; return 1; }

d="$projectdir/Hi-Fi"
if [ ! -d "$d" ]; then
  echo "ERROR: no existe $d"
  return 1
fi

find "$d" -type f -name '*.fastq' > hifi_reads.txt

if [ ! -s hifi_reads.txt ]; then
  echo "ERROR: No .fastq.gz files found in $d"
  cd ..
  return 1
fi

reads='hifi_reads.txt'
prefix='GDB136'
threads=32

mkdir -p "$prefix"
out="$prefix/$prefix"

echo "Running hifiasm with the following parameters:"
echo "Reads file: $reads"
echo "Output prefix: $out"
echo "Threads: $threads"


xargs -a "$reads" hifiasm -t "$threads" -o "$out" > "$out.out" 2> "$out.err"

echo "runned the command: hifiasm -t $threads -o $out > $out.out 2> $out.err"
