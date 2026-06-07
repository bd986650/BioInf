#!/bin/bash

set -euo pipefail

REF="reference/ecoli.fa"
R1="data/SRR2584863_1.fastq"
R2="data/SRR2584863_2.fastq"

OUT="results"

mkdir -p $OUT

echo "=== CHECK INPUT ==="

if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "ERROR: FASTQ files not found!"
    echo "Run: fasterq-dump SRR2584863 --split-files -O data"
    exit 1
fi

echo "=== FASTQC ==="
fastqc $R1 $R2 -o $OUT

echo "=== ALIGNMENT (BWA MEM) ==="
bwa mem $REF $R1 $R2 > $OUT/sample.sam

echo "=== SAM -> BAM ==="
samtools view -bS $OUT/sample.sam > $OUT/sample.bam

echo "=== SORT BAM ==="
samtools sort $OUT/sample.bam -o $OUT/sample.sorted.bam

echo "=== FLAGSTAT ==="
samtools flagstat $OUT/sample.sorted.bam > $OUT/flagstat.txt

echo "=== PARSE MAPPING RATE ==="

MAPPED=$(grep " mapped (" $OUT/flagstat.txt | head -1 | \
sed -E 's/.*\(([0-9.]+)%.*/\1/')

echo "Mapped rate: $MAPPED%"

THRESHOLD=90

echo "=== QC CHECK ==="

if (( $(echo "$MAPPED > $THRESHOLD" | bc -l) )); then
    echo "OK" > $OUT/qc.txt
    echo "RESULT: OK"
else
    echo "not OK" > $OUT/qc.txt
    echo "RESULT: not OK"
fi

echo "=== VARIANT CALLING ==="
freebayes -f $REF $OUT/sample.sorted.bam > $OUT/sample.vcf

echo "DONE ✔"