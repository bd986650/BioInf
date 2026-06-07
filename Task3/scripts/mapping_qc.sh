#!/bin/bash
set -e

REF="reference/ecoli.fa"
R1="data/reads_1.fq"
R2="data/reads_2.fq"
OUT_DIR="results"
TMP_DIR="tmp"

mkdir -p "$OUT_DIR" "$TMP_DIR"

echo "=== [1] FASTQC ==="
fastqc "$R1" "$R2" -o "$OUT_DIR"

echo "=== [2] ALIGNMENT (BWA MEM) ==="
bwa mem "$REF" "$R1" "$R2" > "$OUT_DIR/sample.sam"

echo "=== [3] SAM -> BAM ==="
samtools view -bS "$OUT_DIR/sample.sam" > "$OUT_DIR/sample.bam"

echo "=== [4] SORT BAM ==="
samtools sort "$OUT_DIR/sample.bam" -o "$OUT_DIR/sample.sorted.bam"

echo "=== [5] INDEX BAM ==="
samtools index "$OUT_DIR/sample.sorted.bam"

echo "=== [6] FLAGSTAT ==="
samtools flagstat "$OUT_DIR/sample.sorted.bam" > "$OUT_DIR/flagstat.txt"

echo "DONE ✔"
echo "Results in: $OUT_DIR"