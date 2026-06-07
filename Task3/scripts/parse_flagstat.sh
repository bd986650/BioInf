#!/bin/bash
set -e

FILE=$1

if [ -z "$FILE" ]; then
    echo "Usage: bash parse_flagstat.sh results/flagstat.txt"
    exit 1
fi

echo "Parsing: $FILE"

# total reads
TOTAL=$(grep "in total" "$FILE" | awk '{print $1}')

# mapped reads (first occurrence line with "mapped (")
MAPPED=$(grep "mapped (" "$FILE" | head -n1 | awk '{print $1}')

if [ -z "$TOTAL" ] || [ -z "$MAPPED" ]; then
    echo "ERROR: cannot parse flagstat file"
    exit 2
fi

# percent calculation
PERCENT=$(awk "BEGIN {printf \"%.2f\", ($MAPPED/$TOTAL)*100}")

echo "Total reads   : $TOTAL"
echo "Mapped reads  : $MAPPED"
echo "Mapped %      : $PERCENT%"

# QC rule
THRESHOLD=90

echo "Threshold     : $THRESHOLD%"

if (( $(awk "BEGIN {print ($PERCENT > $THRESHOLD)}") )); then
    echo "QC RESULT     : OK ✔"
    exit 0
else
    echo "QC RESULT     : NOT OK ✖"
    exit 1
fi