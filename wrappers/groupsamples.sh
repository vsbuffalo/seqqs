#!/bin/bash
# groupsamples.sh - group samples across set value from Casava pipeline
# Vince Buffalo <vsbuffaloAAAAAA@gmail.com> (sans poly-A)
set -e
set -u
set -o pipefail

VERSION=0.1

usage() {
    echo -e "\
usage: groupsamples.sh directory [pair]\n\
	directory	directory to search for files\n\
	pair		either 1 or 2, indicating which pair to grab.\n\
Note: all results sorted by name, and paired end reads are checked for consistency.\n\
trim.sh version $VERSION" >&2
    exit 1
}

if [ $# -lt 1 ]; then
    usage;
fi

DIR=$1

if [ $# -eq 1 ]; then
    # single ended mode
    find $DIR -name "*.fastq.gz" | sort
elif [ $# -gt 1 ]; then
    # paired ended mode
    PAIR=$2
    PAIRS_1=$(find $DIR -name "*.fastq.gz" -and -name "*_R1_*" | sort | wc -l)
    PAIRS_2=$(find $DIR -name "*.fastq.gz" -and -name "*_R2_*" | sort | wc -l)
    
    # check same length of paired files
    if [ ! "$PAIRS_1" == "$PAIRS_2" ]; then
	echo "[groupsamples.sh] error: different number of R1 and R2 samples." >&2
	exit 1
    fi
    
    # further checking by taking set variable for each file and
    # checking there are two
    TMP=$(find $DIR -name "*.fastq.gz" -and -name "*_R?_*" | sed 's/.*_\([0-9]*\).fastq.gz/\1/' | sort | uniq -c | sed 's/ *\([0-9]*\) [0-9]*/\1/' | uniq -c  | wc -l)
    if [ $TMP -gt 1 ]; then
	echo "[groupsamples.sh] error: non-matching set entries for R1 and R2 samples." >&2
	exit 1
    fi
    
    # everything checks out, output the correct entry.
    if [ "$PAIR" -eq 1 ]; then 
	find $DIR -name "*.fastq.gz" -and -name "*_R1_*" | sort
    elif [ "$PAIR" -eq 2 ]; then
	find $DIR -name "*.fastq.gz" -and -name "*_R2_*" | sort
    fi
else
    usage
fi