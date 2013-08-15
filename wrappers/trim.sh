#!/bin/bash
# trim.sh - generic, slightly insane paired-end quality trimming script
# Vince Buffalo <vsbuffaloAAAAAA@gmail.com> (sans poly-A)
set -e
set -u
set -o pipefail

VERSION=0.3

usage() {
    echo -e "\
usage: trim.sh sample_name in1.fq in2.fq [adapters prior qual_thresh]\n\
	sample_name	sample name for both paired-end files \n\
	in1.fq.ggz     	paired-end file 1\n\
	in2.fq.gz      	paired-end file 2\n\
\n\
	Optional arguments:\n\
	adapters	adapters file for Scythe (FASTA file, default: illumina_adapters.fa)\n\
	prior		prior contamination rate (default: 0.4)\n\
	qual_thresh	minimum quality for window (default: 20)\ntrim.sh version $VERSION" >&2
    exit 1
}

if [ $# -lt 3 ]; then
    usage;
fi


## pre-config
# replace $1, $2, $3 with your sample names or use command line
# arguments
SAMPLE_NAME=$1
IN1=$2
IN2=$3

## Trimming pipeline presets
# path to illumina adapters file - replace this with a path to your
# adapters file
ADAPTERS=${4:-illumina_adapters.fa}

# scythe's prior
PRIOR=${5:-0.4}

# sickle's quality threshold
QUAL_THRESH=${6:-20}

check_file_exists() {
    if [ ! -f "$1" ]; then 
	echo "[trim.sh] error: file '$1' does not exist." >&2
	exit
    fi
}

## Check that programs are in $PATH
((command -v scythe && command -v sickle && command -v seqqs) > /dev/null) || \
    (echo "[trim.sh] error: either scythe, sickcle, or seqqs is not in your $PATH" && exit 1)

## Check input
check_file_exists $ADAPTERS
check_file_exists $IN1
check_file_exists $IN2

if [ ! ${#SAMPLE_NAME} -gt 0 ]; then
    echo "[trim.sh] error: specify sample name" >&2
    exit
fi

echo "[trim.sh] running sample '$SAMPLE_NAME' with input PE files '$IN1'/'$IN2'..."

# time trimming process process
T="$(date +%s)"

sickle pe -t sanger -q $QUAL_THRESH \
    -f <(seqqs -e -p raw_${SAMPLE_NAME}_R1 "$IN1" | scythe -a "$ADAPTERS" -p $PRIOR - 2> ${SAMPLE_NAME}_R1_scythe.stderr) \
    -r <(seqqs -e -p raw_${SAMPLE_NAME}_R2 "$IN2" | scythe -a "$ADAPTERS" -p $PRIOR - 2> ${SAMPLE_NAME}_R2_scythe.stderr) \
    -o >(seqqs -e -p trimmed_${SAMPLE_NAME}_R1 - | gzip > ${SAMPLE_NAME}_R1_trimmed.fq.gz) \
    -p >(seqqs -e -p trimmed_${SAMPLE_NAME}_R2 - | gzip > ${SAMPLE_NAME}_R2_trimmed.fq.gz) \
    -s >(seqqs -e -p trimmed_${SAMPLE_NAME}_singles - | gzip > ${SAMPLE_NAME}_singles_trimmed.fq.gz) > ${SAMPLE_NAME}_sickle.stderr

T="$(($(date +%s)-T))"
echo "[trim.sh] $SAMPLE_NAME took seconds: ${T}"
