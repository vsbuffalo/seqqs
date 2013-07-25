# Wrappers Scripts

`seqqs` (and the tool `pairs` packaged with it, to interleave and
uninterleave files) is wrapped in a simple Python script. Currently
only paired-end operation is supported.

## qp.py

`qp.py` is a Python script to process paired-end files using `pairs`
(for interleaving and uninterleaving files), `seqqs` (for gathering
read-quality information),
[scythe](https://github.com/vsbuffalo/scythe) (for adapter trimming),
and [seqtk](https://github.com/lh3/seqtk) (for quality trimming).

A sample run would look like this:

    python contrib/qp.py pe -a illumina_adapters.fa 269F6_CGATGT_L008_R1_001.fastq.gz \
	  269F6_CGATGT_L008_R2_001.fastq.gz

Here `269F6_CGATGT_L008_R1_001.fastq.gz` and
`269F6_CGATGT_L008_R2_001.fastq.gz` are the first and second read pair
FASTQ files. The subcommand `pe` indicates that you are processing a
paired-end file (single end support coming soon). The option `-a`
requires an argument: the path to a FASTA file of expected Illumina
adapters.

`qp.py` uses default parameters for Scythe's prior contamination reate
(0.4) and `trimfq` error rate (0.05). Scythe is robust againt
incorrect priors, with higher values usually working better. But, if
these need to be adjusted, they can be by setting `-p` (for Scythe's
prior) and `-e` for `trimfq`'s error rate.

Another option is `-s`, which will tell `qp.py` whether to use `pairs`
to split, or uninterleave the files into two read pair files and a
seperate file full of orphaned single reads (which can happen if a
read's pair was all poor quality and entirely trimmed). Final
splitting with `-s` should be used normally, but some programs do
require interleaved input.

One technical note: pair ordering is maintained because (1) if Scythe
finds an adapter at the first base, it trims everything and leaves a
single N and (2) `trimfq` will trim and leave up to 30bp by
default. So the order will not be changed by the tools in this
pipeline (it is also possible that they change in future versions
however, but this is unlikely). `seqqs` is run in a way in this script
(with its `-s` option so that it is "strict" and will break if reads
are not paired as expected).


