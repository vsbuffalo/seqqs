# Wrappers Scripts

`seqqs` is designed to work with a stream, through pipes or
[process substitution](http://vincebuffalo.org/2013/08/08/the-mighty-named-pipe.html),
and fits well in with other tools I've written or contributed to like
[Scythe](https://github.com/vsbuffalo/scythe) and
[Sickle](https://github.com/najoshi/sickle). Scythe does adapter
trimming (when supplied a FASTA file of adpaters) and Sickle does
paired end quality trimming.

To make quality processing easier, this directory contains `trim.sh`,
a simple Bash-based wrappper script for an entire quality
pipeline. This script takes three command line arguments: the sample
name, the paired-end file 1, and the paired-end file 2. There are
three additional optional positional arugments: the adapters file, the
prior to use with Scythe, and the quality threshold to use with
Sickle. Other specific options can be adjusted in the script.

All options can be seen by running `trim.sh` without arugments. Note
that this assumes Sanger-quality encoded paired-end Illumina sequences
from the CASVAVA pipeline version 1.8 or greater. Also, at this point,
`trim.sh` won't merge many split Illumina files, which will have names
like:

    10_116_CACGAT_L001_R1_001.fastq.gz
	10_116_CACGAT_L001_R1_002.fastq.gz
	# etc
	
These will have to merge yourself, or adapt the script (a version that
does this automatically may be added).

An example run would look like:

    bash trim.sh sample3 sample3_R1_001.fastq.gz sample3_R2_001.fastq.gz
	
This creates `[raw,trimmed]_sample3_[R1,R2]_[qual,nucl,len].txt` which
are the quality statistics output files from `seqqs`. The prefix `raw`
and `trimmed` correspond to before and after trimming, and these can
be quite useful to compare. Additionally, `trim.sh` will create
`sample3_R1_trimmed.fastq.gz` and `sample3_R2_trimmed.fastq.gz`, which
are the trimmed files. Scythe and Sickle's diagnostics, etc, can be
seen in the `*.stderr` file.
