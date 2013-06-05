# seqqs

Seqqs (SEQuence Quality Statistics, pronounced "seeks") is a C library
for quickly gathering quality statistics from sequence files. It's
mostly adapted from [qrqc](http://github.com/vsbuffalo/qrqc), except
it is designed to be run in quality processing pipelines. It can also
be compiled as a dynamic library and called from other programs.

## Requirements and Installation

Seqqs can be compiled using GCC or Clang; compilation during
development used the latter. Seqqs relies on Heng Li's kseq.h and
khash.h, which is bundled with the source.

Seqqs requires Zlib, which can be obtained at <http://www.zlib.net/>.

To install, just run `make` in the `seqqs` directory.

## Usage

Documentation is internal; just compile and run `./seqqs`. Here are
some usage examples.

Without any options, `seqqs` works like so:

    cat in.fq | seqqs -
	# or:
	seqqs in.fq
	
Note that `-` tells `seqqs` to read from standard input. Without any
options, this will create `qual.txt`, `nucl.txt`, and `len.txt`.

`seqqs` is designed to be placed in pipelines and act as a quality
gathering step without disrupting the flow (similar to Unix `tee`). To
enable this, use `-e` (for emit):

    cat in.fq | seqqs -e -
	
For complex quality pipelines, `seqqs` can also take a prefix argument
to prevent overwriting output files. If we wanted to create a complex
workflow that gathers quality on raw input, gathers quality
statistics, then trims using Heng Li's
[seqtk](http://github.com/lh3/seqtk) trimfq command, and then gathers
output statistics, we could use: 

    cat in.fq | seqqs -e -p raw-$(date +%F) - | seqtk trimfq - | \
	  seqqs -e -p trimmed-$(date +%F) > trimmed.fq

`seqqs` can also gather positional k-mers, which can help in
discovering enrichment due to positional contaminants like untrimmed
barcodes and adapters. As a quick aside: you should check for these!
Many sequencing data set are plagued by positional contaminants,
especially as barcoding grows in popularity. The k-mer option is `-k
<n>` where *n* is the k-mer size:

    cat in.fq | seqqs -k 6
	
`seqqs` can also work with **interleaved** paired-end files. The
results are no different, but two output files (one for each set of
reads in a pair) are created. These have the names like the default,
except they have `_1.txt` and `_2.txt` suffixes. Also, `seqqs` will
warn if pairing looks incorrect. If `-s` (strict) is set, `seqqs` will
error out if interleaved pairs do not have the same name (ignoring
`/1` and `/2` and excluding the comment).

## Using Output

qrqc will soon have functions to gather this output and make plots
from it. For now, 
