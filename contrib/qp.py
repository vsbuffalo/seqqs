#!/usr/bin/env python
import os
import re
import sys
import argparse
from subprocess import check_call

import pdb

# note: by default outputs everything in current directory

DEFAULT_SCYTHE_PRIOR = 0.4
DEFAULT_TRIMFQ_ERROR = 0.05

QUAL_PE_CMD = """\
pairs join -t -s {in_1} {in_2} \
| seqqs -k {k} -e -i -s -p {prefix}_raw - \
| scythe -a {adapters} -p {prior} - 2> {prefix}_scythe.stderr \
| seqtk trimfq -q {trim_error} - \
| seqqs -e -k {k} -i -s -p {prefix}_trimmed - """

FINAL_SPLIT_CMD = """\
| pairs split -1 {out_1} -2 {out_2} -u {out_unpaired} -
"""

def split_fastq_ext(filename):
    parts = re.match(r'(.*)\.(fastq(\.gz)?)$', filename, flags=re.I).groups()
    if sum(field is None for field in parts) > 1:
        raise ValueError("cannot parse FASTQ file name")
    return (parts[0], parts[1])


def make_outputfiles(infiles, outext="fastq"):
    outdict = dict()
    for i, filename in enumerate(infiles):
        filepath, ext = split_fastq_ext(os.path.basename(filename))
        outdict["out_%d" % (i+1)] = filepath + "-trimmed" + os.extsep + outext
    unpaired_file = re.sub(r"_R[12]$", "_RS-trimmed", filepath)
    outdict["out_unpaired"] = unpaired_file + os.extsep + outext
    outdict["prefix"] = re.sub(r"_R[12]$", "", filepath)
    return outdict

def process_pe(args):
    argsdict = dict(in_1=args.in1, in_2=args.in2, adapters=args.a, prior=args.p, trim_error=args.e)
    argsdict.update(make_outputfiles((args.in1, args.in2)))
    cmd = QUAL_PE_CMD
    if args.s:
       cmd += FINAL_SPLIT_CMD 
    cmd = cmd.format(**argsdict)
    sys.stderr.write("\nrunning pipeline:\n\t%s\n" % cmd)
    ret = check_call(cmd, shell=True)
    return ret
    
    
def process_se(args):
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='qp', description="""\
qp, run a generic quality processing pipeline on short reads. 

Note that seqqs (the program collecting quality information) does not
hash positional k-mers by default. If you turn this on, it can use a lot of
memory, so use wish caution. Short k-mers (around 6-8bp) are often
enough to detect contaminants like barcoes and adapters.""")
    parser.add_argument('-p', type=float, help='scythe prior contaminant rate', default=DEFAULT_SCYTHE_PRIOR)
    parser.add_argument('-e', type=float, help='trimfq error rate threshold', default=DEFAULT_TRIMFQ_ERROR)
    parser.add_argument('-k', type=int, help='k-mer size for positional k-mer hashing', 
                        default=0)
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_pe = subparsers.add_parser('pe', help='paired-end processing')
    parser_pe.add_argument('-a', type=str, help='adapters FASTA file', required=True)
    parser_pe.add_argument('-s', action="store_true", help="split interleaved file back into pairs and orphan file")
    parser_pe.add_argument('in1', type=str, help='read pair 1 FASTQ file')
    parser_pe.add_argument('in2', type=str, help='read pair 2 FASTQ file')
    parser_pe.set_defaults(func=process_pe)
    
    parser_se = subparsers.add_parser('se', help='single-end processing')
    parser_se.add_argument('-a', type=str, help='adapters FASTA file')
    parser_se.add_argument('in', type=str, help='singeles read FASTQ file')
    parser_se.set_defaults(func=process_se)
    
    args = parser.parse_args()
    args.func(args)
