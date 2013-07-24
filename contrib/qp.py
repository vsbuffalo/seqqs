import sys
import argparse
from subprocess import check_call

QUAL_CMD = """
pairs join -t -s {in_1} {in_2} \
| seqqs -e -i -s -p raw - \
| scythe -a {adapters} -p {prior} - \
| seqtk trimfq -q {trim_error} - \
| seqqs -e -i -s -p trimmed - \
| pairs split -1 {out_1} -2 {out_2} -u {out_unpaired} -
"""

def process_pe(args):
    argdict = {in_1=args.in1, in_2=args.in2, adapters=args.a, prior=args.p, trim_error=args.e}
    cmd = QUAL_CMD.format(**argsdict)
    sys.stderr.write("running pipeline\n\t%s\n" % cmd)
    ret = system.check_call(cmd)
    
    
def process_se(args):
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='qp')
    parser.add_argument('-p', type=float, help='scythe prior contaminant rate', default=DEFAULT_SCYTHE_PRIOR)
    parser.add_argument('-e', type=float, help='trimfq error rate threshold', default=DEFAULT_TRIMFQ_ERROR)
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_pe = subparsers.add_parser('pe', help='paired-end processing')
    parser_pe.add_argument('-a', type=str, help='adapters FASTA file', required=True)
    parser_pe.add_argument('in1', type=str, help='read pair 1 FASTQ file', required=True)
    parser_pe.add_argument('in2', type=str, help='read pair 2 FASTQ file', required=True)

    parser_se = subparsers.add_parser('se', help='single-end processing')
    parser_se.add_argument('-a', type=str, help='adapters FASTA file', required=True)
    parser_se.add_argument('in', type=str, help='singeles read FASTQ file', required=True)

    
