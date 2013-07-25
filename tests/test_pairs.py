# Test that read pair files being processed by pairs are not being
# broken in any way. The best way to ensure this is to join and split
# and ensure the files are exactly the same.

import sys
import os
import pdb
from subprocess import call, check_output, CalledProcessError

devnull = open(os.devnull, 'w')

# some presets
bad_files = ("bad-1.fq", "bad-2.fq") # breaks strict, unmatching headers
ok_files = ("ok-1.fq", "ok-2.fq") # differing only in /1 and /2

def make_outfile(name):
    name, ext = os.path.splitext(name)
    return name + "-out" + "." + ext 

def test_join_split(infiles):
    outfiles = [make_outfile(f) for f in infiles]
    outfiles.append("unpaired.fq")
    args = tuple(list(infiles) + outfiles)
    results = list()
    cmd = "../pairs join -s %s %s | ../pairs split -1 %s -2 %s -u %s - " % args
    print "running:", cmd
    results.append(call(cmd, shell=True) == 0)
    md5s = dict()
    for f in args:
        cmd = "md5 %s" % f
        print "running:", cmd
        md5s[f] = check_output(cmd, shell=True).split(" = ")[1]
    for outfile in outfiles:
        results.append(md5s[outfile.replace("-out", "")] == md5s[outfile])
    return all(results)

def test_join_strict():
    """
    Test that strict fails with unmatching names, but is ok in all
    other cases.
    """
    results = list()
    cmd1 = "../pairs join -s %s %s" % bad_files
    print "running:", cmd1
    results.append(call(cmd1, shell=True, stdout=devnull, stderr=devnull) == 1)
    cmd2 = "../pairs join %s %s" % bad_files
    print "running:", cmd2
    results.append(call(cmd2, shell=True, stdout=devnull, stderr=devnull) == 0)
    cmd3 = "../pairs join -s %s %s" % ok_files
    results.append(call(cmd3, shell=True, stdout=devnull, stderr=devnull) == 0)
    if not all(results):
        pdb.set_trace()
    return all(results)

if __name__ == "__main__":
    if (len(sys.argv) < 3):
        sys.exit("usage: test_pairs.py pair-1.fq pair-2.fq\n")
    infiles = (sys.argv[1], sys.argv[2])
    tests = list()
    tests.append(("tests_join_split", test_join_split(infiles)))
    tests.append(("test_join_strict", test_join_strict()))
    total = 0
    not_passed = 0
    print "results:"
    for name, value in tests:
        total += 1
        not_passed += int(value)
        print "\t%s\t%s" % (name, ["Failed", "Passed"][int(value)])
    if not all(tests):
        sys.exit("%d/%s tests failed!\n" % (not_passed, total))
    sys.stderr.write("%d/%s tests passed.\n" % (not_passed, total))
