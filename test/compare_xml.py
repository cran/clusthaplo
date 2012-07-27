#!/usr/bin/env python

import os
import sys

my_dir = sys.argv[1]
his_dir = sys.argv[2]
my = set([x for x in os.listdir(my_dir) if x.endswith('.xml')])
his = set([x for x in os.listdir(his_dir) if x.endswith('.xml')])
our = sorted(my.intersection(his))


def read_xml(f):
    lines = open(f).readlines()
    probas = [l for l in lines if '<PROBABILITY' in l]
    rest = [l for l in lines if '<PROBABILITY' not in l]
    return (probas, rest)


def compare_xml(my, his):
    my_p, my_r = my
    his_p, his_r = his
    return (len(his_r) == len(my_r),
            filter(lambda x: x[0] != x[1], zip(his_r, my_r)),
            len(his_p) == len(my_p),
            '%3f%%' % (reduce(lambda a, (my, his): my == his and 1 + a or a,
                              zip(my_p, his_p), 0) * 100. / len(my_p)),
            len([x for x in his_p if '1' in x]),
            len([x for x in my_p if '1' in x]),
            len(my_p))


for xml in our:
    print xml,
    c = compare_xml(read_xml(os.path.join(my_dir, xml)),
                    read_xml(os.path.join(his_dir, xml)))
    if not c[0]:
        print "[number of lines differ]",
    if not c[2]:
        print "[number of PROBAs differ]",
    print c[3]
    if c[1]:
        print '\n'.join(' <-> '.join(map(repr, x)) for x in c[1])
    if c[3] != '100%':
        print "Nproba(1)", c[4], c[5], 'out of', c[6]
