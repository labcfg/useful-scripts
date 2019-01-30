#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 18:04:35 2017

@author: Gaiaz Nugmanov
# thanks https://github.com/roippi/kmerclumps/blob/master/kmers.py
"""

import sys, os
import argparse
from tqdm import tqdm
import pandas as pd
from operator import itemgetter
import os.path
from collections import deque, Counter
from itertools import islice
from os.path import join
import distance
import pandas as pd
from Bio import SeqIO
#---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description = 'Calculate stats of kmer' + 
                                 'for sequences (*.fasta and *.txt(table))')
parser.add_argument('-ix', '--input', type = str, help = 'input file')
parser.add_argument('-colname', '--colname',
                    help = 'if *.txt(table): write name of column with read' + 
                    '(default = READ1_BEST); else: dont carry',
                    default = 'READ1_BEST')
parser.add_argument('-ox', '--output', type = str, help = 'output file')
parser.add_argument('-min','--minimum', type = int,
                    help = 'min length of kmer')
parser.add_argument('-max','--maximum', type = int,
                    help = 'max length of kmer')
parser.add_argument('-ham','--hamming', type = int,
                    help = 'max hamming distance for scc of kmer. default = 0',
                    default = 0)
if len(sys.argv)==1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

inputfile = os.path.abspath(args.input)
outputfile = args.output
k_min = args.minimum
k_max = args.maximum
ham = args.hamming
colname = args.colname

if k_min > k_max:
	sys.exit("\nError message: min > max. Please, change it\n")
if not os.path.isfile(inputfile):
	sys.exit('\nError message: Dont find inputfile or directory\n')

def sliding_window(seq, n=2):
    it = iter(seq)
    result = deque(islice(it, n), maxlen=n)
    if len(seq) <= n:
        yield seq
    for elem in it:
        result.append(elem)
        yield ''.join(result)
#---------------------------------------------------------------------------
def collapse_by_hamming(x, maxham):
    x_cllps = deque(sorted(Counter(x).items(), key=itemgetter(1),
                   reverse=True))
    while(len(x_cllps) > 0):
        seq, count = x_cllps.popleft()
        uniq = 1
        #x_cllps_save = deque()
        for aseq, acount in list(x_cllps):
            if distance.hamming(seq, aseq) <= maxham:
                count += acount
                uniq += 1
                x_cllps.remove((aseq, acount))
        yield ((seq, str(len(seq)), str(count), str(uniq)), len(x_cllps))
#---------------------------------------------------------------------------
def total_for_stream(x):
    x_cllps = deque(sorted(Counter(x).items(), key=itemgetter(1),
                    reverse=True))
    return len(x_cllps)
#---------------------------------------------------------------------------
def getkmers(genomelist, maxham, klen):
    kmer_holder = []
    for line in genomelist:
        for kmer in sliding_window(line.rstrip(), klen):
            kmer_holder.append(kmer)
    return (collapse_by_hamming(kmer_holder, maxham),
            total_for_stream(kmer_holder))
#---------------------------------------------------------------------------
def getseq(inputfile, colname):
    filename, ext = os.path.splitext(inputfile)
    if ext == '.txt':
        df = pd.read_table(inputfile)
        try:
            return df[colname]
        except:
            sys.exit('\nError message: Didnt find this column - ' +
                     colname + '. Check colname or extension\n')
    elif ext == '.fasta':
        fasta_sequences = SeqIO.parse(open(inputfile),'fasta')
        return [x.seq.tostring().upper() for x in fasta_sequences]
    else:
        sys.exit('\nError message: Didnt know this extension\n')
#---------------------------------------------------------------------------

if __name__ == "__main__":
    
    handle = open(outputfile, 'w')
    handle.write('\t'.join(['SEQ', 'LEN', 'AMOUNT', 'VAR']) + '\n')
    df = pd.read_table(inputfile)
    
    for klen in range(k_min, k_max+1):
        print('##### kmer_len = ' + str(klen) + ' #####')
        seq_list = getseq(inputfile, colname)
        kmer_stat, total = getkmers(seq_list, maxham=ham, klen=klen)
        tstream = tqdm(total=total, miniters=10)
        t_old = total
        for x, t in kmer_stat:
            tstream.update(t_old-t)
            handle.write('\t'.join(x) + '\n')
            t_old = t
        tstream.close()
    handle.close()