#!/usr/bin/env python3

# Fiona Tamburini
# Sept 2018
# parse cd-hit output and gff file into quasi-gff format for plotting in R

# usage
# ./parseClusters.py cdhit/bmt_crass_subset-Rotated.faa.clstr prodigal/bmt_crass_subset-Rotated.gff

# to do -- check how forward/reverse strand is handeled

import sys
import re

cluster_file = snakemake.input[0]
gff_file = snakemake.input[1]
out_file = snakemake.output[0]

sys.stdout = open(out_file,'wt')

# read cluster file into dictionary
clusters = {}

with open(cluster_file, 'r') as c:
    for line in c:
        if line.startswith('>'):
            curr_cluster = line.rstrip().lstrip('>')
        else:
            info = line.rstrip().split('\t')[1].split(' ')[1]
            node, start, end = info.split('_#_')[0:3]
            node = node.lstrip('>')
            node = re.sub(r'_[0-9]+$', '', node)

            if node in clusters:
                clusters[node][(start, end)] = curr_cluster
            else:
                clusters[node] = {(start, end): curr_cluster}

# parse gff file and add cluster information
with open(gff_file, 'r') as gff:
    for line in gff:

        # skip if info headers
        if not line.startswith('#'):
            seqid, source, type, start, end, score, strand, phase, attributes = line.rstrip().split('\t')

            print('\t'.join([seqid, source, type, start, end, score, strand, phase, clusters[seqid][(start, end)]]))

            # if forward strand
            # if start < end:
            #     print('\t'.join([seqid, source, type, start, end, score, "+", phase, clusters[seqid][(start, end)]]))
            # else:
            #     print('\t'.join([seqid, source, type, start, end, score, '-', phase, clusters[seqid][(start, end)]]))
