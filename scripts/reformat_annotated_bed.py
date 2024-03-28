#!/usr/bin/env python
"""
Map Jaspar motifs to genes
"""

import sys


def main(input_filename, output_filename):
    gene_name_index = 22
    motif_name_index = 3
    with open(output_filename, 'w') as outfile:
        with open(input_filename, 'r') as infile:
            for line in infile:
                line_split = line.split('\t')
                line_split[motif_name_index] = line_split[motif_name_index].split('_')[0]
                line_split[gene_name_index] = line_split[gene_name_index].split('.')[0]
                line_joined = '\t'.join(line_split)
                outfile.write(line_joined + '\n')

if __name__ == '__main__':
    if len(sys.argv) < 3:
        raise ValueError("\n".join(["2 arguments required!", "1: input file", "2: output file"]))
    else:
        main(*sys.argv[1:])
