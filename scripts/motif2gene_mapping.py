#!/usr/bin/env python
"""
Map Jaspar motifs to genes
"""

import sys


def main(gtf_file, motifs_file, output_file):
    protein2gene=dict()
    with open(gtf_file,'r') as infile:
        for l in infile:
            l = l.strip().split()
            if len(l)<3 or l[2] != "gene": continue
            gid = ""
            gname = ""
            for i,x in enumerate(l):
            #    print("##"+x+"##",i)
                if x=="gene_id": gid=i+1
                if x=="gene_name": gname=i+1
            if gid != "" and gname != "":
                x=l[gid].strip(";").split(".")[0]+"\""
                x=x.upper()
                y=l[gname].strip(";")
                if not y in protein2gene: protein2gene[y]=x
    with open(output_file, 'w') as outfile:
        with open(motifs_file,'r') as infile:
            for l in infile:
                if not l.startswith(">"): continue
                motif = l.strip().split()[-1]
                umotif = "\"" + motif.upper() + "\""
                if not umotif in protein2gene: continue
                outfile.write("\t".join(list(map(lambda x:x.strip("\""),[motif,protein2gene[umotif]]))))


if __name__ == '__main__':
    if "snakemake" in locals() or "snakemake" in globals():
        main(snakemake.input.gtf, snakemake.input.pfm, snakemake.output.txt)
    else:
        if len(sys.argv) < 4:
            raise ValueError("\n".join(["3 arguments required!", "1: GTF file", "2: Jaspar motifs file", "3: Output file"]))
        else:
            main(*sys.argv[1:])
