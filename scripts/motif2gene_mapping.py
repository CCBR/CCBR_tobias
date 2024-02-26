import sys
if len(sys.argv) != 3:
    print("2 arguments required!")
    print("1: GTF file")
    print("2: Jaspar motifs file")
    exit()
protein2gene=dict()
with open(sys.argv[1],'r') as infile:
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
with open(sys.argv[2],'r') as infile:
    for l in infile:
        if not l.startswith(">"): continue
        motif = l.strip().split()[-1]
        umotif = "\"" + motif.upper() + "\""
        if not umotif in protein2gene: continue
        print("\t".join(list(map(lambda x:x.strip("\""),[motif,protein2gene[umotif]]))))
