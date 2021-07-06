# for Bio import motifs
import re

def get_motif_format(content):
	""" Get motif format from string of content """
	
	#Estimate input format
	if re.match(r".*MEME version.+", content, re.DOTALL) is not None: # MOTIF\s.+letter-probability matrix.+[\d\.\s]+", content, re.MULTILINE) is not None:
		motif_format = "meme"

	elif re.match(r">.+A.+\[", content, re.DOTALL) is not None:
		motif_format = "jaspar"

	elif re.match(r">.+", content, re.DOTALL) is not None:
		motif_format = "pfm"

	elif re.match(r"AC\s.+", content, re.DOTALL) is not None:
		motif_format = "transfac"
	
	else:
		motif_format = "unknown"

	return(motif_format)

def get_TFs(path):

	content = open(path).read()
	file_format = get_motif_format(content)
	print(file_format)
	if file_format != "jaspar" and file_format != "pfm":
		exit("ERROR: motifs need to be in PFM (preferred) format!")

	f=open(path).readlines()
	TFs_tmp=list(map(lambda x:x.strip().replace(">",""),filter(lambda x:x.startswith(">"),f)))
	if file_format == "pfm" and len(TFs_tmp[0].split()) != 1:
		exit("ERROR: motifs are not in valid PFM (preferred) format!")
	if file_format == "pfm":
		TFs=list(map(lambda x:x+"_"+x,TFs_tmp))
	return TFs

path="/data/kopardevn/GitRepos/CCBR_tobias/resources/HOCOMOCOv11_HUMAN_core_motifs.txt"
TFs=get_TFs(path)
print(TFs)
