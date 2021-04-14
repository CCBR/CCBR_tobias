import shutil

CONFIGFILE = str(workflow.overwrite_configfiles[0])

#
MEMORY="100"
# get working dir from config
WORKDIR = config["workdir"]

# get resources folder
try:
	RESOURCESDIR = config["resourcesdir"]
except KeyError:
	RESOURCESDIR = join(WORKDIR,"resources")

# get scripts folder
try:
	SCRIPTS = config["scripts"]
except KeyError:
	SCRIPTS = join(WORKDIR,"scripts")

## Load tools from YAML file
with open(config["tools"]) as f:
	TOOLS = yaml.safe_load(f)

GENOME = config["genome"]
REFFA = config["reffa"][GENOME]
AUTOCORRECT_EXTRA_PARAMS = config["autocorrect_extra_params"]
FOOTPRINTING_EXTRA_PARAMS = config["footprinting_extra_params"]
BINDETECT_EXTRA_PARAMS = config["bindetect_extra_params"]
BLACKLIST = config["blacklist"]
MOTIFS = config["motifs"]
if BLACKLIST != "":
	BLACKLIST = "--blackist "+BLACKLIST
# PEAKS = join(WORKDIR,"peaks.bed")
# shutil.copyfile(config["peaks"],PEAKS)
PEAKS = config["peaks"]

# print(config["data"])

#Check if there is at least one condition with bamfiles
if len(config["data"]) > 0:
	for condition in config["data"]:
		if len(config["data"][condition]) == 0:
			print("ERROR: Could not find any bamfiles in \"{0}\" in configfile {1}".format(":".join(("data", condition)), CONFIGFILE))
else:
	print("ERROR: Could not find any conditions (\"data:\{condition\}\") in configfile {0}".format(CONFIGFILE))
	sys.exit()

#Files related to experimental data (bam)
input_files = []
CONDITION_IDS = list(config["data"].keys())
for condition in CONDITION_IDS:
	if not isinstance(config["data"][condition], list):
		config['data'][condition] = [config['data'][condition]]

	cond_input = []
	for f in config['data'][condition]:
		globbed = glob.glob(f)
		if len(globbed) == 0:
			exit("ERROR: Could not find any files matching filename/pattern: {0}".format(f))
		else:
			cond_input.extend(globbed)

	config["data"][condition] = list(set(cond_input))						#remove duplicates
	input_files.extend(config['data'][condition])


contrasts=list(config["contrasts"].keys())
if len(contrasts) != len(set(contrasts)):
	exit("ERROR: Contrast names should be unique!")
CONTRASTS=list()
CONTRASTS2CONDITIONSBW=dict()
CONTRASTS2CONDITIONS=dict()
CC1=list()
CC2=list()
for contrast in contrasts:
	conditions=config["contrasts"][contrast]
	if not isinstance(conditions, list):
		exit("Error: contrast {0} should be a list".format(contrast))
	if len(conditions) != 2:
		exit("Error: contrast {0} needs 2 valid conditions".format(contrast))
	contrastname="_vs_".join(conditions)
	for cond in conditions:
		if not cond in CONDITION_IDS:
			exit("ERROR: Contrast {0} has a invalide condition {1}".format(contrast,cond))
		CC1.append(contrastname)
		CC2.append(cond)
	CONTRASTS.append(contrastname)
	CONTRASTS2CONDITIONSBW[contrastname]=list()
	CONTRASTS2CONDITIONS[contrastname]=conditions
	for c in conditions:
		fname=c+"_footprints.bw"
		CONTRASTS2CONDITIONSBW[contrastname].append(join(WORKDIR, "footprinting", fname))


# read TF list

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

TFs=get_TFs(MOTIFS)

output_files = []

# id2bam = {condition:{} for condition in CONDITION_IDS}
id2bam = dict()
for condition in CONDITION_IDS:
	config_bams = config['data'][condition]
	sampleids = [os.path.splitext(os.path.basename(bam))[0] for bam in config_bams]
	id2bam[condition] = {sampleids[i]:config_bams[i] for i in range(len(sampleids))}	# Link sample ids to bams

# print(id2bam)
output_files.extend(expand(os.path.join(WORKDIR, "coverage", "{condition}_coverage.bw"), condition=CONDITION_IDS))
output_files.extend(expand(os.path.join(WORKDIR, "footprinting", "{condition}_footprints.bw"), condition=CONDITION_IDS))
output_files.extend(expand(os.path.join(WORKDIR, "overview", "all_{condition}_bound.bed"), condition=CONDITION_IDS))
# print(output_files)
# sys.exit()