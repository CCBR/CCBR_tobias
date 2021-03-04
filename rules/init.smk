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
PEAKS = join(WORKDIR,"peaks.bed")
shutil.copyfile(config["peaks"],PEAKS)


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