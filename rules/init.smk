#######################################################
import json
import shutil
#######################################################

print("#"*100)
print("""
#   ____ ____ ____  ____  
#  / ___/ ___| __ )|  _ \ 
# | |  | |   |  _ \| |_) |
# | |__| |___| |_) |  _ < 
#  \____\____|____/|_| \_\ 
#  _____ ___  ____ ___    _    ____  
# |_   _/ _ \| __ )_ _|  / \  / ___| 
#   | || | | |  _ \| |  / _ \ \___ \ 
#   | || |_| | |_) | | / ___ \ ___) |
#   |_| \___/|____/___/_/   \_\____/ 
#  ____  _            _ _            
# |  _ \(_)_ __   ___| (_)_ __   ___ 
# | |_) | | '_ \ / _ \ | | '_ \ / _ \ 
# |  __/| | |_) |  __/ | | | | |  __/
# |_|   |_| .__/ \___|_|_|_| |_|\___|
#         |_|                        
""")
print("#"*100)


#########################################################
# FILE-ACTION FUNCTIONS 
#########################################################
def check_existence(filename):
  if not os.path.exists(filename):
    exit("# File/Folder: %s does not exists!"%(filename))

def check_readaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.R_OK):
    exit("# File: %s exists, but cannot be read!"%(filename))

def check_writeaccess(filename):
  check_existence(filename)
  if not os.access(filename,os.W_OK):
    exit("# File: %s exists, but cannot be read!"%(filename))

def get_file_size(filename):
    filename=filename.strip()
    if check_readaccess(filename):
        return os.stat(filename).st_size

#########################################################
# MOTIF related FUNCTIONS 
#########################################################

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
    # print(file_format)
    if file_format != "jaspar" and file_format != "pfm":
        print("ERROR: cannot read {0} file".format(path))
        exit("ERROR: motifs need to be in PFM (preferred) format!")

    f=open(path).readlines()
    TFs_tmp=list(map(lambda x:x.strip().replace(">",""),filter(lambda x:x.startswith(">"),f)))
    if file_format == "pfm" and len(TFs_tmp[0].split()) != 1:
        exit("ERROR: motifs are not in valid PFM (preferred) format!")
    if file_format == "pfm":
        TFs=list(map(lambda x:x+"_"+x,TFs_tmp))
    return TFs

#########################################################
# DEFINE CONFIG FILE AND READ IT
#########################################################
CONFIGFILE = str(workflow.overwrite_configfiles[0])

# set memory limit 
# used for sambamba sort, etc
# MEMORY="100"

# read in various dirs from config file
WORKDIR = config["workdir"]
# get resources folder
try:
    RESOURCESDIR = config["resourcesdir"]
except KeyError:
    RESOURCESDIR = join(WORKDIR,"resources")
check_existence(RESOURCESDIR)

# get scripts folder
try:
    SCRIPTSDIR = config["scriptsdir"]
except KeyError:
    SCRIPTSDIR = join(WORKDIR,"scripts")
check_existence(SCRIPTSDIR)

# genome fasta
GENOME = config["genome"]
REFFA = config[GENOME]["reffa"]
check_readaccess(REFFA)
GTF = config[GENOME]["gtf"]
check_readaccess(GTF)

# extra Tobias parameters
AUTOCORRECT_EXTRA_PARAMS = config["autocorrect_extra_params"]
FOOTPRINTING_EXTRA_PARAMS = config["footprinting_extra_params"]
BINDETECT_EXTRA_PARAMS = config["bindetect_extra_params"]

# containers
CONTAINERS = config["containers"]

# blacklist BED, if any
BLACKLIST = config["blacklist"]
if BLACKLIST != "":
    check_readaccess(BLACKLIST)
    BLACKLIST = "--blackist "+BLACKLIST

# motifs to use for testing
MOTIFS = config["motifs"]
check_readaccess(MOTIFS)

# filtered peaks
PEAKS = config["peaks"]
# PEAKS = join(WORKDIR,"peaks.bed")
# shutil.copyfile(config["peaks"],PEAKS)
# print(config["data"])

#########################################################
# READ IN TOOLS REQUIRED BY PIPELINE
# THESE INCLUDE LIST OF BIOWULF MODULES (AND THEIR VERSIONS)
# MAY BE EMPTY IF ALL TOOLS ARE DOCKERIZED
#########################################################

## Load tools from YAML file
try:
    TOOLSYAML = config["tools"]
except KeyError:
    TOOLSYAML = join(WORKDIR,"tools.yaml")
check_readaccess(TOOLSYAML)
with open(TOOLSYAML) as f:
    TOOLS = yaml.safe_load(f)

#########################################################
# UROPA BASE CONFIG TEMPLATE
#########################################################

try:
    UROPABASEYAML = config["uropabaseyaml"]
except KeyError:
    UROPABASEYAML = join(WORKDIR,"uropa_base_config.yaml")
check_readaccess(UROPABASEYAML)

#########################################################
# READ CLUSTER PER-RULE REQUIREMENTS
#########################################################

## Load cluster.json
try:
    CLUSTERJSON = config["clusterjson"]
except KeyError:
    CLUSTERJSON = join(WORKDIR,"cluster.json")
check_readaccess(CLUSTERJSON)
with open(CLUSTERJSON) as json_file:
    CLUSTER = json.load(json_file)
## Create lambda functions to allow a way to insert read-in values
## as rule directives
getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER and "threads" in CLUSTER[rname] else int(CLUSTER["__default__"]["threads"])
getmemg=lambda rname:CLUSTER[rname]["mem"] if rname in CLUSTER else CLUSTER["__default__"]["mem"]
getmemG=lambda rname:getmemg(rname).replace("g","G")
#########################################################


#########################################################
# MAKE A LIST OF CONDITIONS
#########################################################

#Check if there is at least one condition with bamfiles
if len(config["data"]) > 0:
    for condition in config["data"]:
        if len(config["data"][condition]) == 0:
            print("ERROR: Could not find any bamfiles in \"{0}\" in configfile {1}".format(":".join(("data", condition)), CONFIGFILE))
else:
    print("ERROR: Could not find any conditions (\"data:\{condition\}\") in configfile {0}".format(CONFIGFILE))
    sys.exit()

# Files related to experimental data (dedup bam)

input_files = [] # list of all input bam files
CONDITION_IDS = list(config["data"].keys()) # list of all conditions
CONDITION2BAMS = dict()
# force each condition to be a list
for condition in CONDITION_IDS:
    if not condition in CONDITION2BAMS:
        CONDITION2BAMS[condition] = list()

    if not isinstance(config["data"][condition], list):
        config['data'][condition] = [config['data'][condition]]

    cond_input = []
    for f in config['data'][condition]:
        globbed = glob.glob(f)
        if len(globbed) == 0:
            exit("ERROR: Could not find any files matching filename/pattern: {0}".format(f))
        else:
            cond_input.extend(globbed)

    config["data"][condition] = list(set(cond_input))					#remove duplicates
    CONDITION2BAMS[condition] = config["data"][condition]
    input_files.extend(config['data'][condition])

#########################################################
# MAKE A LIST OF CONTRASTS
#########################################################

# contrasts related global variables
contrasts=list(config["contrasts"].keys())
if len(contrasts) != len(set(contrasts)):
    exit("ERROR: Contrast names should be unique!")
CONTRASTS=list() # all contrasts to be run
CONTRASTS2CONDITIONSBW=dict() # contrast name to condition bigwigs in that contrast
CONTRASTS2CONDITIONS=dict() # contrast name to conditions in that contrast
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
            print("ERROR: Contrast {0} has a invalid condition {1}".format(contrast,cond))
            exit("ERROR: Known conditions are: {0}".format(str(CONDITION_IDS)))
        CC1.append(contrastname)
        CC2.append(cond)
    CONTRASTS.append(contrastname)
    CONTRASTS2CONDITIONSBW[contrastname]=list()
    CONTRASTS2CONDITIONS[contrastname]=conditions 
    for c in conditions:
        fname=c+"_footprints.bw"
        CONTRASTS2CONDITIONSBW[contrastname].append(join(WORKDIR, "footprinting", fname))

#########################################################
# LOAD IN THE MOTIFS FROM THE TF DATABASE
#########################################################
# read in the motifs
TFs=get_TFs(MOTIFS)
output_files = []

# id2bam = {condition:{} for condition in CONDITION_IDS}
id2bam = dict()
for condition in CONDITION_IDS:
    config_bams = config['data'][condition]
    sampleids = [os.path.splitext(os.path.basename(bam))[0] for bam in config_bams]
    id2bam[condition] = {sampleids[i]:config_bams[i] for i in range(len(sampleids))}	# Link sample ids to replicate bams

# print(id2bam)
output_files.extend(expand(os.path.join(WORKDIR, "coverage", "{condition}_coverage.bw"), condition=CONDITION_IDS))
output_files.extend(expand(os.path.join(WORKDIR, "footprinting", "{condition}_footprints.bw"), condition=CONDITION_IDS))
output_files.extend(expand(os.path.join(WORKDIR, "overview", "all_{condition}_bound.bed"), condition=CONDITION_IDS))
# print(output_files)
# sys.exit()

#########################################################
# PRINT PIPELINE PARAMETERS
#########################################################

print("# Pipeline Parameters:")
print("#"*100)
print("# Working dir :",WORKDIR)
# print("# Results dir :",RESULTSDIR)
print("# Scripts dir :",SCRIPTSDIR)
print("# Resources dir :",RESOURCESDIR)
print("# Tools YAML: ", TOOLSYAML)
print("# Uropa Base Template YAML: ",UROPABASEYAML)
print("# Cluster JSON: ",CLUSTERJSON)
print("# Genome :",GENOME)
print("# Reference Fasta :",REFFA)
if BLACKLIST != "": 
    print("# Blacklist Bed :",BLACKLIST)
print("# Motifs database: ",MOTIFS)
print("# Peaks Bed :",PEAKS)
print("# Conditions: ")
for i,c in enumerate(CONDITION_IDS):
    print("# \t {0}. {1} ".format(i+1,c))
print("# Constrasts: ")
for i,c in enumerate(CONTRASTS2CONDITIONS.keys()):
    print("# \t {0}. {1}:{2} ".format(i+1,c,CONTRASTS2CONDITIONS[c]))
print("#"*100)

# exit("Done")

