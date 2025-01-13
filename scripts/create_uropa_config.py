import json
import argparse
parser = argparse.ArgumentParser(description='create uropa input json file')
parser.add_argument('-j',dest="json",required=True,help="base template JSON absolute path")
parser.add_argument('-g',dest="gtf",required=True,help="GTF file absolute path")
parser.add_argument('-b',dest="bed",required=True,help="peaks BED file absolute path")
args = parser.parse_args()
config = json.load(open(args.json))
config["gtf"] = args.gtf
config["bed"] = args.bed

string_config=json.dumps(config, indent=4)
print(string_config)
