import json
config = dict()
config = {"queries":[
			{"feature":"gene", "feature.anchor":"start", "distance":[10000,1000], "filter_attribute":"gene_biotype", "attribute_values":"protein_coding", "name":"protein_coding_promoter"},
			{"feature":"gene", "distance":1, "filter_attribute":"gene_biotype", "attribute_values":"protein_coding", "internals":0.1, "name":"protein_coding_internal"},
			{"feature":"gene", "feature.anchor":"start", "distance":[10000,1000], "name":"any_promoter"},
			{"feature":"gene", "distance":1, "internals":0.1, "name":"any_internal"},
			{"feature":"gene", "distance":[50000, 50000], "name":"distal_enhancer"},
			],
		"show_attributes":["gene_biotype", "gene_id", "gene_name"],
		"priority":"True"
		}

config["gtf"] = "/data/CCBR_Pipeliner/db/PipeDB/Indices/GTFs/mm10/gencode.vM21.annotation.gtf"
config["bed"] = "/data/CCBR/projects/ccbr872/atacseq/test/peaks/genrich/tn5knicks/query_regions.bed"

string_config = json.dumps(config, indent=4)
config_file = open("uropa.config", "w")
config_file.write(string_config)
config_file.close()