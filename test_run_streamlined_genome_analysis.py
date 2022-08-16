from test_file_parameters import *
from variables import *
import sys
sys.path.append("/home/tlin/mt_code/brex/")
import map_genes_to_cisregs
import mapping_genes

try:
    import pandas as pd
except ImportError as e:
    sys.exit("Module pandas not found.")
try:
    import pyranges as pr
except ImportError as e:
    sys.exit("Module pyranges not found.")
try:
    from brex import region_extractor as rex
    from brex import helper as h
except ImportError as e:
    print(e)
    sys.exit("Module brex not found. Try running setup.py")

run_name = "whole_maf_indel_"
target_species = "mm10"
output_folder = intermediary_folder

try:
    indel_length = int(sys.argv[1])
except:
    raise "Error: argument for indel size needs to be integer"

#Todo: change this indel thing
reg_ex = rex.RegionExtractor(run_name, intermediary_folder, indel_length)
reg_ex.find_indels_from_maf(input_maf, target_species, query_species, indel_file_name, adapted_dict)
#print("finished step 1")
#mapping_genes.map_genes_to_cis_regs()
map_genes_to_cisregs.run_analysis(cis_reg_indel_folder, indel_file_name)
genes_to_cis_reg = map_genes_to_cisregs.determine_cis_regs_genes(gene_input_file)
map_genes_to_cisregs.summarize_cis_regs(genes_to_cis_reg, cis_reg_indel_folder).to_csv(gene_result_csv)
#print("finished step 4")
