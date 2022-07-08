from file_parameters import *
from variables import *
import sys
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

reg_ex = rex.RegionExtractor(run_name, intermediary_folder, indel_length)
reg_ex.loop_find_indels_for_query_from_maf(number_of_maf_parts, input_maf_part, target_species, query_species, in_group, indel_file, adapted_dict, number_of_cores)
print("finished step 1")
#mapping_genes.map_genes_to_cis_regs()
map_genes_to_cisregs.main()
print("finished step 4")
