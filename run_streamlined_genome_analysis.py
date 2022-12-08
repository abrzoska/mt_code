import cut_maf
from file_parameters import *
from variables import *
import sys
import count_indels
import gene_mapper
sys.path.append("/home/tlin/mt_code/brex/")

try:
    import pandas as pd
except ImportError as e:
    sys.exit("Module pandas not found.")
try:
    import pyranges as pr
except ImportError as e:
    sys.exit("Module pyranges not found.")

from brex import region_extractor as region_extractor

try:
    indel_length = int(sys.argv[1])
except:
    raise "Error: argument for indel size needs to be integer"

run_name = f"indels{indel_length}"
run_variables = RunVariables(run_name)

if needs_preprocessing:
    #cut_maf.main()
    region_extractor.loop_find_indels_for_query_from_maf(number_of_maf_parts, input_maf_part, target_species,
                                                         query_species, indel_file, adapted_dict, number_of_cores)
    gene_mapper.map_genes_to_cis_regs()

count_indels.run_parrallel_analysis(run_variables.cis_reg_indel_folder, indel_file, indel_length, number_of_maf_parts, cis_reg_file)
genes_to_cis_reg = count_indels.determine_cis_regs_genes(gene_input_file, gene_to_cis_reg_dic)
count_indels.summarize_cis_regs(genes_to_cis_reg, run_variables.cis_reg_indel_folder, run_variables.merged_df_pickle, run_variables.cis_reg_pickle,
                                run_variables.result_pickle).to_csv(run_variables.gene_result_csv)
count_indels.summarize_genes(run_variables.gene_result_csv, run_variables.gene_result_summary_txt)
