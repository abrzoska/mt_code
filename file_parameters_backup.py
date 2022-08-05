import os

run_name = "wrn_test"

#TODO: make run directories

def create_dic_if_non_existent(path):
    exists = os.path.exists(path)
    if not exists:
        os.makedirs(path)

create_dic_if_non_existent(f"/mnt/moreAddspace/mt_code_data/intermediary_data/{run_name}")
create_dic_if_non_existent(f"/mnt/moreAddspace/mt_code_data/debug/{run_name}")
create_dic_if_non_existent(f"/home/tlin/mt_code_data/intermediary_data/{run_name}/indels/")
create_dic_if_non_existent(f"/home/tlin/mt_code_data/intermediary_data/{run_name}/cis_reg_indels/")

#inputs
gene_input_file = "/mnt/moreAddspace/mt_code_data/data/Cross_all_final_200421.xlsx"

#data
gene_promoter_file = "/mnt/moreAddspace//mnt/moreAddspace/mt_code_data/data/gene_promoters.tsv"
cis_reg_file = "/mnt/moreAddspace/mt_code_data/data/encode9_trimmed.bed"
gene_file = "/mnt/moreAddspace/mt_code_data/data/mart_export2.txt"
input_maf = "/mnt/moreAddspace/abrzoska/mm10/multi.anno.maf"
input_maf_part = "/home/tlin/maf_split/maf_part_NUMBER.maf"
maf_split = "/home/tlin/maf_split/"

#intermdiary
indel_file = f"/home/tlin/mt_code_data/intermediary_data/{run_name}/indels/maf_part_NUMBER.bed"
indel_file_name = f"/home/tlin/mt_code_data/intermediary_data/{run_name}/all_chromosomes_0.bed"
#bed_out = f"/mnt/moreAddspace/mt_code_data/intermediary_data/{run_name}/indels/all_chromosomes"
cis_reg_indel_folder = f"/home/tlin/mt_code_data/intermediary_data/{run_name}/cis_reg_indels/"
intermediary_folder = f"/home/tlin/mt_code_data/intermediary_data/{run_name}/"

#result
gene_result_csv = f"/mnt/moreAddspace/mt_code_data/mt_code_results/{run_name}_genes_csv"

#pickles
merged_df_pickle = f"/mnt/moreAddspace/mt_code_data/debug/{run_name}/merged_dfs.pkl"
cis_reg_pickle = f"/mnt/moreAddspace/mt_code_data/debug/{run_name}/genes_to_cis_reg.pkl"
result_pickle = f"/mnt/moreAddspace/mt_code_data/debug/{run_name}/result.pkl"
gene_to_cis_reg_dic = f"/mnt/moreAddspace/mt_code_data/debug/{run_name}/gene_to_cis_reg.pkl"

#scripts
splitter_script = "/home/timlin/projects/mt_code/chromosome_splitter.sh"

#test
motif_file = "/mnt/moreAddspace/mt_code_data/test_data/wrn_fimo.tsv"
wrn_cis_reg_motifs = "/mnt/moreAddspace/mt_code_data/test_data/wrn_result.tsv"

input_maf = "/mnt/moreAddspace/mt_code_data/test_data/ENSMUSG00000025475-indels.maf"
maf_split = "/mnt/moreAddspace/mt_code_data/test_data/maf_split/"
input_maf_part = "/mnt/moreAddspace/mt_code_data/test_data/maf_split/maf_part_NUMBER.maf"


#meta #TODO: computationally determine this and include the cutter script in the pipeline
#number_of_maf_parts = 790
#number_of_maf_parts = 4
number_of_maf_parts = 27