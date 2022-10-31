import os

run_name = "indels2"

#TODO: make run directories

def create_dic_if_non_existent(path):
    exists = os.path.exists(path)
    if not exists:
        os.makedirs(path)

create_dic_if_non_existent(f"/mnt/moreAddspace/mt_code_data/intermediary_data/{run_name}")
create_dic_if_non_existent(f"/mnt/moreAddspace/mt_code_data/debug/{run_name}")
create_dic_if_non_existent(f"/home/tlin/mt_code_data/intermediary_data/{run_name}/indels/")
create_dic_if_non_existent(f"/home/tlin/mt_code_data/intermediary_data/{run_name}/cis_reg_indels/")
create_dic_if_non_existent(f"/home/tlin/mt_code_data/debug/{run_name}/")

#inputs
gene_input_file = "/mnt/moreAddspace/mt_code_data/data/Cross_all_final_200421.xlsx"

#data
gene_promoter_file = "/mnt/moreAddspace/mt_code_data/data/gene_promoters.tsv"
cis_reg_file = "/mnt/moreAddspace/mt_code_data/data/encode9_trimmed.bed"
gene_file = "/mnt/moreAddspace/mt_code_data/data/mart_export.txt"
input_maf = "/mnt/moreAddspace/abrzoska/mm10/multi.anno.maf"
input_maf_part = "/home/tlin/maf_split/maf_part_NUMBER.maf"
maf_split = "/home/tlin/maf_split/"
#cis_reg_pickle = f"/home/tlin/mt_code_data/debug/genes_to_cis_reg.pkl"
gene_to_cis_reg_dic = f"/home/tlin/mt_code_data/data/gene_to_cis_reg.pkl"

#intermdiary
indel_file = f"/home/tlin/mt_code_data/intermediary_data/maf_test/indels/maf_part_NUMBER.bed"
indel_file_name = f"/home/tlin/mt_code_data/intermediary_data/maf_test/all_chromosomes_0.bed"
#bed_out = f"/mnt/moreAddspace/mt_code_data/intermediary_data/{run_name}/indels/all_chromosomes"
cis_reg_indel_folder = f"/home/tlin/mt_code_data/intermediary_data/{run_name}/cis_reg_indels/"
cis_reg_indel_file = f"/home/tlin/mt_code_data/intermediary_data/{run_name}/cis_reg_indels/maf_part_NUMBER_analysis.csv"
intermediary_folder = f"/home/tlin/mt_code_data/intermediary_data/{run_name}/"

#result
gene_result_csv = f"/home/tlin/mt_code_data/result/{run_name}_genes.csv"
gene_result_summary_txt = f"/home/tlin/mt_code_data/result/{run_name}_summary.txt"

#pickles
merged_df_pickle = f"/home/tlin/mt_code_data/debug/{run_name}/merged_dfs.pkl"
result_pickle = f"/home/tlin/mt_code_data/debug/{run_name}/result.pkl"
cis_reg_pickle = f"/home/tlin/mt_code_data/debug/{run_name}/genes_to_cis_reg.pkl"

#test
motif_file = "/mnt/moreAddspace/mt_code_data/test_data/wrn_fimo.tsv"
wrn_cis_reg_motifs = "/mnt/moreAddspace/mt_code_data/test_data/wrn_result.tsv"

#input_maf = "/mnt/moreAddspace/mt_code_data/test_data/wrn_promoter.maf"
#maf_split = "/mnt/moreAddspace/mt_code_data/test_data/wrn_maf_split/"
#input_maf_part = "/mnt/moreAddspace/mt_code_data/test_data/wrn_maf_split/maf_part_NUMBER.maf"


#meta #TODO: computationally determine this and include the cutter script in the pipeline
number_of_maf_parts = 790
#number_of_maf_parts = 4
#number_of_maf_parts = 27
