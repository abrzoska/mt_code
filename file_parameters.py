import os

run_name = "indels_3"

#TODO: make run directories

def create_dic_if_non_existent(path):
    exists = os.path.exists(path)
    if not exists:
        os.makedirs(path)

create_dic_if_non_existent(f"/lustre/project/m2_jgu-nannospalax/data/mt_code_intermediary_data/{run_name}")
create_dic_if_non_existent(f"/lustre/project/m2_jgu-nannospalax/data/mt_code_debug/{run_name}")

#inputs
gene_input_file = "/lustre/project/m2_jgu-nannospalax/data/mt_code_data/Cross_all_final_200421.xlsx"

#data
gene_promoter_file = "/lustre/project/m2_jgu-nannospalax/data/mt_code_data/gene_promoters.tsv"
cis_reg_file = "/lustre/project/m2_jgu-nannospalax/data/mt_code_data/encode9_trimmed.bed"
gene_file = "/lustre/project/m2_jgu-nannospalax/data/mt_code_data/mart_export2.txt"
input_maf = "/lustre/project/m2_jgu-nannospalax/data/mt_code_data/multi.anno.maf"
input_maf_part = "/lustre/project/m2_jgu-nannospalax/data/mt_code_data/maf_split/maf_part_NUMBER.maf"
maf_split = "/lustre/project/m2_jgu-nannospalax/data/mt_code_data/maf_split/"

#intermdiary
indel_file = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_intermediary_data/{run_name}/indels/maf_part_NUMBER.bed"
indel_file_name = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_intermediary_data/{run_name}/all_chromosomes_0.bed"
#bed_out = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_intermediary_data/{run_name}/indels/all_chromosomes"
cis_reg_indel_folder = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_intermediary_data/{run_name}/cis_reg_indels/"
intermediary_folder = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_intermediary_data/{run_name}/"

#result
gene_result_csv = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_results/{run_name}_genes_csv"

#pickles
merged_df_pickle = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_debug/{run_name}/merged_dfs.pkl"
cis_reg_pickle = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_debug/{run_name}/genes_to_cis_reg.pkl"
result_pickle = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_debug/{run_name}/result.pkl"
gene_to_cis_reg_dic = f"/lustre/project/m2_jgu-nannospalax/data/mt_code_debug/{run_name}/gene_to_cis_reg.pkl"

#scripts
splitter_script = "/home/timlin/projects/mt_code/chromosome_splitter.sh"

#test
motif_file = "/lustre/project/m2_jgu-nannospalax/data/mt_code_test_data/wrn_fimo.tsv"
wrn_cis_reg_motifs = "/lustre/project/m2_jgu-nannospalax/data/mt_code_test_data/wrn_result.tsv"
#input_maf = "/lustre/project/m2_jgu-nannospalax/data/mt_code_test_data/ENSMUSG00000025475-indels.maf"

#meta #TODO: computationally determine this and include the cutter script in the pipeline
number_of_maf_parts = 790
#number_of_maf_parts = 4