import pandas as pd
import pickle
from file_parameters import *
from variables import *
from memory_profiler import profile
from joblib import Parallel, delayed
import list_util


def determine_cis_regs_genes(gene_file):
    with open(gene_to_cis_reg_dic, 'rb') as f:
        gene_to_cis_reg = pickle.load(f)
    df = pd.read_excel(gene_file, index_col=None)
    target_genes = df["GeneID_Mouse"]
    return dict((k, gene_to_cis_reg[k]) for k in target_genes.tolist() if k in gene_to_cis_reg)

#Todo parrallel?
@profile
def summarize_cis_regs(genes_to_cis_reg, out_folder):
    first = True
    files = os.listdir(out_folder)
    result = []
    for file in files:
        file_path = os.path.join(out_folder, file)
        df = pd.read_csv(file_path)
        df = df.groupby(['cis_reg_id'], as_index=False).sum()
        if first:
            merged_dfs = df
            first = False
        else:
            merged_dfs = pd.merge(merged_dfs, df, how='outer')
    print(merged_dfs)
    with open(merged_df_pickle, 'wb') as handle:
        pickle.dump(merged_dfs, handle)
    with open(cis_reg_pickle, 'wb') as handle:
        pickle.dump(genes_to_cis_reg, handle)
    merged_dfs["cis_reg_id"] = merged_dfs["cis_reg_id"].str.strip("'")
    for gene in genes_to_cis_reg.keys():
        sum_lines = merged_dfs[merged_dfs["cis_reg_id"].isin(list_util.flatten_list(genes_to_cis_reg[gene].tolist()))].sum()
        result.append([gene] + sum_lines.to_list()[1:])
    with open(result_pickle, 'wb') as handle:
        pickle.dump(result, handle)
    return pd.DataFrame(result, columns=analysis_header)


def summarize_genes(filename, summary_file_name):
    genes_of_interest = pd.read_csv(filename)
    gene_summary_file = open(summary_file_name, "w+")
    for key in ["query_species_only_count_del", "query_species_only_count_in"] + list_util.flatten_list([[f"{adapted_label}_del", f"{adapted_label}_in"] for adapted_label in adapted_labels]):
        sum =  genes_of_interest[genes_of_interest[key]!=0].sum()[key]
        gene_summary_file.write(f"# of genes with {key}: {sum}\n")
    gene_summary_file.close()


def get_and_analyze_indels(cis_reg_start, cis_reg_end, cis_reg_id, indel_group, min_indel_size):
    adapted_species_dicts = {}
    associated_indels = indel_group.loc[(indel_group['Start'] > cis_reg_start) & (indel_group['End'] < cis_reg_end) & indel_group['blockCount'] >= min_indel_size]
    associated_indels_group = associated_indels.groupby(by="Batch")
    query_species_only_count_del = 0
    query_species_only_count_in = 0
    has_indels = False
    for adapted_label in adapted_labels:
            adapted_species_dicts[f"{adapted_label}_del"] = 0
            adapted_species_dicts[f"{adapted_label}_in"] = 0
    for key, group in associated_indels_group:
        lines_with_this_end_number = group
        ## 1) InDel ist NUR in Spalax, sonst keiner Spezies
        if len(lines_with_this_end_number) == 1:
            only_query_species_with_this_end_number = lines_with_this_end_number.iloc[0]
            if query_species in only_query_species_with_this_end_number.Name:
                has_indels = True
                type = only_query_species_with_this_end_number.Name.split(".")[-2]
                if type == "del":
                    query_species_only_count_del += 1
                else:
                    query_species_only_count_in += 1
        ##InDel ist in Spalax und ...
        #TODO: maybe do this a bn
        else:
            #Note: this means that if one species is in another group this is reported
            names_list = lines_with_this_end_number['Name'].tolist()
            type = names_list[0].split(".")[-2]
            ## 3) InDel ist in Spalax UND einer Hypoxie adaptierten Spezies, egal ob Ingroup oder Outgroup
            ## 4) InDel ist in Spalax UND einer langlebigen Spezies , egal ob Ingroup oder   Outgroup
            for adapted_label in adapted_labels:
                if any(adapted_label in s for s in names_list):
                    has_indels = True
                    if type == "del":
                        adapted_species_dicts[f"{adapted_label}_del"] += 1
                        break
                    else:
                        adapted_species_dicts[f"{adapted_label}_in"] += 1
                        break
    if not has_indels:
        return ""
    return ",".join(list(map(lambda x: str(x), [cis_reg_id, query_species_only_count_del, query_species_only_count_in] + list(adapted_species_dicts.values())))) + "\n"

@profile
def run_analysis(result_folder, indel_file, min_indel_size):
    columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand,", "ThickStart", "ThickEnd", "ItemRGB", "blockCount"]
    cis_regs_groups = pd.read_csv(cis_reg_file, delimiter='\t', names=["chr", "start", "end", "element_id", "rgb"])
    cis_regs_groups = cis_regs_groups.groupby("chr")
    indels = pd.read_csv(indel_file, sep="\t", skiprows=0, header=None, names=columns)
    indels = indels.drop_duplicates(subset=['Name'])
    indels["Batch"] = indels["Name"].str.split(".")
    indels["Batch"] = indels["Batch"].apply(lambda x: x[-1] if len(x) > 0 else 0)
    indels = indels.groupby("Chromosome")
    all_stats = []
    for key, indel_group in indels:
        if key in cis_regs_groups.groups.keys():
            cis_reg_group = cis_regs_groups.get_group(key)
            stats = [get_and_analyze_indels(cis_reg[0], cis_reg[1], cis_reg[2], indel_group, min_indel_size) for cis_reg in
                     zip(cis_reg_group["start"], cis_reg_group["end"], cis_reg_group["element_id"])]
            stats = list(filter(lambda x: x != "", stats))
            all_stats += stats
        else:
            print(f"Key not found: {key}")
    result_name = indel_file.split("/")[-1].replace(".bed", "_analysis.csv")
    analysis_file = open(f"{result_folder}{result_name}", "w")
    analysis_file.write(",".join(analysis_header) + "\n")
    analysis_file.writelines(all_stats)
    analysis_file.close()


def run_parrallel_analysis(result_folder, indel_file, min_indel_size):
    Parallel(n_jobs=number_of_cores)(delayed(run_analysis)(result_folder, indel_file.replace("NUMBER", str(i))
                                                           , min_indel_size) for i in range(number_of_maf_parts))


@profile
def main(min_indel_size):
    #run_parrallel_analysis(cis_reg_indel_folder, indel_file, min_indel_size)
    genes_to_cis_reg = determine_cis_regs_genes(gene_input_file)
    summarize_cis_regs(genes_to_cis_reg, cis_reg_indel_folder).to_csv(gene_result_csv)
