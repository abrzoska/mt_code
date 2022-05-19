import pandas as pd
import pickle
import os
import list_util
from variables import *

adapted = ["micOch1","criGri1","mesAur1","HLcasCan1","hetGla2","HLfukDam1","HLmarMar1","triMan1","chrAsi1","ornAna2","vicPac2","HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosMut1","panHod1","HLenhLut1","odoRosDiv1","lepWed1","neoSch1","conCri1"]
adapted_label = "hx"
adapted2 = ["hg38","pantro5","gorgor5","panpan2","ponabe2","hlcascan1","hetgla2","hlfukdam1","loxafr3","hlturtru3","hlorcorc1","hldelleu1","lipvex1","phycat1","balacu1","hlbalmys1","hlrhisin1","hlhiparm1","eptfus1","myodav1","myobra1","myoluc2","hlminnat1","hldesrot1"]
adapted_label2 = "ll"
adapted_dict = {
  adapted_label: adapted,
  adapted_label2: adapted2,
}
adapted_labels = list(adapted_dict.keys())

upstream_region = 50000
#Todo Adjust these variables
query_species = "nanGal1"
# useful cut command cut -d'    ' -f-5,8
# filter only for target genes
indel_groups =['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']


#analysis_header = ",".join(["cis_reg_id", "query_species_only_count_del", "query_species_only_count_in", "indel_non_adapted_included_count_del", "indel_non_adapted_included_count_in"] + [item for sublist in [[f"{i}_del", f"{i}_in"] for i in adapted_labels] for item in sublist])
analysis_header = ["cis_reg_id", "query_species_only_count_del", "query_species_only_count_in", "indel_non_adapted_included_count_del", "indel_non_adapted_included_count_in", "hx", "ll"]

cis_regs = pd.read_csv(cis_reg_file, delimiter='\t', names=["chr", "start", "end", "element_id", "rgb"])
df = cis_regs.groupby(by="chr")

def determine_cis_regs_genes(gene_file):
    with open('gene_to_cis_reg.pkl', 'rb') as f:
        gene_to_cis_reg = pickle.load(f)
    df = pd.read_excel(gene_file, index_col=None)
    target_genes = df["GeneID_Mouse"]
    return dict((k, gene_to_cis_reg[k]) for k in target_genes.tolist() if k in gene_to_cis_reg)


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
    with open("/home/Timo/pickles/merged_dfs", 'wb') as handle:
        pickle.dump(merged_dfs, handle)
    with open("/home/Timo/pickles/genes_to_cis_reg", 'wb') as handle:
        pickle.dump(genes_to_cis_reg, handle)
    merged_dfs["cis_reg_id"] = merged_dfs["cis_reg_id"].str.strip("'")
    for gene in genes_to_cis_reg.keys():
        sum_lines = merged_dfs[merged_dfs["cis_reg_id"].isin(list_util.flatten_list(genes_to_cis_reg[gene].tolist()))].sum()
        result.append([gene] + sum_lines.to_list()[1:])
    with open("/home/Timo/pickles/result", 'wb') as handle:
        pickle.dump(result, handle)
    return pd.DataFrame(result, columns=analysis_header)


def get_count_of_cis_reg_per_gene(folder, gene_to_cis_reg):
    result = []
    for indel_group in indel_groups:
        indel_analysis = pd.read_csv(f"{folder}/{indel_group}_analysis.csv")
        for gene in gene_to_cis_reg.keys:
            indels = gene_to_cis_reg[gene]
            selected_rows = indel_analysis.loc[indel_analysis["cis_reg_id"] in indels]
            row = [gene]
            for key in ["query_species_only_count_del", "query_species_only_count_in", "indel_non_adapted_included_count"] + adapted_labels:
                row.append(selected_rows[key].sum())
            result.append(row)
    return result


def get_and_analyze_indels(cis_reg_start, cis_reg_end, cis_reg_id, indel_group):
    adapted_species_dicts = {}
    associated_indels = indel_group.loc[(indel_group['Start'] > cis_reg_start) & (indel_group['End'] < cis_reg_end)]
    associated_indels_group = associated_indels.groupby(by="Batch")
    query_species_only_count_del = 0
    query_species_only_count_in = 0
    indel_non_adapted_included_count_del = 0
    indel_non_adapted_included_count_in = 0
    for adapted_label in adapted_labels:
            adapted_species_dicts[f"{adapted_label}_del"] = 0
            adapted_species_dicts[f"{adapted_label}_in"] = 0
    for key, group in associated_indels_group:
        lines_with_this_end_number = group
        ## 1) InDel ist NUR in Spalax, sonst keiner Spezies
        if len(lines_with_this_end_number) == 1:
            only_query_species_with_this_end_number = lines_with_this_end_number.iloc[0]
            if query_species in only_query_species_with_this_end_number.Name:
                type = only_query_species_with_this_end_number.Name.split(".")[-2]
                if type == "del":
                    query_species_only_count_del += 1
                else:
                    query_species_only_count_in += 1
        ##InDel ist in Spalax und ...
        #TODO: maybe do this a bn
        else:
            names_list = lines_with_this_end_number['Name'].tolist()
            if any(query_species in s for s in names_list):
                ## 2) InDel ist in Spalax UND einer Spezies aus der Outgroup (egal ob long lived/adaptiert)
                type = names_list[0].split(".")[-2]
                if type == "del":
                    indel_non_adapted_included_count_del += 1
                else:
                    indel_non_adapted_included_count_in += 1
                ## 3) InDel ist in Spalax UND einer Hypoxie adaptierten Spezies, egal ob Ingroup oder Outgroup
                ## 4) InDel ist in Spalax UND einer langlebigen Spezies , egal ob Ingroup oder   Outgroup
                for adapted_label in adapted_labels:
                    if adapted_label in adapted_species_dicts.keys():
                        if type == "del":
                            adapted_species_dicts[f"{adapted_label}_del"] += 1
                        else:
                            adapted_species_dicts[f"{adapted_label}_in"] += 1
    return ",".join(list(map(lambda x: str(x) for x in [cis_reg_id, query_species_only_count_del, query_species_only_count_in,
                     indel_non_adapted_included_count_del, indel_non_adapted_included_count_in] + list(adapted_species_dicts.values())))) + "\n"


def run_analysis(folder, cis_regs_groups):
    columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand,", "ThickStart", "ThickEnd", "ItemRGB"]
    for key, cis_reg_group in cis_regs_groups:
        if key in indel_groups:
            indel_file_name = f"{folder}/{key}_indels.bed"
            indels = pd.read_csv(indel_file_name, sep="\t", skiprows=0, header=None, names=columns)
            indels = indels.drop_duplicates(subset=['Name'])
            indels["Batch"] = indels["Name"].str.split(".")
            indels["Batch"] = indels["Batch"].apply(lambda x: x[-1] if len(x) > 0 else 0)
            stats = [get_and_analyze_indels(cis_reg[0], cis_reg[1], cis_reg[2], indels) for cis_reg in zip(cis_reg_group["start"],
                                                                    cis_reg_group["end"], cis_reg_group["element_id"])]
            analysis_file = open(f"{folder}/{key}_analysis.csv", "w")
            analysis_file.write(analysis_header + "\n")
            analysis_file.writelines(stats)
            analysis_file.close()
        else:
            print("Key not found: " + key)


def main():
    #run_analysis("/moreAddSpace/tlin", df)
    genes_to_cis_reg = determine_cis_regs_genes(gene_input_file)
    summarize_cis_regs(genes_to_cis_reg, cis_reg_indel_folder).to_csv("/moreAddSpace/tlin/results/genes_csv")

main()
