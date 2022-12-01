import sys
import pandas as pd
import subprocess
import pickle

import file_parameters

import list_util
from file_parameters import *
import numpy as np
import os
from variables import *

def get_relevant_indels_for_cis_reg(cis_reg_start, cis_reg_end, cis_reg_id, indel_file, min_indel_size, gene_id):
    adapted_species_dicts = {}
    columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand,", "ThickStart", "ThickEnd", "ItemRGB", "blockCount"]
    indels = pd.read_csv(indel_file, sep="\t", skiprows=0, header=None, names=columns)
    indels = indels.drop_duplicates(subset=['Name'])
    indels["Batch"] = indels["Name"].str.split(".")
    indels["Batch"] = indels["Batch"].apply(lambda x: x[-1] if len(x) > 0 else 0)
    indel_group = indels
    associated_indels = indel_group.loc[(indel_group['Start'] > cis_reg_start) & (indel_group['End'] < cis_reg_end) & (indel_group['blockCount'] >= min_indel_size)]
    associated_indels_group = associated_indels.groupby(by="Batch")
    has_indels = False
    indel_lines = []
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
                indel_lines.append(list(only_query_species_with_this_end_number) + [cis_reg_id, gene_id, "Spalax Specific"])
        ##InDel ist in Spalax und ...
        #TODO: maybe do this a bn
        else:
            spalax_line = lines_with_this_end_number.iloc[0]
            #Note: this means that if one species is in another group this is reported
            names_list = lines_with_this_end_number['Name'].tolist()
            ## 3) InDel ist in Spalax UND einer Hypoxie adaptierten Spezies, egal ob Ingroup oder Outgroup
            ## 4) InDel ist in Spalax UND einer langlebigen Spezies , egal ob Ingroup oder   Outgroup
            for adapted_label in adapted_labels:
                if any(adapted_label in s for s in names_list):
                    has_indels = True
                    indel_lines.append(list(spalax_line) + [cis_reg_id, gene_id, f"Adapted group: {adapted_label}"])
    if not has_indels:
        return []
    with open(f"{cis_reg_id}_called_indels_only.bed", "w") as f:
        for indel_line in indel_lines:
            f.write("\t".join(indel_line) + "\n")
    return indel_lines

script_path_maf_in_region = "/home/abrzoska/abrzoska_backup/bin/x86_64/mafsInRegion"
script_path_indels_for_cis_reg = "/mnt/spalax2/get_indels.sh"

gene = sys.argv[1]
cut_maf = sys.argv[2]
subprocess.run(["mkdir", gene])
os.chdir(gene)

try:
    min_indel_length = sys.argv[3]
except:
    print("Change in script usage: please specify min indel length as last parameter (e.g. 2)")
    exit(1)

indel_for_gene_file = f"{gene}_indels.bed"

run_name = f"indels{min_indel_length}"
run_variables = RunVariables(run_name)

gene_promoters = pd.read_csv(gene_promoter_file, sep="\t")
gene_promoters = gene_promoters.drop_duplicates(subset=["gene_id"])
gene_prom = gene_promoters.loc[gene_promoters["gene_id"] == gene]

try:
    gene_prom_chr = gene_prom[["chr"]].iloc[0]["chr"]
except:
    print(f"gene not found: {gene}")

gene_prom_start = gene_prom[["prom_start"]].iloc[0]["prom_start"]
gene_prom_end = gene_prom[["prom_end"]].iloc[0]["prom_end"]

with open(gene_to_cis_reg_dic, 'rb') as f:
    gene_to_cis_reg = pickle.load(f)

ensembl = list_util.flatten_list(gene_to_cis_reg[gene].tolist())

first = True
for i in range(number_of_maf_parts):
    cis_reg_part = pd.read_csv(run_variables.cis_reg_indel_file.replace("NUMBER", str(i)))
    cis_reg_part_filtered = cis_reg_part[cis_reg_part["cis_reg_id"].isin(ensembl)]
    if first:
        first = False
        cis_regs_all = cis_reg_part_filtered
    else:
        cis_regs_all = pd.concat([cis_regs_all, cis_reg_part_filtered])

cis_regs_all.to_csv(f"{gene}_cis_regs.csv")

first = True
for i in range(number_of_maf_parts):
    columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand,", "ThickStart", "ThickEnd", "ItemRGB", "blockCount"]
    indels_part_df = pd.read_csv(indel_file.replace("NUMBER", str(i)), sep="\t", skiprows=0, header=None, names=columns)
    indels_part = indels_part_df.loc[(indels_part_df['Start'] > gene_prom_start) & (indels_part_df['End'] < gene_prom_end)
                                     & (indels_part_df["Chromosome"] == gene_prom_chr)]
    if first:
        first = False
        indels_all = indels_part
    else:
        indels_all = pd.concat([indels_all, indels_part])


np.savetxt(indel_for_gene_file, indels_all, delimiter='\t', fmt='%s')

all_cis_regs = pd.read_csv(cis_reg_file, delimiter='\t', names=["chr", "start", "end", "element_id", "rgb"])
all_called_indels_for_gene = []
for cis_reg_id in ensembl:
    subprocess.run([script_path_indels_for_cis_reg, cis_reg_id, indel_for_gene_file])
    cis_reg_line = all_cis_regs.loc[all_cis_regs["element_id"] == cis_reg_id]
    all_called_indels_for_gene += (get_relevant_indels_for_cis_reg(int(cis_reg_line["start"]), int(cis_reg_line["end"]), cis_reg_id, indel_for_gene_file, int(min_indel_length), gene))

with open(f"{gene}_called_indels_only.bed", "w") as f:
    for indel_line in all_called_indels_for_gene:
        f.write("\t".join(indel_line) + "\n")

if cut_maf == "with_maf":
    bed_file = f"{gene}.bed"
    with open(bed_file, "w") as gene_bed:
        gene_bed.write(f"{gene_prom_chr}\t{gene_prom_start}\t{gene_prom_end}")
    subprocess.run([script_path_maf_in_region, bed_file, f"{gene}_promoter_region.maf", input_maf])
    subprocess.run(["rm", bed_file])

# find promoter region
# find indels
# make list with indels and cisregs
# cut maf?
