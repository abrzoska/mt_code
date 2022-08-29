import sys
import pandas as pd
import subprocess
import pickle
import list_util
from file_parameters import *
import numpy as np

maf_to_bed_path = "/home/abrzoska/abrzoska_backup/bin/x86_64/mafsInRegion"

gene = sys.argv[1]
cut_maf = sys.argv[2]

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

#print(gene_to_cis_reg)
ensembl = gene_to_cis_reg[gene]
#print(ensembl)
first = True
for i in range(number_of_maf_parts):
    cis_reg_part = pd.read_csv(cis_reg_indel_file.replace("NUMBER", str(i)))
    #print(cis_reg_indel_file.replace("NUMBER", str(i)))
    cis_reg_part_filtered = cis_reg_part[cis_reg_part["cis_reg_id"].isin(list_util.flatten_list(ensembl.tolist()))]
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

np.savetxt(f'{gene}_indels.bed', indels_all, delimiter='\t', fmt='%s')

if cut_maf != "no_maf":
    bed_file = f"{gene}.bed"
    with open(bed_file, "w") as gene_bed:
        gene_bed.write(f"{gene_prom_chr}\t{gene_prom_start}\t{gene_prom_end}")
    subprocess.run([maf_to_bed_path, bed_file, f"{gene}_promoter_region.maf", input_maf])
    subprocess.run(["rm", bed_file])

# find promoter region
# find indels
# make list with indels and cisregs
# cut maf?