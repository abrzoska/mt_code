import pickle
import pandas as pd
from file_parameters import *
from list_util import *


def get_motifs_regs(cis_reg_name, cis_reg_start, cis_reg_end, motif_df):
    associated_cis_regs = motif_df.loc[
        (motif_df['start'] > cis_reg_start) & (motif_df['stop'] < cis_reg_end), ['motif_alt_id']]
    print(associated_cis_regs.to_numpy())
    return "\t".join(cis_reg_name + associated_cis_regs.to_numpy())


chosen_gene_id = "ENSMUSG00000031583"
coordinate_correction = 33385526

with open(gene_to_cis_reg_dic, "rb") as f:
    gene_to_cis_reg_dic = pickle.load(f)

cis_regs = gene_to_cis_reg_dic[chosen_gene_id]

cis_regs_df = pd.read_csv(cis_reg_file, delimiter='\t', names=["chr", "start", "end", "element_id", "rgb"])

selected_cis_regs_df = cis_regs_df[cis_regs_df["element_id"].isin(flatten_list(cis_regs))]

motif_df = pd.read_csv(motif_file, delimiter='\t')
motif_df['start'] += coordinate_correction
motif_df['stop'] += coordinate_correction

motifs = [get_motifs_regs(row[0], row[1], row[2], motif_df) for row in
            zip(selected_cis_regs_df['element_id'], selected_cis_regs_df['start'], selected_cis_regs_df['end'])]

with open(wrn_cis_reg_motifs, "w+") as f:
    for line in motifs:
        f.write(f"{motif_df}\n")
