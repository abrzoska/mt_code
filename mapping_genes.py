import pandas as pd
import pickle

from variables import *


def get_promoter(gene, strand, start, end, chr):
    if strand == 1:
        return [gene, strand, start - 50000, start, chr]
    return [gene, strand, end, end + 50000, chr]


def get_cis_regs(gene, prom_start, prom_end, cis_regs):
    associated_cis_regs = cis_regs.loc[
        (cis_regs['start'] > prom_start) & (cis_regs['end'] < prom_end), ['element_id']]
    return gene, associated_cis_regs.to_numpy()


genes = pd.read_csv(gene_file, delimiter='\t')
result = [get_promoter(row[0], row[1], row[2], row[3], row[4]) for row in
          zip(genes['Gene stable ID'], genes['Strand'], genes['Gene start (bp)'], genes['Gene end (bp)'],
              genes['Chromosome/scaffold name'])]
genes = pd.DataFrame(result, columns=["gene_id", "strand", "prom_start", "prom_end", "chr"])
genes["chr"] = genes["chr"].apply(lambda x: "chr" + str(x))
genes.to_csv("../gene_promoters.tsv", sep='\t')

cis_regs = pd.read_csv(cis_reg_file, delimiter='\t', names=["chr", "start", "end", "element_id", "rgb"])
df = cis_regs.groupby(by="chr")

genes = genes.groupby(by="chr")
gene_to_cis_reg = {}
for k, gb in genes:
    groups = df.groups
    if k in groups:
        cis_regs_group = df.get_group(k)
        cis_regs = [get_cis_regs(row[0], row[2], row[3], cis_regs_group) for row in
                                    zip(gb['gene_id'], gb['strand'], gb['prom_start'], gb['prom_end'], gb['chr'])]
        for gene_id, cis_associated in cis_regs:
            gene_to_cis_reg[gene_id] = cis_associated

with open('gene_to_cis_reg.pkl', 'wb+') as f:
    pickle.dump(gene_to_cis_reg, f)