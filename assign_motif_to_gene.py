import pandas as pd
from file_parameters import *
import pickle

#https://pythonspeed.com/articles/chunking-pandas/
motif_file = "/home/tlin/mt_code_data/data/familial_binding_sites_19102021.bed"
gene = "" #Todo add gene id

with open(gene_to_cis_reg_dic, 'rb') as f:
    gene_to_cis_reg = pickle.load(f)



def process(df):
    for cis_reg in cis_regs:
        associated_cis_regs = df.loc[(df['start'] > prom_start) & (df['end'] < prom_end), ['element_id']]


chunksize = 10 ** 6
ensembl = gene_to_cis_reg[gene]
for chunk in pd.read_csv(motif_file, chunksize=chunksize): #Todo: add more read parameters
    process(chunk)


