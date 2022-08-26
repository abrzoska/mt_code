import pandas as pd
from file_parameters import *

def process(df):
    for cis_reg in cis_regs
    associated_cis_regs = cis_regs.loc[(cis_regs['start'] > prom_start) & (cis_regs['end'] < prom_end), ['element_id']]]

motif_file = "/home/tlin/mt_code_data/data/familial_binding_sites_19102021.bed"


motifs = open(motif_file)

chunksize = 10 ** 6
with pd.read_csv(motif_file, chunksize=chunksize) as reader:
    for chunk in reader:
        process(chunk)


