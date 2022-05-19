import pickle

with open("/home/Timo/pickles/merged_dfs", 'rb') as handle:
    merged_dfs = pickle.load(handle)

with open("/home/Timo/pickles/genes_to_cis_reg", 'rb') as handle:
    genes_to_cis_reg = pickle.load(handle)

with open("/home/Timo/pickles/result", 'rb') as handle:
    result = pickle.load(handle)
