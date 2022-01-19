import csv
from pathlib import Path
import sys
import pandas as pd
import os
import pickle
import pandas as pd


mouse_dict_path = os.getcwd()+"/"+"mouse_genes_lt_0.5.pkl"
dict_file = open(mouse_dict_path, "rb")
mouse_dict_of_genes = pickle.load(dict_file)
b_out = open("mouse_genes_list.txt", "w")
for key, value in mouse_dict_of_genes.items():
    if len(value) > 3: 
        print(f'{key}: {value}')    
        b_out.write("{0}\n".format(key))
b_out.close()

def create_dict_for_mouse():
    translation_of_genes = {}
    df = pd.read_excel (r'../listen/Cross_all_final_200421.xlsx')
    #iterate over rows with iterrows()
    for index, row in df.iterrows():
         # access data using column names
        #print(index, row['GeneID_Rat'], row['GeneID_Mouse'])
        translation_of_genes[row['GeneID_Rat']] = row['GeneID_Mouse']
    dict_path = os.getcwd()+"/"+"genes_lt_0.5.pkl"
    dict_file = open(dict_path, "rb")
    dict_of_genes = pickle.load(dict_file)

    mouse_dict_path = os.getcwd()+"/"+"mouse_genes_lt_0.5.pkl"
    mouse_dict_of_genes = {}
    try:
        for key in dict_of_genes.keys():
            mouse_name = translation_of_genes[key]
            mouse_dict_of_genes[mouse_name]=dict_of_genes[key]
    except KeyError:
        print("Key {0} was not found.".format(key))
    mouse_dict_file = open(mouse_dict_path, "wb")
    pickle.dump(mouse_dict_of_genes,mouse_dict_file)
def is_float(value):
  try:
    float(value)
    return True
  except:
    return False
def create_dict_from_csv():
    path = sys.argv[1]
    organ = sys.argv[2]
    if Path(path).is_file():
        n = "".join(path.split(".")[0:-1])
        c = n.split("/")[-1]    
        output_file = "{0}.txt".format(c)
        dict_path = os.getcwd()+"/"+"genes_lt_0.5.pkl"
        if os.path.exists(dict_path) and os.path.getsize(dict_path) != 0:
            print("Using found dictionary {0}".format(dict))
            dict_file = open(dict_path, "rb")
            dict_of_genes = pickle.load(dict_file)
            dict_file.close()      
        else:
            print("No dictionary with name {0} found. Creating new dictionary.".format(dict))
            dict_of_genes = {}
                
        with open(path, newline='') as csvfile:
            r = csv.reader(csvfile, delimiter=',')        
            for row in r:
                #print(row[-1])
                if is_float(row[-1]) and float(row[-1]) < 0.05:
                    gene_name = row[0]
                    if gene_name in dict_of_genes:
                        if organ not in dict_of_genes[gene_name]:
                            dict_of_genes[gene_name].append(organ)
                    else:   
                        dict_of_genes[gene_name] = [organ]
                    print(row[0])
            print(dict_of_genes)
            dict_file = open(dict_path, "wb")
            pickle.dump(dict_of_genes,dict_file)