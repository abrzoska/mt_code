import json
from Bio import Entrez
import sys
import os
import hashlib
import pickle
#TODO email übergeben
Entrez.email = "abrzoska@students.uni-mainz.de"

def get_list_of_aliases(gene_name):
    if(len(gene_name)) > 0:
        preferred_entry = Entrez.esearch(db="gene",term="{gene_name} [Preferred Symbol] AND 9606 [Taxonomy ID]".format(gene_name=gene_name),retmode="json")
        preferred_entry_json = json.loads(preferred_entry.read())
        ids = preferred_entry_json['esearchresult']['idlist']
        #print(len(ids))
        if ids == []:
            print("Not Found with Preferred Symbol.")
            all_ids_entries =  Entrez.esearch(db="gene",term="{gene_name} and 9606 [taxonomy id]".format(gene_name=gene_name),retmode="json")
            all_ids_json = json.loads(all_ids_entries.read())
            all_ids = all_ids_json['esearchresult']['idlist']
            print("\tOther candidates:")
            print(all_ids)
            # print all_ids
            for id in all_ids:
                record_with_aliases = Entrez.efetch(db="gene",id=id,retmode="json")
                aliases = []
                for line in record_with_aliases:
                    if line.startswith("official symbol:"):
                        aliases.append(line.split("and name")[0].split(":")[1])
                    if line.startswith("other aliases:"):
                        for x in [y.strip() for y in [x.strip() for x in line.split(":")[1:]][0].split(",")]:
                            aliases.append(x)
                        # print aliases
                        if gene_name in aliases:
                            return aliases[0]
                        elif gene_name.replace(" ","") in aliases:
                            return aliases[0]
                        elif gene_name.replace("-","") in aliases:
                            return aliases[0]
                        else:
                            continue
            return "not found"
        else:
            # 
            #print("________________")
            #print("ORIGINAL Gene ID:"+ str(ids))
            #print("Gene ID:"+ str(ids))
            #print("Gene Name:"+ gene_name)
            for id in ids:
                record_with_aliases = Entrez.efetch(db="gene",id=id,retmode="json")
                #print(record_with_aliases.read())
                for line in record_with_aliases:
                    if line.startswith("Other Aliases:"):
                        unclean_aliases = line.split(":")[1:]
                        
                        aliases = unclean_aliases[0].strip().split(",")
                        #print(aliases)
                        return(aliases)
            
def get_dict_from_gene_list(file, dict):
    #
    dict_path = os.getcwd()+"/"+"{0}.pkl".format(dict)
    print(dict_path)
    if os.path.exists(dict_path) and os.path.getsize(dict_path) != 0:
        print("Using found dictionary {0}".format(dict))
        dict_file = open(dict_path, "rb")
        dict_of_tfs = pickle.load(dict_file)
        dict_file.close()      
    else:
        print("No dictionary with name {0} found. Creating new dictionary.".format(dict))
        dict_of_tfs = {}
    #      
    with open(file, "r")as lf:
        for gene_name_unclean in lf:
            if gene_name_unclean.isspace():
                continue
            gn = gene_name_unclean.split()
            gene_name = gn[0].strip()
            print("Checking {0}".format(gene_name))
            if gene_name in dict_of_tfs:
                print("Gene entry already exists: {0}.".format(gene_name))
            elif not len(gene_name) > 0:
                print("Empty entry: {0}.".format(gene_name))
            else:                
                #print("This should not be none: {0}", test)
                res_unclean = get_list_of_aliases(gene_name)
                if(len(res_unclean)>0):
                    res = list(map(str.strip,res_unclean))
                    dict_of_tfs["_{0}".format(gene_name)] = gene_name
                    res = res  + [gene_name]
                    for i in range(len(res)):                    
                        split_alias_list = res[:i] + res[i+1:]
                        #print('Main: {0} Rest ->   {1}'.format(res[i], split_alias_list))
                        dict_of_tfs[res[i]] = split_alias_list
                        #dict_of_tfs[res[i]] = res
                        #print('Main: {0} Rest ->   {1}'.format(res[i], res))
    #print(dict_of_tfs.keys())
    dict_file = open(dict_path, "wb")
    pickle.dump(dict_of_tfs,dict_file)
    lf.close()
    
def remove_duplicates(in_path, out_path):
    exists_hash = set()
    outfile = open(out_path, "w")
    for line in open(in_path, "r"):
        hash = hashlib.md5(line.strip().encode('utf-8')).hexdigest()
        if hash not in exists_hash:
            outfile.write(line)
            exists_hash.add(hash)
            
def get_tfmotifview_tfs_from_bed(tfmv_in, tfmv_out):
    outfile = open(tfmv_out, "w")
    for line in open(tfmv_in, "r"):        
        line = (line.split()[3]).split(".")[0]
        outfile.write(line+"\n")

def check_which_tfs_are_in_common(dict_name,tfmv_out, cmp_result):
    dict_path = os.getcwd()+"/"+"{0}.pkl".format(dict_name)
    dict_file = open(dict_path, "rb")
    dict = pickle.load(dict_file)
    only_in_dict = []
    only_in_predictions = []
    present_in_both = []
    list_of_predicted_tfs = []
    for gene_name in open(tfmv_out, "r"):        
        gn = gene_name.strip()
        list_of_predicted_tfs.append(gn)
        if gn in dict:            
            present_in_both.append(gn)
        else:
            only_in_predictions.append(gn)
    for key in dict.keys():
        if key.startswith("_"):            
            val = dict[key] 
            if val not in list_of_predicted_tfs:
                only_in_dict.append(val)
    print("Writing output TP/FP/FN.")
    outfile = open(cmp_result, "w")

    outfile.write("\nPresent in both sets (True positive)\n\n")    
    for i in present_in_both:
        outfile.write("{0} ({1})\n".format(i, ",".join(dict.get("{0}".format(i))))) 
        
    outfile.write("\nOnly in dictionary, based on experiments (False Negative)\n\n")
    for i in only_in_dict:
        outfile.write("{0} ({1})\n".format(i,",".join(dict.get("{0}".format(i)))))  
        
    outfile.write("\nOnly in predictions, based on algorithms (possible False Positive)\n\n")    
    for i in only_in_predictions:
        outfile.write("{0}\n".format(i))      

        

    outfile.close()
def main(infile, outfile, dict, tfmv_in, tfmv_tf, tfmv_out, comparison_result):  
    remove_duplicates(infile, outfile)
    get_dict_from_gene_list(outfile, dict)
    get_tfmotifview_tfs_from_bed(tfmv_in, tfmv_tf)
    remove_duplicates(tfmv_tf, tfmv_out)
    check_which_tfs_are_in_common(dict, tfmv_out, comparison_result) 
        
if __name__ == "__main__":  
    main("wrn_encode.txt", "wrn_encode_clean.txt", "arc", "wrn_TF.txt", "wrn_TF_.txt","wrn_TF_clean.txt", "results.txt")    