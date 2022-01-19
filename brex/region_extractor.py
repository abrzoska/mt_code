import logging
import requests
import os
import sys
import pyranges as pr
from brex import helper as h
from pathlib import Path
import re
import subprocess
from difflib import SequenceMatcher

class RegionExtractor:
    MIN_INDEL_SIZE = 4
    TMP_FOLDER = "/tmp"
    SIMILARITY = 0.8
    #TODO ddefine a overlap to feather below case
    # A-----A
    # -------
    def __init__(self):
        self.helper = h.Helper()
        logging.basicConfig(filename='brex.log', level=logging.DEBUG)
        logging.info("RegionExtractor instantiated. ")
    ##   find candidates within range for bfile 
    def find_bed_candidates(self, bfile, bfile_name, scaffold, start, end, output_dir):
        gr = self.helper.read_bedfile(bfile)
        ##s = self.helper.remove_quotation_marks(scaffold)
        s=scaffold
        print("Scaffold {0}".format(scaffold))
        cand = gr[s,start:end].as_df()
        output = self.get_output_name(bfile_name, scaffold, start, end, output_dir)
        self.write_bed_cand_to_output(cand, output)
        return output        
    ##   write results from dataframe to file    
    def write_bed_cand_to_output(self, cand, path):
        with open(path, 'w') as f:            
            candAsString = cand.to_string(header = True, index = False)
            f.write(candAsString)  
    ##   puts out a bed and a maf (query and target species only) that only  
    ##   include the indels contained in the maf
    ##   TODO indicate in-group
    def find_indels_for_query_from_maf(self, maf, target, query, in_group, adapted_species, maf_out, bed_out):
        m_out = open(maf_out, "w")
        b_out = open(bed_out, "w")
        del_number = 1
        in_number = 1        
        scaffold = ""    
        is_target_in_in_group = target in in_group
        b_out.write(f'track name="Indels for {query}" itemRgb="On"\n')
        target_line = ""
        query_line = ""
        score_line = ""   
        out_block_lines = []
        in_block_lines = []
        out_block_seqs = []
        in_block_seqs = []        
        adapted_species_lines = []
        species_list = ""
        with open(maf, "r")as mf:
            for line in mf:           
                if line.startswith("##"):
                    m_out.write(line)
                if line.startswith("a"):
                    score_line = line
                if line.startswith("s"):
                    species = self.get_species_from_maf_line(line)
                    if species != -1: 
                        if species == query:
                            query_line = line                          
                        elif species == target:
                            target_line = line                          
                        elif species in in_group:
                            in_block_lines.append(line)
                        else:
                            out_block_lines.append(line) 
                    else:
                        header = header + line             
                ## block is over.          
                if line == "\n":
                    if len(query_line) > 0 and len(target_line) > 0 and self.is_indel(target_line, query_line):   
                    ##If both lines are filled (= the block contains the query species), proceed                 
                    ## if the block contains query species and query species contains indel
                        ##add in target line
                        if is_target_in_in_group:
                            in_block_lines.append(target_line)
                            in_block_seqs.append(self.get_sequence_from_maf_line(target_line))
                        else: 
                            out_block_lines.append(target_line)   
                        ##
                        for line in in_block_lines:
                            in_block_seqs.append(self.get_sequence_from_maf_line(line))
                        for idx, line in enumerate(out_block_lines):
                            out_block_seqs.append(self.get_sequence_from_maf_line(line))
                            species = self.get_species_from_maf_line(line)
                            if species in adapted_species:
                                adapted_species_lines.append(idx)
                        if len(scaffold) == 0:
                            scaffold = self.get_scaffold_from_maf_line(target_line)
                        # if len(adapted_species_lines) > 0:
                            # print("adapted_species_lines")# print(adapted_species_lines) 
                        ##pseudo code
                        ##get gap locations
                        query_seq = self.get_sequence_from_maf_line(query_line)
                        query_gaps = self.get_gap_locations(query_seq)
                        is_indel_present_in_in_group = False
                        for i in query_gaps: 
                            gapsize = abs(i[0]-i[1])
                            if gapsize > self.MIN_INDEL_SIZE:
                            ##check if it present in ingroup
                                for line in in_block_lines:
                                    seq = self.get_sequence_from_maf_line(line)
                                    ##if this gap is shared by the in-group we do not care 
                                    if self.get_similarity(query_seq[i[0]:i[1]], seq[i[0]:i[1]]) > self.SIMILARITY:
                                        is_indel_present_in_in_group = True
                                        break        
                                ##check if present in outgroup#
                                if not is_indel_present_in_in_group:
                                    is_adapted = False
                                    for line in out_block_lines:
                                        seq = self.get_sequence_from_maf_line(line)
                                        if self.get_similarity(query_seq[i[0]:i[1]], seq[i[0]:i[1]]) > self.SIMILARITY:    
                                            ##COOL
                                            ##get bed formation?
                                            ## collect all species !!
                                            current_species = self.get_species_from_maf_line(line)
                                            #species_list = species_list + " ... "+ current_species
                                            #print(line)
                                            ##TODO IS ADAPTED?
                                            if current_species in adapted_species:
                                                is_adapted = True
                                                #print(line)
                                            b_out.write(self.get_bed_line(target_line,line, "del", i[0], i[1], del_number, is_adapted))
                                        is_adapted = False
                                is_indel_present_in_in_group = False
                                #if len(species_list) > 0 :
                                    #print("species_list {0}".format(species_list))
                                #species_list = ""
                                #for line in out_block_lines:
                                    #seq = self.get_sequence_from_maf_line(line)
                        #write to maf file
                            del_number = del_number + 1
                        m_out.write(score_line)
                        m_out.write(target_line)
                        m_out.write(query_line)
                        m_out.write("\n")
                    scaffold = ""
                    target_line = ""
                    query_line = ""
                    out_block_lines = []
                    in_block_lines = []
                    out_block_seqs = []
                    in_block_seqs = []
                    adapted_species_lines = []
                    in_group_species_w_same_indel = []                     
        mf.close()    
        m_out.close()
        b_out.close()   
    def get_bed_line(self, target, line, indel, s, e, number, is_adapted):
        scaffold = self.get_scaffold_from_maf_line(target)
        block_start = self.get_start_from_maf_line(target)
        strand = self.get_strand_from_maf_line(target)
        start = s + int(block_start)
        end = e + int(block_start)       
        name = self.get_species_from_maf_line(line) +"."+ indel + "."+ str(number)
        color = self.get_color(indel, is_adapted)
        return f'{scaffold}\t{start}\t{end}\t{name}\t0\t{strand}\t{start}\t{end}\t{color}\n'
    
    def run(self, start, end, scaffold, maf_in, target_species, query_species, adapted_species, in_group, m_out, b_out):
        ## check all files + log
        complete_path = os.getcwd()+self.TMP_FOLDER 
        print("Location of tmp folder: {0}".format(complete_path))
        if not Path(complete_path).is_dir():
            subprocess.call(['mkdir', complete_path])
            print("Created Folder {0}.".format(self.TMP_FOLDER))
        else:
            print("tmp folder {0} exists.".format(complete_path))            
        ## create bed from start, end
        name = "{0}_{1}_{2}".format(scaffold,start,end)
        bed = complete_path+"/"+name+".bed"
        maf_out = complete_path+"/"+name+".maf"
        if not Path(bed).is_file():   
            print("Writing bed to file .. {0}".format(bed))            
            with open(bed, 'a') as bedfile:
                bedfile.write("{0}\t{1}\t{2}\t{3}".format(scaffold,start,end,name))
            bedfile.close()
        else:
            print("bed file {0} exists.".format(bed))
        if not Path(maf_out).is_file():
            print("Starting mafsInRegion .. {0}".format(maf_out))
            subprocess.call(['mafsInRegion', bed, maf_out, maf_in])
            print("Finished mafsInRegion.")
        else:
            print("maf file {0} exists.".format(maf_out)) 
        ####maf, target, query, in_group, adapted_species, maf_out, bed_out
        
        self.find_indels_for_query_from_maf(maf_out,target_species,query_species, in_group, adapted_species, m_out, b_out)
        print("Maf file containing Indels {0}".format(m_out))
        print("Bed file containing Indels {0}".format(b_out))    
    ##   given a query bed and a list of additions .bed files, get a listing of which elements from the query .bed 
    ##   have overlaps with elements from the target beds      
    ##   track name=HbVar type=bedDetail description="HbVar custom track"
    def annotate_bed_with_overlap(self, query_bed, target_beds):
        target_bed_files = []        
        for file in target_beds:
            tmp_file_name =file+"-tmp"
            self.helper.cut_down_bed_file(file,tmp_file_name , 4)
            target_bed_files.append(self.helper.read_bedfile(tmp_file_name))            
        query_file = open(query_bed, "r")    
        #        
        for line in query_file:
            if(not line.startswith("tr") and not line.startswith("Chr")):   
                scaffold = self.get_scaffold_from_row(line)
                startEnd = self.get_range_from_row(line)
                name = self.get_name_from_row(line)
                frames = []
                if(scaffold != -1 and startEnd != -1 and name != -1):
                    start = startEnd[0]
                    end = startEnd[1]
                    for file in target_bed_files:          
                        cand = file[scaffold,start:end].as_df()
                        #if cand.columns
                        column_name = cand.columns
                        print("sdlsadklsd -{0}-".format(column_name))
                        #cand['Notes'] = cand. + "_" + name
                    ##########
                    ##print to command line
                    ##########
                        # if not cand.empty:
                            # print(str(start) + "-" + str(end)) 
                            # print(cand)
                            # print("#####################")
                            # frames.append(cand)
                    # if len(frames) > 0:
                        # print("#####################")
                        # print(line[0:40])
                        # if len(frames) >= 2:
                            # result= frames[0].append(frames[1])
                        # else:
                            # result=frames[0]
                        # print(result)
                        # print("#####################")
        query_file.close()          
    ## GETTER    
    def get_gap_locations(self, seq):
        gap_locations = []
        matches = list(re.finditer('-+', seq))
        for i, gap in enumerate(matches, 1):
            #gap_locations.append([gap.start(),gap.end()-1])
            gap_locations.append([gap.start(),gap.end()])
        return gap_locations                     
    def get_size(self, maf_line):
        pattern = re.compile("[0-9]")
        l = maf_line.split()
        if(len(l) > 2 and pattern.match(l[3])):        
            return int(l[3])            
    def get_color(self, indel, is_adapted):
        if indel == "del":
            if is_adapted:
                return "128,0,128"
            else:
                return "235,64,52"
        elif indel == "in":
            return "60,222,49"
        else:
            return "133,133,133"            
    def get_range_from_row(self,row):
        ar = row.split()
        if(len(ar)>2 and not ar[0].startswith("tr") and not ar[0].startswith("Chr")):
            return [int(ar[1]),int(ar[2])]
        return -1
    def get_scaffold_from_row(self, row):
        ar = row.split()
        if(len(ar)>0 and not ar[0].startswith("tr") and not ar[0].startswith("Chr")):
            return ar[0]
        return -1    
    def get_name_from_row(self, row):
        ar = row.split()
        if(len(ar)>2 and not ar[0].startswith("tr") and not ar[0].startswith("Chr")):
            return ar[3]
        return -1                          
    ##   returns a list of all species included in the alignment inside the maf 
    ##   and writes them to the outfile 
    def get_species_in_maf(self, maf, out):        
        species = {}
        with open(maf, 'r') as rf:
            for line in rf:
                res=self.get_species_from_maf_line(line)
                if res != -1:
                    species[res] = 1
        list_species = list(species.keys())
        outfile = open(out, "w")
        for s in list_species:
            outfile.write(s + "\n")
        outfile.close()        
    ##   returns the species name from a single line of a maffile
    ##   returns -1 if it is not a line containing a species     
    def get_species_from_maf_line(self,line):
        pattern = re.compile("[a-zA-Z0-9]*\.[a-zA-Z0-9]*")
        l = line.split()
        if(len(l) > 1 and pattern.match(l[1])):        
            return l[1].split(".")[0]
        else:
            return -1       
    def get_strand_from_maf_line(self,line):
        pattern = re.compile("[+-]")
        l = line.split()
        #print(l)
        if(len(l) > 1 and pattern.match(l[4])):        
            return l[4]
        else:
            return -1   
    def get_sequence_from_maf_line(self,line):
        pattern = re.compile("[a-zA-Z0-9]*")
        l = line.split()
        if(len(l) > 1 and pattern.match(l[6])):   
            return l[6]
        else:
            print("AAAAH : ".format(l))
            return -1                 
    def get_scaffold_from_maf_line(self,line):
        pattern = re.compile("[a-zA-Z0-9]*\.[a-zA-Z0-9]*")
        l = line.split()
        if(len(l) > 1 and pattern.match(l[1])):        
            return l[1].split(".")[1]
        else:
            return -1       
    def get_start_from_maf_line(self,line):
        pattern = re.compile("[0-9]")
        l = line.split()
        if(len(l) > 1 and pattern.match(l[2])):        
            return l[2]
        else:
            return -1   
    def get_size_from_maf_line(self,line):
        pattern = re.compile("[0-9]")
        l = line.split()
        if(len(l) > 1 and pattern.match(l[3])):        
            return l[3]
        else:
            return -1   
    def get_similarity(self, a, b):
        return SequenceMatcher(None, a, b).ratio()            
    def get_output_name(self, bfile_name, scaffold, exon_start, upstream_start, output_dir):
        output=output_dir+"/" + bfile_name + "_" + scaffold + "_" + repr(upstream_start) + "_" +  repr(exon_start)
        return output +".bed"      
    ##### FILE TRANSFORMATIONS     
    ##    returns the maf containing only the species inside listfile
    def reduce_maf_via_list(self, listfile, maffile, reduced_maffile):
        species_reduced = {}
        with open(listfile, "r")as lf:
            for species in lf:
                print(species)
                species_reduced[species.strip()] = 1
        lf.close()    
        print(list(species_reduced.keys()))        
        rmf = open(reduced_maffile, "w")                       
        with open(maffile, "r")as mf:
            for line in mf:
                l = self.get_species_from_maf_line(line)                
                if l == -1 or l in species_reduced:
                    print(l)
                    rmf.write(line)
        mf.close()
        rmf.close()    
    ##    returns the bed version of a maf line        
    def is_indel(self, target, query):
        return self.get_size(target) != self.get_size(query)
    def is_gaps_only(self, str_array):
        pattern = re.compile("-+")
        for str in str_array:
            if(not pattern.match(str_array[0])):
                #print(str)
                return False
        return True 
   