import logging
import requests
import os
import sys
import pyranges as pr
from brex import helper as h
from pathlib import Path
import re

class RegionExtractor:
    MIN_INDEL_SIZE = 4
    def __init__(self):
        self.helper = h.Helper()
        logging.basicConfig(filename='brex.log', level=logging.DEBUG)
        logging.info("RegionExtractor instantiated. ")
    ##
    ##   find candidates within range for bfile 
    ##
    def find_bed_candidates(self, bfile, bfile_name, scaffold, start, end, output_dir):
        gr = self.helper.read_bedfile(bfile)
        ##s = self.helper.remove_quotation_marks(scaffold)
        s=scaffold
        print("Scaffold {0}".format(scaffold))
        cand = gr[s,start:end].as_df()
        output = self.get_output_name(bfile_name, scaffold, start, end, output_dir)
        self.write_bed_cand_to_output(cand, output)
        return output        
    ##
    ##   write results from dataframe to file 
    ##        
    def write_bed_cand_to_output(self, cand, path):
        with open(path, 'w') as f:            
            candAsString = cand.to_string(header = True, index = False)
            f.write(candAsString)  
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
    ##   puts out a bed and a maf (query and target species only) that only  
    ##   include the indels contained in the maf
    ##   TODO indicate in-group
    def find_indels_for_query_from_maf(self, maf, target, query, in_group ,maf_out, bed_out):
        m_out = open(maf_out, "w")
        b_out = open(bed_out, "w")
        indel_number = 1
        scaffold = ""    
        is_target_in_in_group = target in in_group
        #print("is_target_in_in_group {0}".format(is_target_in_in_group))
        b_out.write(f'track name="Indels for {query}" itemRgb="On"\n')
        with open(maf, "r")as mf:
            target_line = ""
            query_line = ""
            score_line = ""
            out_block_lines = []
            in_block_lines = []
            in_group_species_w_same_indel = []
            #header = ""
            #header_written = False
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
                        print(line)
                        print(header)
                ##
                ##   block is over. If both lines are filled (= the block contains the query species), proceed
                ##                
                if line == "\n":         
                    if is_target_in_in_group:
                        in_block_lines.append(target_line)
                    else: 
                        out_block_lines.append(target_line)                
                    ## if the block contains query species and query species contains indel
                    if len(query_line) > 0 and len(target_line) > 0 and self.is_indel(target_line, query_line):
                        if len(scaffold) == 0:
                            scaffold = self.get_scaffold_from_maf_line(target_line)
                        target_gaps = self.get_gap_locations(self.get_sequence_from_maf_line(target_line))
                        for i in target_gaps:
                            if abs(i[0]-i[1]) > self.MIN_INDEL_SIZE:
                                all_seq = []
                                out_block = ""
                                # target_seq = self.get_sequence_from_maf_line(target_line)
                                # print("t: {0}".format(target_seq[i[0]:i[1]]))
                                #out_block = out_block + "t: {0}\n".format(target_seq[i[0]:i[1]])
                                ##check if this gap is present in spalax
                                query_seq = self.get_sequence_from_maf_line(query_line)
                                #print("q: {0}".format(query_seq[i[0]:i[1]]))
                                out_block = out_block + "q: {0}\n".format(query_seq[i[0]:i[1]])
                                all_seq.append(query_seq[i[0]:i[1]])
                                ##check if it present in ingroup
                                for line in in_block_lines:
                                    seq = self.get_sequence_from_maf_line(line)
                                    #print("i: {0}".format(seq[i[0]:i[1]]))
                                    out_block = out_block + "i: {0}\n".format(seq[i[0]:i[1]])
                                    all_seq.append(query_seq[i[0]:i[1]])
                                ##check if present in outgroup
                                for line in out_block_lines:
                                    seq = self.get_sequence_from_maf_line(line)
                                    #print("o: {0}".format(seq[i[0]:i[1]]))
                                    out_block = out_block + "o: {0}\n".format(seq[i[0]:i[1]])
                                    all_seq.append(query_seq[i[0]:i[1]])
                                if not self.is_gaps_only(all_seq):
                                    print("gap location: {0} - {1}".format(i[0],i[1]))
                                    print(out_block)
                                    print("_______")
                        query_gaps = self.get_gap_locations(query_line)
                        #print("query gaps")
                        #print(query_gaps)
                        ##does in-group have same indel
                        ##does out-group have same indel   
                        ##special alert for certain species
                        #write to maf file
                        m_out.write(score_line)
                        m_out.write(target_line)
                        m_out.write(query_line)
                        # print("################################\n\n")
                        # #print(score_line)
                        # print(target_line)
                        # print("#####gaps positions#######")                        
                        # #print(self.get_gap_locations(target_line))                        
                        # print("#####gaps positions#######")                        
                        # print(query_line)
                        # print("#### in_group  ######")
                        # print(in_block_lines)
                        # print("#### out group ######")
                        # print(out_block_lines)
                        # print("################################\n\n")
                        m_out.write("\n")
                        #write to bed file
                        b_out.write(self.get_bed_format_from_maf_lines(target_line,query_line,indel_number,scaffold,query))
                        #count up for naming purposes
                        indel_number = indel_number + 1
                    scaffold = ""
                    target_line = ""
                    query_line = ""
                    out_block_lines = []
                    in_block_lines = []
        mf.close()    
        m_out.close()
        b_out.close()
        
    def is_indel(self, target, query):
        return self.get_size(target) != self.get_size(query)
        
    def is_gaps_only(self, str_array):
        pattern = re.compile("-+")
        for str in str_array:
            if(not pattern.match(str_array[0])):
                #print(str)
                return False
        return True
                
    ## GETTER    
    def get_gap_locations(self, seq):
        gap_locations = []
        matches = list(re.finditer('-+', seq))
        for i, gap in enumerate(matches, 1):
            #gap_locations.append([gap.start(),gap.end()-1])
            gap_locations.append([gap.start(),gap.end()])
        return gap_locations            
    def get_in_or_del(self, targetsize, querysize):
        if targetsize > querysize:
            return "del"
        elif targetsize < querysize:
            return "in"
        else:
            return "m"            
    def get_size(self, maf_line):
        pattern = re.compile("[0-9]")
        l = maf_line.split()
        if(len(l) > 2 and pattern.match(l[3])):        
            return int(l[3])            
    def get_color(self, indel):
        if indel == "del":
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
    def get_sequence_from_maf_line(self,line):
        pattern = re.compile("[a-zA-Z0-9]*")
        l = line.split()
        if(len(l) > 1 and pattern.match(l[6])):   
            return l[6]
        else:
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
    def get_output_name(self, bfile_name, scaffold, exon_start, upstream_start, output_dir):
        output=output_dir+"/" + bfile_name + "_" + scaffold + "_" + repr(upstream_start) + "_" +  repr(exon_start)
        return output +".bed"      
    ##### FILE TRANSFORMATIONS     
    ##
    ##    returns the maf containing only the species inside listfile
    ##
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
                    rmf.write(line)
        mf.close()
        rmf.close()    
    ##
    ##    returns the bed version of a maf line
    ##        
    def get_bed_format_from_maf_lines(self, target, query, number, scaffold, query_species):
        t_size = self.get_size(target)
        q_size = self.get_size(query)        
        start = self.get_start_from_maf_line(target)
        end = int(start) + t_size
        indel = self.get_in_or_del(t_size, q_size)
        name = query_species + "."+ str(number) +"."+ indel 
        #print(indel)
        color = self.get_color(indel)
        #print(color)
        return f'{scaffold}\t{start}\t{end}\t{name}\t0\t.\t{start}\t{end}\t{color}\n'     
##        