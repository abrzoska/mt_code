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
    #define a overlap to feather below case
    # A-----A
    # -------    
    GAP_SIMILARITY = 0.8
    NON_GAP_SIMILARITY = 0.8
    INSERT_GAP_TOLERANCE = 0.2
    GAP_TOLERANCE = 0.02
    NON_GAP_TOLERANCE = 0.9
    BAR_SIZE=5
    def __init__(self):
        self.helper = h.Helper()
        logging.basicConfig(filename='brex.log', level=logging.DEBUG)
        logging.info("RegionExtractor instantiated. ")
    ##   find candidates within range for bfile 
    def find_bed_candidates(self, bfile, bfile_name, scaffold, start, end, output_dir):
        gr = self.helper.read_bedfile(bfile)
        print("Scaffold {0}".format(scaffold))
        cand = gr[scaffold,start:end].as_df()
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
    def find_indels_for_query_from_maf(self, maf, target, query, in_group, adapted_species, maf_out, bed_out, b_header):
        m_out = open(maf_out, "w")
        b_out = open(bed_out, "w")
        del_number = 1
        in_number = 1        
        scaffold = ""    
        is_target_in_in_group = target in in_group
        b_out.write(b_header)
        target_line = ""
        query_line = ""
        score_line = ""   
        out_block_lines = []
        in_block_lines = []       
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
                            target_species = self.get_species_from_maf_line(line)
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
                    ##if the block contains query species and query species contains indel
                        ##add in target line -> NO presence is implicit.
                        if is_target_in_in_group:
                            in_block_lines.append(target_line)
                        else: 
                            out_block_lines.append(target_line)   
                        for idx, line in enumerate(out_block_lines):
                            species = self.get_species_from_maf_line(line)
                        if len(scaffold) == 0:
                            scaffold = self.get_scaffold_from_maf_line(target_line) 
                        query_seq = self.get_sequence_from_maf_line(query_line)
                        target_seq = self.get_sequence_from_maf_line(target_line)
                        #########################
                        ##        GAPS          #
                        ##      IN MOUSE        #
                        #########################     
                        target_gaps = self.get_gap_locations(target_seq)
                        #query_gaps = self.get_non_gap_locations(query_seq)
                        #for i in query_gaps:
                        for i in target_gaps:
                            gapsize = abs(i[0]-i[1])
                            if gapsize > self.MIN_INDEL_SIZE:    
                                is_indel_present_in_in_group = False
                                query_indel = query_seq[i[0]:i[1]]
                                ##both are gaps
                                if self.is_mostly_gaps(query_indel):
                                    if target_species in in_group:
                                        is_indel_present_in_in_group = True
                                    else:
                                        for line in in_block_lines:
                                            seq = self.get_sequence_from_maf_line(line)
                                            if self.is_mostly_gaps(seq[i[0]:i[1]]):
                                                is_indel_present_in_in_group = True
                                                break;
                                        if not is_indel_present_in_in_group:
                                            for line in out_block_lines:
                                                seq = self.get_sequence_from_maf_line(line)
                                                current_species = self.get_species_from_maf_line(line)
                                                ##if theyre gap -> shared mutation
                                                if self.is_mostly_gaps(seq[i[0]:i[1]]):
                                                    b_out.write(self.get_bed_line(target_line, line, "del", i[0], i[1], in_number, current_species in adapted_species, False))
                                                ##if theyre non gap -> unique mutation
                                                elif self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                                    b_out.write(self.get_bed_line(target_line, query_line, "del", i[0], i[1], in_number, False, True))
                                elif self.is_mostly_non_gaps(query_indel):
                                    ##target is not gap 
                                    for line in in_block_lines:
                                        seq = self.get_sequence_from_maf_line(line)
                                        if self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                            is_indel_present_in_in_group = True
                                            break;
                                    if not is_indel_present_in_in_group:
                                        for line in out_block_lines:
                                            seq = self.get_sequence_from_maf_line(line)
                                            current_species = self.get_species_from_maf_line(line)
                                            #if self.is_mostly_gaps(seq):
                                            ##if theyre non gap -> shared mutation
                                            if self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                                # print("#######################")
                                                # print("mostly gap t seq?\t{0}\n".format(target_seq[i[0]:i[1]]))
                                                # print("mostly non gap q seq?\t{0}\n".format(query_indel))
                                                # print("mostly non gap out seq?\t{0}\n".format(seq[i[0]:i[1]]))
                                                b_out.write(self.get_bed_line(target_line, line, "in", i[0], i[1], in_number, current_species in adapted_species, False))
                                            ##if theyre gap -> unique mutation
                                            elif self.is_mostly_gaps(seq[i[0]:i[1]]):
                                                b_out.write(self.get_bed_line(target_line, query_line, "in", i[0], i[1], in_number, False, True))
                        #########################
                        ##      NON GAPS        #
                        ##     IN MOUSE         #
                        #########################
                        target_non_gaps = self.get_non_gap_locations(target_seq)
                        for i in target_non_gaps:
                            gapsize = abs(i[0]-i[1])
                            if gapsize > self.MIN_INDEL_SIZE:    
                                is_indel_present_in_in_group = False    
                                query_indel = query_seq[i[0]:i[1]]
                                ##both are non gaps
                                if self.is_mostly_non_gaps(query_indel):
                                    #print("mostly gap query seq?\t{0}\n###".format(query_indel))
                                    if target_species in in_group:
                                        is_indel_present_in_in_group = True
                                    else:
                                        for line in in_block_lines:
                                            seq = self.get_sequence_from_maf_line(line)
                                            if self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                                #print("mostly gap in seq?\t{0}\n###".format(seq[i[0]:i[1]]))
                                                #print("mostly in line gap seq?\n{0}\n###".format(seq))
                                                is_indel_present_in_in_group = True
                                                break;
                                        if not is_indel_present_in_in_group:
                                            for line in out_block_lines:
                                                #print("mostly gap out seq?\t{0}\n###".format(seq[i[0]:i[1]]))
                                                seq = self.get_sequence_from_maf_line(line)
                                                current_species = self.get_species_from_maf_line(line)
                                                ##if theyre non gap -> shared mutation
                                                if self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                                    b_out.write(self.get_bed_line(target_line, line, "del", i[0], i[1], in_number, current_species in adapted_species, False))
                                                ##if theyre gap -> unique mutation
                                                elif self.is_mostly_gaps(seq[i[0]:i[1]]):
                                                    b_out.write(self.get_bed_line(target_line, query_line, "del", i[0], i[1], in_number, False, True))
                                if self.is_mostly_gaps(query_indel):
                                    #print("mostly non gap t seq?\t{0}\n".format(target_seq[i[0]:i[1]]))
                                    #print("mostly non gap q seq?\t{0}\n".format(query_indel))
                                    ##query is gap, target is non gap 
                                    for line in in_block_lines:
                                        seq = self.get_sequence_from_maf_line(line)
                                        if self.is_mostly_gaps(seq[i[0]:i[1]]):
                                            is_indel_present_in_in_group = True
                                            break;
                                    if not is_indel_present_in_in_group:
                                        for line in out_block_lines:
                                            seq = self.get_sequence_from_maf_line(line)
                                            #print("mostly non gap out seq?\t{0}\n".format(seq[i[0]:i[1]]))
                                            current_species = self.get_species_from_maf_line(line)
                                            #if self.is_mostly_gaps(seq):
                                            ##if line gap -> shared mutation
                                            if self.is_mostly_gaps(seq[i[0]:i[1]]):
                                                b_out.write(self.get_bed_line(target_line, line, "del", i[0], i[1], in_number, current_species in adapted_species, False))
                                            ##if line non gap -> unique mutation
                                            elif self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                                b_out.write(self.get_bed_line(target_line, query_line, "del", i[0], i[1], in_number, False, True))   
                            del_number = del_number + 1
                            in_number = in_number + 1
                        #write to maf file                         
                        m_out.write(score_line)
                        m_out.write(target_line)
                        m_out.write(query_line)
                        for line in out_block_lines:
                            seq = self.get_sequence_from_maf_line(line)
                            m_out.write(line)   
                        for line in in_block_lines:
                            seq = self.get_sequence_from_maf_line(line)
                            m_out.write(line)                               
                        m_out.write("\n")
                    scaffold = ""
                    target_line = ""
                    query_line = ""
                    out_block_lines = []
                    in_block_lines = []                 
        mf.close()    
        m_out.close()
        b_out.close()   
    def get_bed_line(self, target, line, indel, s, e, number, is_adapted, is_unique_to_query):
        #print("target \t{0}".format(target))
        #print("line \t{0}".format(line))

        scaffold = self.get_scaffold_from_maf_line(target)
        start = self.get_start_from_maf_line(target)
        size = self.get_size_from_maf_line(target)
        end = int(start)+int(size)
        #print("start {0}, size {1} end {2} --- s {3}, e {4}".format(start, size, end, s, e))
        strand = self.get_strand_from_maf_line(target)
        #start = int(block_start)-self.BAR_SIZE
        #end = int(block_start)+self.BAR_SIZE      

        ##
        #tseq = self.get_sequence_from_maf_line(target)
        #seq = self.get_sequence_from_maf_line(line)
        #print("target seq \t{0}".format(tseq))
        #print("target snip \t{0}".format(tseq[s:e]))
        #print("line seq \t{0}".format(seq))
        #print("blockstart \t{0} + {1}".format(block_start, s))
        #print("blockend \t{0} + {1}".format(block_start, e))
        #print("###")
        ##
        
        name = self.get_species_from_maf_line(line) +"."+ indel + "."+ str(number)
        color = self.get_color(indel, is_adapted, is_unique_to_query)
        return f'{scaffold}\t{start}\t{end}\t{name}\t0\t{strand}\t{start}\t{end}\t{color}\n'
    
    def run(self, start, end, scaffold, maf_in, target_species, query_species, adapted_species, in_group, m_out, b_out, b_header):
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
        self.find_indels_for_query_from_maf(maf_out,target_species,query_species, in_group, adapted_species, m_out, b_out, b_header)
        print("Maf file containing Indels {0}".format(m_out))
        print("Bed file containing Indels {0}".format(b_out))    
    ##   given a query bed and a list of additions .bed files, 
    ##   get a listing of which elements from the query .bed 
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
                        #print("sdlsadklsd -{0}-".format(column_name))
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
            if(abs(gap.start() - gap.end())>self.MIN_INDEL_SIZE):
                gap_locations.append([gap.start(), gap.end()])
        return gap_locations 
    def get_non_gap_locations(self, seq):
        gap_locations = []
        matches = list(re.finditer('[ACGTacgt]+', seq))       
        for i, gap in enumerate(matches, 1):
            #gap_locations.append([gap.start(),gap.end()])
            if(abs(gap.start()- gap.end())>self.MIN_INDEL_SIZE):
                gap_locations.append([gap.start(), gap.end()])            
        #print("seq: {0} -> matches {1}".format(seq,gap_locations))
        return gap_locations         
    def get_size(self, maf_line):
        pattern = re.compile("[0-9]")
        l = maf_line.split()
        if(len(l) > 2 and pattern.match(l[3])):        
            return int(l[3])            
    def get_color(self, indel, is_adapted, is_unique_to_query):
        if is_unique_to_query:
            return "128,128,0"
        if indel == "del":
            if is_adapted:
                return "128,0,128"
            else:
                return "235,64,52"
        elif indel == "in":
            if is_adapted:
                return "73,116,165"
            else:
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
    def is_insertion(self, query, seq):
        complete_gap = "-" * len(seq)
        a = self.get_similarity(seq, complete_gap) #needs to be mostly gap
        b = self.get_similarity(query, complete_gap) #needs to be mostly non gap        
        #print("{0} \n {1}\n->\t {2}".format(query, seq, a))
        return a < self.GAP_SIMILARITY and b < self.INSERT_GAP_TOLERANCE

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
    def is_gaps_only(self, str):
        pattern = re.compile("-+")
        if(not pattern.match(str)):
                return False
        return True 
    def has_no_gaps(self, str):
        pattern = re.compile('[ACGTacgt]+')
        if(not pattern.match(str)):
                return False
        return True 
    def is_mostly_gaps(self, str):
        if self.is_gaps_only(str):
            return True
        no_of_gaps = str.count('-')
        perc = no_of_gaps/len(str)
        if perc < self.GAP_TOLERANCE:
            #print("is_mostly_gaps FALSE? {0} ({1}) -> {2} vs {3}".format(str, no_of_gaps,perc, self.GAP_TOLERANCE))
            return False
        #print("is_mostly_gaps TRUE? {0} ({1}) -> {2} vs {3}".format(str, no_of_gaps,perc, self.GAP_TOLERANCE))
        return True 
    def is_mostly_non_gaps(self,str):
        if self.has_no_gaps(str):
            return True
        no_of_gaps = str.count('-')
        perc = (len(str)-no_of_gaps)/len(str)
        if perc < self.NON_GAP_TOLERANCE:
            #print("is_mostly_non_gaps FALSE? {0} ({1}) -> {2} vs {3}".format(str, no_of_gaps,perc, self.NON_GAP_TOLERANCE))        
            return False
        #print("is_mostly_non_gaps TRUE? {0} ({1}) -> {2} vs {3}".format(str, no_of_gaps,perc, self.NON_GAP_TOLERANCE))
        return True
   