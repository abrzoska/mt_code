import logging
import requests
import os
import sys
import bamread
import pysam
import pyranges as pr
from brex import helper as h
from pathlib import Path
import numpy as np

class RegionExtractor:
    def __init__(self):
        self.helper = h.Helper()
        logging.basicConfig(filename='brex.log', level=logging.DEBUG)
        logging.info("RegionExtractor instantiated. ")

    def get_output_name(self, scaffold, exon_start, upstream_start, output_dir):
        output=output_dir+"/" + scaffold + "_" + repr(upstream_start) + "_" +  repr(exon_start)
        return output +".bed"
            
    def find_bed_candidates(self, bfile,scaffold, start, end,output_dir):
        gr = self.helper.read_bedfile(bfile)
        s = self.helper.remove_quotation_marks(scaffold)
        cand = gr[s,self.start:self.end].as_df()
        output = self.get_output_name(scaffold, start, end, output_dir)
        self.write_bed_cand_to_output(cand, output)
        
    def write_bed_cand_to_output(self, cand, path):
        with open(path, 'w') as f:
            #header="TODO"
            #candAsString = cand.to_string(header = True, index = False)
            candAsString = cand.to_string(header = True, index = False)
            print(candAsString)
            f.write(candAsString)
            
    #track name=HbVar type=bedDetail description="HbVar custom track"
    def annotate_bed_with_overlap(self, query_bed, target_beds):
        target_bed_files = []
        
        #todo extract from shit
        scaffold = "chr8"
        for file in target_beds:
            target_bed_files.append(self.helper.read_bedfile(file))            
        query_file = open(query_bed, "r")    
#        cand = gr[s,self.start:self.end].as_df()        
        for line in query_file:
            if(not line.startswith("tr")):                
                [start,end] = self.get_range_from_row(line)
                frames = []
                #i = 0
                for file in target_bed_files:                    
                    cand = file[scaffold,start:end].as_df()
                    if not cand.empty:
                        #print(str(start) + "-" + str(end)) 
                        #print(cand)
                        #print("#####################")
                        frames.append(cand)
                #TODO
                if len(frames) > 0:
                #if len(frames) >= 2:
                    print("#####################")
                    print(line[0:40])
                    #print("______")                
                    #print(len(frames))
                    #print(frames[0])
                    if len(frames) >= 2:
                        #print(frames[1])
                        #TODO 
                        result= frames[0].append(frames[1])
                    else:
                        result=frames[0]
                    print(result)
                    print("#####################")

                #i = i + len(frames)                      
                #print(str(start) +" "+ str(end))
        query_file.close()
        
    def get_range_from_row(self,row):
        #print("rrrr "+row)
        ar = row.split()
        if(len(ar)>2 and not ar[0].startswith("tr")):
            return [int(ar[1]),int(ar[2])]
            #return [ar[2],ar[3]]
        #return [1,1]