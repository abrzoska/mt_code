import logging
import os
import pandas as pd
from brex import helper as h
from pathlib import Path
import re
import subprocess
from difflib import SequenceMatcher
from time import perf_counter
from memory_profiler import profile
from joblib import Parallel, delayed
import functions as tools

class RegionExtractor:
    #define an overlap to feather below case
    # A-----A
    # -------    
    GAP_SIMILARITY = 0.8
    INSERT_GAP_TOLERANCE = 0.2
    GAP_TOLERANCE = 0.8 #old value: 0.02
    NON_GAP_TOLERANCE = 0.9
    def __init__(self, run_name, folder, min_indel_size):
        self.helper = h.Helper()
        name_log = "RE_{0}.log".format(run_name)
        logging.basicConfig(filename=name_log, level=logging.DEBUG)
        logging.info("RegionExtractor instantiated")
        self.TMP_FOLDER = folder + "tmp"
        self.MIN_INDEL_SIZE = min_indel_size

    ##   find candidates within range for bfile 
    def find_bed_candidates(self, bfile, bfile_name, scaffold, start, end, output_dir):
        gr = self.helper.read_bedfile(bfile)
        print("Scaffold {0}".format(scaffold))
        logging.info("Scaffold {0}".format(scaffold))
        cand = gr[scaffold,start:end].as_df()
        output = self.get_output_name(bfile_name, scaffold, start, end, output_dir)
        self.write_bed_cand_to_output(cand, output)
        return output
    ##   write results from dataframe to file
    def write_bed_cand_to_output(self, cand, path):
        with open(path, 'w') as f:
            candAsString = cand.to_string(header = True, index = False)
            f.write(candAsString)

    def find_indels_from_maf(self, maf, target, query, bed_out, adapted_dict):
        b_out = open(bed_out, "w")
        species_lines = []
        target_line = ""
        query_line = ""
        batch_lines_max = 100000
        batch_current_line = 0
        batch_lines = []
        in_number = 0
        with open(maf, "r") as mf:
            for line in mf:
                if line.startswith("s"):
                    species = self.get_species_from_maf_line(line)
                    if species != -1:
                        if species == query:
                            query_line = line
                        elif species == target:
                            target_line = line
                        else:
                            species_lines.append(line)
                elif line == "\n":
                    query_seq = tools.get_sequence_from_maf_line(query_line)
                    target_seq = tools.get_sequence_from_maf_line(target_line)

                    if target_seq != -1 and query_seq != -1:
                        target_gaps = tools.get_gap_locations(target_seq)
                        for target_gap in target_gaps:
                            query_indel = query_seq[target_gap[0]:target_gap[1]]
                            if tools.get_insertion_size(query_indel) > 0:
                                indel_length = target_gap[1] - target_gap[0]
                                scaffold, start, end, strand, debug_start, debug_end = tools.calculate_genomic_position(target_line, target_gap[0], target_gap[1], True)
                                batch_lines.append(tools.get_bed_line(query_line, "in", in_number, adapted_dict, scaffold, start, end, strand, indel_length, debug_start, debug_end))
                                batch_current_line += 1
                                for species_line in species_lines:
                                    seq = tools.get_sequence_from_maf_line(species_line)
                                    if seq != -1 and tools.is_mostly_non_gaps(seq[target_gap[0]:target_gap[1]]):
                                        batch_lines.append(tools.get_bed_line(species_line, "in", in_number, adapted_dict, scaffold, start, end, strand, indel_length, debug_start, debug_end))
                                        batch_current_line += 1
                                in_number += 1
                        query_gaps = tools.get_gap_locations(query_seq)

                        for query_gap in query_gaps:
                            target_indel = target_seq[query_gap[0]:query_gap[1]]
                            if tools.is_mostly_non_gaps(target_indel):
                                indel_length = query_gap[1] - query_gap[0]
                                scaffold, start, end, strand, debug_start, debug_end = tools.calculate_genomic_position(target_line, query_gap[0], query_gap[1], False)
                                batch_lines.append(tools.get_bed_line(query_line, "del", in_number, adapted_dict, scaffold, start, end, strand, indel_length, debug_start, debug_end))
                                batch_current_line += 1
                                for species_line in species_lines:
                                    seq = tools.get_sequence_from_maf_line(species_line)
                                    if seq != -1 and tools.is_mostly_gaps(seq[query_gap[0]:query_gap[1]]):
                                        batch_lines.append(tools.get_bed_line(species_line, "del", in_number, adapted_dict, scaffold, start, end, strand, indel_length, debug_start, debug_end))
                                        batch_current_line += 1
                                in_number += 1
                    target_line = ""
                    query_line = ""
                    species_lines = []
                if batch_current_line >= batch_lines_max:
                    b_out.writelines(batch_lines)
                    batch_lines = []
                    batch_current_line = 0

        b_out.writelines(batch_lines)
        mf.close()
        b_out.close()

    ##   puts out a bed and a maf (query and target species only) that only
    ##   include the indels contained in the maf
    ##   TODO indicate in-group
    @profile
    def find_indels_for_query_from_maf(self, maf, target, query, in_group, bed_out, adapted_dict):
        b_out = open(bed_out, "w")
        del_number = 1
        in_number = 1
        adapted_species = adapted_dict.values()
        target_line = ""
        query_line = ""
        species_lines = []
        batch_lines_max = 100000
        batch_current_line = 0
        batch_lines = []
        with open(maf, "r") as mf:
            for line in mf:
                if line.startswith("s"):
                    species = self.get_species_from_maf_line(line)
                    if species != -1:
                        if species == query:
                            query_line = line
                        elif species == target:
                            target_line = line
                            #Todo: schauen ob man nächste zeile überhaupt braucht
                            species_lines.append(target_line)
                        else:
                            species_lines.append(line)
                ## block is over.
                if line == "\n":
                    if len(query_line) > 0 and len(target_line) > 0 and self.is_indel(target_line, query_line):
                    ##If both lines are filled (= the block contains the query species), proceed
                    ##if the block contains query species and query species contains indel
                        query_seq = self.get_sequence_from_maf_line(query_line)
                        target_seq = self.get_sequence_from_maf_line(target_line)
                        if query_seq != -1 and target_seq != -1:
                            ##add in target line -> NO presence is implicit.
                            #########################
                            ##        GAPS          #
                            ##      IN MOUSE        #
                            #########################
                            target_gaps = self.get_gap_locations(target_seq)
                            for i in target_gaps:
                                query_indel = query_seq[i[0]:i[1]]
                                    ##query is gap
                                if self.is_mostly_gaps(query_indel):
                                    #Todo: correct this double line
                                    for line in species_lines:
                                        current_species = self.get_species_from_maf_line(line)
                                        if current_species in in_group and not current_species in adapted_species:
                                            continue
                                        else:
                                            #print("in: mouse gaps, spalax gaps")
                                            seq = self.get_sequence_from_maf_line(line)
                                            in_or_out = "in" if current_species in in_group else "out"
                                            ##if the seq is also gaps, no additional info is gained
                                            if seq != -1 and self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                                batch_lines.append(self.get_bed_line(target_line, query_line, "del", i[0], i[1], in_number, in_or_out, adapted_dict))
                                                batch_current_line += 1
                                elif self.is_mostly_non_gaps(query_indel):
                                    #Todo: correct this double line
                                    for line in species_lines:
                                        current_species = self.get_species_from_maf_line(line)
                                        if current_species in in_group and not current_species in adapted_species:
                                            continue
                                        else:
                                            seq = self.get_sequence_from_maf_line(line)
                                            in_or_out = "in" if current_species in in_group else "out"
                                            ##if the seq is also non-gaps, but target is, there might be something
                                            if seq != -1 and self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                                batch_lines.append(self.get_bed_line(target_line, line, "in", i[0], i[1], in_number, in_or_out, adapted_dict))
                                                batch_current_line += 1
                                            ##insertion alert
                                            elif seq != -1 and self.is_mostly_gaps(seq[i[0]:i[1]]):
                                                batch_lines.append(self.get_bed_line(target_line, query_line, "in", i[0], i[1], in_number, in_or_out, adapted_dict))
                                                batch_current_line += 1
                            #########################
                            ##      NON GAPS        #
                            ##     IN MOUSE         #
                            #########################
                            target_non_gaps = self.get_non_gap_locations(target_seq)
                            for i in target_non_gaps:
                                gapsize = abs(i[0]-i[1])
                                if gapsize > self.MIN_INDEL_SIZE:
                                    query_indel = query_seq[i[0]:i[1]]
                                    ##both target and query are non gaps
                                    if self.is_mostly_non_gaps(query_indel):
                                        # Todo: correct this double line
                                        for line in species_lines:
                                            current_species = self.get_species_from_maf_line(line)
                                            if current_species in in_group and not current_species in adapted_species:
                                                continue
                                            else:
                                                seq = self.get_sequence_from_maf_line(line)
                                                in_or_out = "in" if current_species in in_group else "out"
                                                ##if seq also non gaps, no additional info is gained
                                                if seq != -1 and self.is_mostly_gaps(seq[i[0]:i[1]]):
                                                    batch_lines.append(self.get_bed_line(target_line, line, "del", i[0], i[1], in_number, in_or_out, adapted_dict))
                                                    batch_current_line += 1
                                    elif self.is_mostly_gaps(query_indel):
                                        current_species = self.get_species_from_maf_line(line)
                                        if current_species in in_group and not current_species in adapted_species:
                                            continue
                                        else:
                                            seq = self.get_sequence_from_maf_line(line)
                                            in_or_out = "in" if current_species in in_group else "out"
                                            if seq != -1 and self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                                batch_lines.append(self.get_bed_line(target_line, query_line, "del", i[0], i[1], in_number, in_or_out, adapted_dict))
                                                batch_current_line += 1
                                            elif seq != -1 and self.is_mostly_gaps(seq[i[0]:i[1]]):
                                                batch_lines.append(self.get_bed_line(target_line, line, "del", i[0], i[1], in_number, in_or_out, adapted_dict))
                                                batch_current_line += 1



                            ##TODO count numbers up only if present.
                                del_number = del_number + 1
                                in_number = in_number + 1
                        if batch_current_line >= batch_lines_max:
                            b_out.writelines(batch_lines)
                            batch_lines = []
                            batch_current_line = 0
                    ############################################################
                    target_line = ""
                    query_line = ""
                    species_lines = []
        b_out.writelines(batch_lines)
        mf.close()
        b_out.close()

    def loop_find_indels_for_query_from_maf(self, maf_pieces, maf_part, target, query, in_group, bed_out, adapted_dict, number_cores):
        Parallel(n_jobs=number_cores)([delayed(self.find_indels_for_query_from_maf)(maf_part.replace("NUMBER", str(i)),
            target, query, in_group, bed_out.replace("NUMBER", str(i)), adapted_dict) for i in range(maf_pieces)])

    ##create bed line from detected INDEL for output incl in/out and adapted value
    def get_bed_line(self, target, line, indel, s, e, number, in_or_out_group, adapted_dict):
        scaffold = self.get_scaffold_from_maf_line(target)
        start = int(self.get_start_from_maf_line(target)) + s
        species = self.get_species_from_maf_line(line)
        end = int(start) + e
        adapted_label = ""
        is_adapted = False
        for k in adapted_dict.keys():
            if species in adapted_dict[k]:
                is_adapted = True
                adapted_label = adapted_label + "." + k
        strand = self.get_strand_from_maf_line(target)
        if(is_adapted):
            name = species+"."+ in_or_out_group + adapted_label +"."+ indel + "."+ str(number)
        else:
            name = species+"."+ in_or_out_group+"."+ indel + "."+ str(number)
        color = self.get_color(indel, is_adapted)
        return f'{scaffold}\t{start}\t{end}\t{name}\t0\t{strand}\t{start}\t{end}\t{color}\n'
    ##run analysis for folder of bed files (output to same folder)
    @profile
    def run_analysis(self, folder, adapted_labels, query_species, run_name):
    ##iterate over files in folder     
        list_of_final_elements = []
        header = "filename"+","+"query_species_only_count"+","+"indel_non_adapted_included_count"
        for x in adapted_labels:
            header = header + "," + x
        #print(header)
        ##adapted labels
        stats = []
        for filename in os.listdir(folder):
            analysis_output_file = f"{folder}{filename}_analysis.csv"
            f = os.path.join(folder, filename)
            query_species_only_names = []
            query_species_only_count = 0
            indel_non_adapted_included_names = []
            indel_non_adapted_included_count = 0
            adapted_species_dicts = {}
            gene_line_start_counter = perf_counter()
            if os.path.isfile(f) and f.endswith(".bed") and not f.endswith("_analysis.bed"):
                columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand,", "ThickStart", "ThickEnd", "ItemRGB"]
                df = pd.read_csv(f, sep="\t", skiprows=1, header=None, names=columns)
                df = df.drop_duplicates(subset=['Name'])
                ##analysis
                max_number = int(list(df.tail(1)['Name'])[0].split(".")[-1])
                if max_number > 1:
                    for i in range(1, max_number):
                        line_ending_in = "."+str(i)
                        lines_with_this_end_number = df[df.Name.str.endswith(line_ending_in)]

                        if len(lines_with_this_end_number) == 1:
                            only_query_species_with_this_end_number = lines_with_this_end_number.iloc[0]
                            if query_species in only_query_species_with_this_end_number.Name:
                                query_species_only_names.append(only_query_species_with_this_end_number.Name)
                                query_species_only_count = query_species_only_count + 1
                        ##InDel ist in Spalax und ...
                        else:
                            names_list = lines_with_this_end_number['Name'].tolist()
                            if any(query_species in s for s in names_list):# whatever
                                ## 2) InDel ist in Spalax UND einer Spezies aus der Outgroup (egal ob long lived/adaptiert)
                                ##PUT ALL IN
                                indel_non_adapted_included_names = indel_non_adapted_included_names + names_list
                                indel_non_adapted_included_count = indel_non_adapted_included_count + len(names_list)
                            ## 3) InDel ist in Spalax UND einer Hypoxie adaptierten Spezies, egal ob Ingroup oder Outgroup
                            ## 4) InDel ist in Spalax UND einer langlebigen Spezies , egal ob Ingroup oder   Outgroup
                                for adapted_label in adapted_labels:
                                    adapted_species = [s for s in names_list if "."+adapted_label+"." in s]
                                    if adapted_label in adapted_species_dicts.keys():
                                        adapted_species_dicts[adapted_label] = adapted_species_dicts[adapted_label] + adapted_species
                                    else:
                                        adapted_species_dicts[adapted_label] = adapted_species
                print(f"Runtime for loop: {perf_counter() - gene_line_start_counter}")
                list_of_final_elements = list_of_final_elements + query_species_only_names + indel_non_adapted_included_names
                row = filename + "," +str(query_species_only_count)+ "," +str(indel_non_adapted_included_count)
                for x in adapted_labels:
                    if x in adapted_species_dicts.keys():
                        list_of_final_elements = list_of_final_elements + adapted_species_dicts[x]
                        row = row+","+str(len(adapted_species_dicts[x]))
                stats.append(row)
                final_elements = df[df['Name'].isin(list_of_final_elements)]
                final_elements.to_csv(analysis_output_file, index=False)
                print(f"Runtime total: {perf_counter() - gene_line_start_counter}")
        analysis_file = open(folder + run_name + "_analysis.csv", "w")
        analysis_file.write(header + "\n")
        analysis_file.writelines(stats)
        analysis_file.close()
    ##convert csv with rgb values to bed
    ##TODO remove files.
    def convert_csv_to_bed(self, headerline, analysis_output_file,analysis_output_file_bed):
        print(analysis_output_file)
        if os.path.isfile(analysis_output_file):
            #analysis_output_file_bed = analysis_output_file.with_suffix('.bed')
            a_file = open(analysis_output_file, "r")
            lines = a_file.readlines()
            a_file.close()
            new_file = open(analysis_output_file_bed, "w")
            new_file.write(headerline+"\n")
            for line in lines:
                if not line.strip("\n").startswith("Chromosome"):
                    ##this kills the rgb values -- workaround for now
                    #print("####")
                    #print(line)
                    line = line.strip("\n")
                    if line.endswith('"'):
                        ##contains rgb
                        line_rgb = re.findall('"([^"]*)"', line)[0]
                        line = ('\t'.join(line.split(",")[0:-3])) + '\t' + line_rgb + '\n'
                        new_file.write(line)
                    else:
                        line = line.replace(',','\t')
                        new_file.write(line)
            new_file.close()
        else:
            print("File does not exist: "+ analysis_output_file)
    ##from maf and bed find indels within upstream region
    def start_create_bed_and_maf_for_list_of_genes(self, output_folder, gene_list, bio_mart_csv, upstream_region, input_maf, target_species, query_species, in_group, adapted_dict):
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        gene_file = open(gene_list, "r")
        df = pd.read_csv(bio_mart_csv)
        for line in gene_file:
            gene_line_start_counter = perf_counter()
            line = line.strip()
            print("Running {0}".format(line))
            #Todo: if slow, then use dic instead of df
            gene = df[df['Gene_stable_ID'] == line]
            if not gene.empty:#len(gene.index) > 0:
                #print(gene)
                #Todo remove duplicate variable instantiation
                column = gene["Chromosome_scaffold_name"]
                column = gene["Transcript_start"]
                latest_transcription_start = column.max()
                earliest_transcription_start = column.min()
                upstream_cutoff = earliest_transcription_start-upstream_region
                print("Upstream cutoff {0} - Latest Transcription start {1}".format(upstream_cutoff, latest_transcription_start))
                s = gene.loc[gene.index[0], 'Chromosome_scaffold_name']
                if s.isnumeric():
                    scaffold = "chr{0}".format(s)
                else:
                    scaffold = s
                print("scaffold {0}".format(scaffold))
                b_header = f'track name="Indels for {query_species} {line} ({scaffold}:{upstream_cutoff}-{latest_transcription_start})" itemRgb="On"\n'
                [output_bed, output_maf] = self.get_output_files_name(line, output_folder)
                print("output_bed {0}".format(output_bed))
                self.create_bed_and_maf_for_alignment(upstream_cutoff, latest_transcription_start, scaffold, input_maf, target_species, query_species, in_group, output_maf, output_bed, b_header, adapted_dict)
                print("...done")
            else:
                #Todo log which genes did not make it
                print("Gene {0} not found in BioMart file, skipping".format(line))
            gene_line_stop_counter = perf_counter()
            time = gene_line_stop_counter-gene_line_start_counter
            print('Elapsed time in seconds:", %.4f' % time)
        gene_file.close()
    def get_output_files_name(self, line, output_folder):
        name_bed = "{0}/{1}-indels.bed".format(output_folder,line)
        name_maf = "{0}/{1}-indels.maf".format(output_folder,line)
        return [name_bed, name_maf]
    ##! if maf extract does not exist it will be created which can take multiple hours
    ## in case of unstable connection to server run via tmux (see start_genome_session.sh)
    def create_bed_and_maf_for_alignment(self, start, end, scaffold, maf_in, target_species, query_species, in_group, m_out, b_out, b_header, adapted_dict):
        ## check all files + log
        complete_path = os.getcwd()+self.TMP_FOLDER
        print("\tLocation of tmp folder: {0}".format(complete_path))
        logging.info("\tLocation of tmp folder: {0}".format(complete_path))
        if not Path(complete_path).is_dir():
            subprocess.call(['mkdir', complete_path])
            print("\tCreated Folder {0}.".format(self.TMP_FOLDER))
            logging.info("\tCreated Folder {0}.".format(self.TMP_FOLDER))
        else:
            print("\t tmp folder {0} exists.".format(complete_path))
            logging.info("\ttmp folder {0} exists.".format(complete_path))
        ## create bed from start, end
        name = "{0}_{1}_{2}".format(scaffold,start,end)
        bed = complete_path+"/"+name+".bed"
        maf_out = complete_path+"/"+name+".maf"
        if not Path(bed).is_file():
            print("\tWriting bed to file .. {0}".format(bed))
            logging.info("\tWriting bed to file .. {0}".format(bed))
            with open(bed, 'a') as bedfile:
                bedfile.write("{0}\t{1}\t{2}\t{3}".format(scaffold,start,end,name))
            bedfile.close()
        else:
            print("\tbed file {0} exists.".format(bed))
            logging.info("\tbed file {0} exists.".format(bed))
        indel_start_counter = perf_counter()
        self.find_indels_for_query_from_maf(maf_out,target_species,query_species, in_group, m_out, b_out, b_header,adapted_dict)
        indel_stop_counter = perf_counter()
        time = indel_stop_counter-indel_start_counter
        #print("\tElapsed time in seconds: (find_indels_for_query_from_maf)", time.isoformat(timespec='seconds'))
        print('\tElapsed time: %.4f s' % time)
        logging.info('\tElapsed time: %.4f s' % time)
        print("\tMaf file containing Indels {0}".format(m_out))
        logging.info("\tMaf file containing Indels {0}".format(m_out))
        print("\tBed file containing Indels {0}".format(b_out))
        logging.info("\tBed file containing Indels {0}".format(b_out))
    ##   given a query bed and a list of additions .bed files, 
    ##   get a listing of which elements from the query .bed
    ##   have overlaps with elements from the target beds
    ##   track name=HbVar type=bedDetail description="HbVar custom track"
    ##TODO
    def filter_bed_with_overlap(self, query_bed, target_beds):
        target_bed_files = []
        for file in target_beds:
            tmp_file_name =file+"-tmp"
            self.helper.cut_down_bed_file(file,tmp_file_name , 4)
            target_bed_files.append(self.helper.read_bedfile(tmp_file_name))
        print(len(target_bed_files))
        query_file = open(query_bed, "r")
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
                        #print(file)
                        cand = file[scaffold,start:end].as_df()
                        print(cand)
                        column_name = cand.columns
                        #cand['Notes'] = cand. + "_" + name
        query_file.close()

    def get_gap_locations(self, seq):
        gap_locations = []
        matches = list(re.finditer('-+', seq))
        for i, gap in enumerate(matches, 1):
            if abs(gap.start() - gap.end()) >= self.MIN_INDEL_SIZE:
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
        if len(l) > 2 and pattern.match(l[3]):
            return int(l[3])

    def get_color(self, indel, is_adapted):
        if indel == "del":
            if is_adapted:
                return "221,106,76"
            return "233,156,93"
        elif indel == "in":
            if is_adapted:
                return "40,151,136"
            return "37,67,80"
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
        if len(l) > 1 and pattern.match(l[1]):
            return l[1].split(".")[0]
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
        #print("OOOHH : ".format(line))
        pattern = re.compile("[a-zA-Z0-9]*")
        l = line.split()
        if(len(l) > 1 and pattern.match(l[6])):
            return l[6]
        else:
            #print("AAAAH : ".format(line))
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

    #Todo: remove or improve this method: this is true as long as first position is "-"
    def is_gaps_only(self, string):
        pattern = re.compile("-+")
        if pattern.match(string):
            return True
        return False

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
