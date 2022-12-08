from memory_profiler import profile
from joblib import Parallel, delayed
import functions as tools
import new_functions as tools2
from datetime import datetime

"""
analyses the MAF file and finds InDels. One loop finds t
"""
@profile
def find_indels_from_maf(maf, target, query):
    json_indel_file = maf.replace(".maf", "_indels.json")
    json_block_file = maf.replace(".maf", "_all_blocks.json")
    relevant_json_block_file = maf.replace(".maf", "_relevant_blocks.json")
    indel_out = open(json_indel_file, "w")
    block_out = open(json_block_file, "w")
    relevant_block_out = open(relevant_json_block_file, "w")
    species_lines = []
    target_line = ""
    query_line = ""
    batch_lines_max = 100000
    batch_current_line = 0
    batch_lines_indels = []
    batch_lines_blocks = []
    batch_relevant_lines_block = []
    maf_id = maf.split("_")[2].split(".")[0]
    block_number = 0
    maf_score = ""
    maf_block_line_start = 0
    current_line = 0
    indel_out.write("[")
    block_out.write("[")
    relevant_block_out.write("[")
    with open(maf, "r") as mf:
        for line in mf:
            current_line += 1
            if line.startswith("s"):
                species = tools.get_species_from_maf_line(line)
                if species != -1:
                    if species == query:
                        query_line = line
                    elif species == target:
                        target_line = line
                    else:
                        species_lines.append(line)
            elif line.startswith("a"):
                maf_score = tools.get_maf_score(line)
                maf_block_line_start = current_line
            elif line == "\n":
                chr_name = tools.get_scaffold_from_maf_line(target_line)
                mm10_start = tools.get_start_from_maf_line(target_line)
                mm10_length_no_gaps = len(tools.get_sequence_from_maf_line(target_line))
                block_number += 1
                query_seq = tools.get_sequence_from_maf_line(query_line)
                target_seq = tools.get_sequence_from_maf_line(target_line)
                batch_lines_blocks.append(tools2.make_block_json_entry(
                    block_number, maf_id, maf_score, maf_block_line_start, current_line, target_seq.replace("-", ""), mm10_length_no_gaps, mm10_start, chr_name
                ))
                is_relevant = False
                if query_seq != -1:
                    target_gaps = tools.get_gap_locations(target_seq)
                    for target_gap in target_gaps:
                        query_indel_length_with_gaps = target_gap[1] - target_gap[0]
                        query_indel = query_seq[target_gap[0]:target_gap[1]]
                        insertion_size = tools.get_insertion_size(query_indel)
                        species_agree = []
                        species_disagree = []
                        scaffold, start, end, strand, debug_start, debug_end = tools.calculate_genomic_position(
                            target_line, target_gap[0], target_gap[1], True)
                        if insertion_size > 0:
                            is_relevant = True
                            for species_line in species_lines:
                                species = tools.get_species_from_maf_line(species_line)
                                seq = tools.get_sequence_from_maf_line(species_line)
                                if seq != -1 and tools.is_mostly_non_gaps(seq[target_gap[0]:target_gap[1]]):
                                    species_agree.append(species)
                                else:
                                    species_disagree.append(species)
                            batch_lines_indels.append(tools2.make_indel_json_entry(
                                scaffold, start, end, query_indel_length_with_gaps, insertion_size, "in",  block_number, maf_id,
                                species_agree, species_disagree
                            ))
                            batch_current_line += 1

                    query_gaps = tools.get_gap_locations(query_seq)

                    for query_gap in query_gaps:
                        target_indel = target_seq[query_gap[0]:query_gap[1]]
                        if tools.is_mostly_non_gaps(target_indel):
                            is_relevant = True
                            length_with_gaps = query_gap[1]- query_gap[0]
                            indel_length = tools.get_insertion_size(target_indel)
                            scaffold, start, end, strand, debug_start, debug_end = tools.calculate_genomic_position(
                                target_line, query_gap[0], query_gap[1], False)

                            species_agree = []
                            species_disagree = []
                            for species_line in species_lines:
                                seq = tools.get_sequence_from_maf_line(species_line)
                                species = tools.get_species_from_maf_line(species_line)
                                if tools.is_mostly_gaps(seq[query_gap[0]:query_gap[1]]):
                                    species_agree.append(species)
                                else:
                                    species_disagree.append(species)
                            batch_lines_indels.append(tools2.make_indel_json_entry(
                                scaffold, start, end, length_with_gaps, indel_length, "del",  block_number, maf_id,
                                species_agree, species_disagree
                            ))
                    if is_relevant:
                        batch_relevant_lines_block.append(tools2.make_block_json_entry(
                            block_number, maf_id, maf_score, maf_block_line_start, current_line,
                            target_seq.replace("-", ""), mm10_length_no_gaps, mm10_start, chr_name,
                            len(species_agree) + len(species_disagree) + 2
                        ))
                target_line = ""
                query_line = ""
                species_lines = []

            if batch_current_line >= batch_lines_max:
                tools2.clean_and_write_entries(batch_lines_indels, indel_out)
                tools2.clean_and_write_entries(batch_lines_blocks, block_out)
                tools2.clean_and_write_entries(batch_relevant_lines_block, relevant_block_out)
                batch_lines_indels = []
                batch_lines_blocks = []
                batch_relevant_lines_block = []
                batch_current_line = 0
    tools2.clean_and_write_entries(batch_lines_indels, indel_out, True)
    tools2.clean_and_write_entries(batch_relevant_lines_block, relevant_block_out, True)
    tools2.clean_and_write_entries(batch_lines_blocks, block_out, True)



def loop_find_indels_for_query_from_maf(maf_pieces, maf_part, target, query, bed_out, adapted_dict, number_cores):
    Parallel(n_jobs=number_cores)([delayed(find_indels_from_maf)(maf_part.replace("NUMBER", str(i)),
                                                                      target, query, bed_out.replace("NUMBER", str(i)),
                                                                      adapted_dict) for i in range(maf_pieces)])

def main():
    maf_files = ["maf_part_9_test1000.maf", "maf_part_9_test100.maf"]

    for maf_file in maf_files:
        start = datetime.now()
        find_indels_from_maf(maf_file, "mm10", "nanGal1")
        runtime = datetime.now() - start
        print(f"Runtime: {runtime} Seconds")


if __name__ == '__main__':
    main()
