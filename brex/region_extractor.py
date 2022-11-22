from memory_profiler import profile
from joblib import Parallel, delayed
import functions as tools

"""
analyses the MAF file and finds InDels. One loop finds t
"""
@profile
def find_indels_from_maf(maf, target, query, bed_out, adapted_dict):
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
                species = tools.get_species_from_maf_line(line)
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
                            scaffold, start, end, strand, debug_start, debug_end = tools.calculate_genomic_position(
                                target_line, target_gap[0], target_gap[1], True)
                            batch_lines.append(
                                tools.get_bed_line(query_line, "in", in_number, adapted_dict, scaffold, start, end,
                                                   strand, indel_length, debug_start, debug_end))
                            batch_current_line += 1
                            for species_line in species_lines:
                                seq = tools.get_sequence_from_maf_line(species_line)
                                if seq != -1 and tools.is_mostly_non_gaps(seq[target_gap[0]:target_gap[1]]):
                                    batch_lines.append(
                                        tools.get_bed_line(species_line, "in", in_number, adapted_dict, scaffold, start,
                                                           end, strand, indel_length, debug_start, debug_end))
                                    batch_current_line += 1
                            in_number += 1
                    query_gaps = tools.get_gap_locations(query_seq)

                    for query_gap in query_gaps:
                        target_indel = target_seq[query_gap[0]:query_gap[1]]
                        if tools.is_mostly_non_gaps(target_indel):
                            indel_length = query_gap[1] - query_gap[0]
                            scaffold, start, end, strand, debug_start, debug_end = tools.calculate_genomic_position(
                                target_line, query_gap[0], query_gap[1], False)
                            batch_lines.append(
                                tools.get_bed_line(query_line, "del", in_number, adapted_dict, scaffold, start, end,
                                                   strand, indel_length, debug_start, debug_end))
                            batch_current_line += 1
                            for species_line in species_lines:
                                seq = tools.get_sequence_from_maf_line(species_line)
                                if seq != -1 and tools.is_mostly_gaps(seq[query_gap[0]:query_gap[1]]):
                                    batch_lines.append(
                                        tools.get_bed_line(species_line, "del", in_number, adapted_dict, scaffold,
                                                           start, end, strand, indel_length, debug_start, debug_end))
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


def loop_find_indels_for_query_from_maf(maf_pieces, maf_part, target, query, bed_out, adapted_dict, number_cores):
    Parallel(n_jobs=number_cores)([delayed(find_indels_from_maf)(maf_part.replace("NUMBER", str(i)),
                                                                      target, query, bed_out.replace("NUMBER", str(i)),
                                                                      adapted_dict) for i in range(maf_pieces)])
