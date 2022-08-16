#in zweiter loop detektiert insertions aber ist an und für sich egal

if seq != -1 and tools.is_mostly_non_gaps(seq[query_gap[0]:query_gap[1]]):
    batch_lines.append(
        tools.get_bed_line(species_line, "in", in_number, adapted_dict, scaffold, start, end, strand, indel_length,
                           debug_start, debug_end))
    batch_current_line += 1

    print(target_gaps)
    print(target_line)
    print(query_line)
    print(tools.get_bed_line(query_line, "in", in_number, adapted_dict, scaffold, start, end, strand, indel_length,
                             debug_start, debug_end))




    print(tools.get_bed_line(species_line, "in", in_number, adapted_dict, scaffold, start, end, strand, indel_length,
                             debug_start, debug_end))
    print(species_line)