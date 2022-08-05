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
                if species == target:
                    target_line = line
                else:
                    species_lines.append(line)
        ## block is over.
        if line == "\n":
            if len(query_line) > 0:
                ##If both lines are filled (= the block contains the query species), proceed
                ##if the block contains query species and query species contains indel
                query_seq = self.get_sequence_from_maf_line(query_line)
                target_gaps = self.get_gap_locations(query_seq)
                for i in target_gaps:
                    query_indel = query_seq[i[0]:i[1]]
                    for line in species_lines:
                        seq = self.get_sequence_from_maf_line(line)
                        in_or_out = "in" if current_species in in_group else "out"
                        if seq != -1 and self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                            batch_lines.append(self.get_bed_line(target_line, query_line, "del", i[0], i[1], in_number, in_or_out, adapted_dict))
                            batch_current_line += 1

                    elif self.is_mostly_non_gaps(query_indel):
                        # Todo: correct this double line
                        for line in species_lines:
                            current_species = self.get_species_from_maf_line(line)
                            if current_species in in_group and not current_species in adapted_species:
                                continue
                            else:
                                seq = self.get_sequence_from_maf_line(line)
                                in_or_out = "in" if current_species in in_group else "out"
                                ##if the seq is also non-gaps, but target is, there might be something
                                if seq != -1 and self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                    batch_lines.append(
                                        self.get_bed_line(target_line, line, "in", i[0], i[1], in_number, in_or_out,
                                                          adapted_dict))
                                    batch_current_line += 1
                                ##insertion alert
                                elif seq != -1 and self.is_mostly_gaps(seq[i[0]:i[1]]):
                                    batch_lines.append(
                                        self.get_bed_line(target_line, query_line, "in", i[0], i[1], in_number,
                                                          in_or_out, adapted_dict))
                                    batch_current_line += 1
                    #########################
                    ##      NON GAPS        #
                    ##     IN MOUSE         #
                    #########################
                    target_non_gaps = self.get_non_gap_locations(target_seq)
                    for i in target_non_gaps:
                        gapsize = abs(i[0] - i[1])
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
                                            batch_lines.append(
                                                self.get_bed_line(target_line, line, "del", i[0], i[1], in_number,
                                                                  in_or_out, adapted_dict))
                                            batch_current_line += 1
                            elif self.is_mostly_gaps(query_indel):
                                current_species = self.get_species_from_maf_line(line)
                                if current_species in in_group and not current_species in adapted_species:
                                    continue
                                else:
                                    seq = self.get_sequence_from_maf_line(line)
                                    in_or_out = "in" if current_species in in_group else "out"
                                    if seq != -1 and self.is_mostly_non_gaps(seq[i[0]:i[1]]):
                                        batch_lines.append(
                                            self.get_bed_line(target_line, query_line, "del", i[0], i[1], in_number,
                                                              in_or_out, adapted_dict))
                                        batch_current_line += 1
                                    elif seq != -1 and self.is_mostly_gaps(seq[i[0]:i[1]]):
                                        batch_lines.append(
                                            self.get_bed_line(target_line, line, "del", i[0], i[1], in_number,
                                                              in_or_out, adapted_dict))
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
mf.close()
b_out.close()