import re

MIN_INDEL_SIZE = 1
NON_GAP_TOLERANCE = 0.8
GAP_TOLERANCE = 0.8


def calculate_genomic_position(target, s, e, is_reference_species):
    actual_gap_length = 1
    if not is_reference_species:
        actual_gap_length = calculate_genomic_distance(target, s, e)
    scaffold = get_scaffold_from_maf_line(target)
    start = int(get_start_from_maf_line(target)) + s
    strand = get_strand_from_maf_line(target)
    end = int(start) + actual_gap_length
    return scaffold, start, end, strand


def calculate_genomic_distance(target, indel_start, indel_end):
    target_sequence = get_sequence_from_maf_line(target)
    indel_sequence = target_sequence[indel_start:indel_end]
    number_of_bases = len(indel_sequence) - indel_sequence.count("-") + 1
    return number_of_bases


def get_bed_line(line, indel, number, adapted_dict, scaffold, start, end, strand):
    species = get_species_from_maf_line(line)
    adapted_label = ""
    is_adapted = False
    for k in adapted_dict.keys():
        if species in adapted_dict[k]:
            is_adapted = True
            adapted_label = adapted_label + "." + k
    if is_adapted:
        name = species + adapted_label + "." + indel + "." + str(number)
    else:
        name = species + "." + indel + "." + str(number)
    color = get_color(indel, is_adapted)
    return f'{scaffold}\t{start}\t{end}\t{name}\t0\t{strand}\t{start}\t{end}\t{color}\n'


def get_scaffold_from_maf_line(line):
    pattern = re.compile("[a-zA-Z0-9]*\.[a-zA-Z0-9]*")
    l = line.split()
    if len(l) > 1 and pattern.match(l[1]):
        return l[1].split(".")[1]
    else:
        return -1


def get_start_from_maf_line(line):
    pattern = re.compile("[0-9]")
    l = line.split()
    if len(l) > 1 and pattern.match(l[2]):
        return l[2]
    else:
        return -1


def get_strand_from_maf_line(line):
    pattern = re.compile("[+-]")
    l = line.split()
    if len(l) > 1 and pattern.match(l[4]):
        return l[4]
    else:
        return -1


def get_species_from_maf_line(line):
    pattern = re.compile("[a-zA-Z0-9]*\.[a-zA-Z0-9]*")
    l = line.split()
    if len(l) > 1 and pattern.match(l[1]):
        return l[1].split(".")[0]
    return -1


def get_color(indel, is_adapted):
    if indel == "del":
        if is_adapted:
            return "221,106,76"
        return "233,156,93"
    elif indel == "in":
        if is_adapted:
            return "40,151,136"
        return "37,67,80"
    return "133,133,133"


def get_sequence_from_maf_line(line):
    pattern = re.compile("[a-zA-Z0-9]*")
    l = line.split()
    if len(l) > 1 and pattern.match(l[6]):
        return l[6]
    else:
        return -1


def get_gap_locations(seq):
    gap_locations = []
    matches = list(re.finditer('-+', seq))
    for i, gap in enumerate(matches, 1):
        if abs(gap.start() - gap.end()) >= MIN_INDEL_SIZE:
            gap_locations.append([gap.start(), gap.end()])
    return gap_locations


def is_mostly_non_gaps(sequence):
    if len(sequence) == 0:
        print("Warning: length is 0")
        return False
    no_of_gaps = sequence.count('-')
    perc = (len(sequence) - no_of_gaps) / len(sequence)
    if perc < NON_GAP_TOLERANCE:
        return False
    return True


def is_mostly_gaps(sequence):
    if len(sequence) == 0:
        return False
    no_of_gaps = sequence.count('-')
    perc = no_of_gaps / len(sequence)
    if perc < GAP_TOLERANCE:
        return False
    return True
