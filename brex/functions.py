import re

#setup paramets
NON_GAP_TOLERANCE = 0.8
GAP_TOLERANCE = 0.8 #old value: 0.02


def calculate_genomic_position(target, start, end, is_reference_species):
    actual_gap_length = 1
    target_sequence = get_sequence_from_maf_line(target)
    actual_start = calculate_actual_start(target_sequence, start)
    if not is_reference_species:
        actual_gap_length = calculate_genomic_distance(target_sequence, start, end)
    scaffold = get_scaffold_from_maf_line(target)
    start = int(get_start_from_maf_line(target)) + start
    actual_start = int(get_start_from_maf_line(target)) + actual_start
    actual_end = int(actual_start) + actual_gap_length
    size = end - start
    strand = get_strand_from_maf_line(target)
    end = int(start) + size
    return scaffold, actual_start, actual_end, strand, start, end


"""
Calcualate the actual genomic start (subtracting gaps)
"""
def calculate_actual_start(target_sequence, start):
    sequence_to_start = target_sequence[0:start]
    gaps = sequence_to_start.count("-")
    return start - gaps


"""
The actual genomic distance in the reference genome is dependent on the number of gaps 
(not relevant when deletions in the mouse genomes/ inserts in spalax are detected)
"""
def calculate_genomic_distance(target_sequence, indel_start, indel_end):
    indel_sequence = target_sequence[indel_start:indel_end]
    number_of_bases = len(indel_sequence) - indel_sequence.count("-") + 1
    return number_of_bases


"""
Calculates the insertion size, e.g., mouse has gap, but Spalax has an insertion, it only has to be one bp to be relevant:
mouse:  -----
spalax: --A--
still an insertion of size 1
"""
def get_insertion_size(sequence):
    return len(sequence) - sequence.count("-")


def get_bed_line(line, indel, number, adapted_dict, scaffold, start, end, strand, indel_size, debug_start, debug_end):
    species = get_species_from_maf_line(line)
    adapted_label = ""
    is_adapted = False
    for k in adapted_dict.keys():
        if species in adapted_dict[k]:
            is_adapted = True
            if adapted_label == "":
                adapted_label = k
            else:
                adapted_label += f".{k}"
    if is_adapted:
        name = f"{species}.{adapted_label}.{indel}.{number}"
    else:
        name = f"{species}.{indel}.{number}"
    color = get_color(indel, is_adapted)
    return f'{scaffold}\t{start}\t{end}\t{name}\t0\t{strand}\t{debug_start}\t{debug_end}\t{color}\t{indel_size}\n'


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
        if abs(gap.start() - gap.end()) >= 1:
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
