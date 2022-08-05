from file_parameters import *

#MAX_BLOCKS = 144000
MAX_BLOCKS = 144
batch_number = 0
current_blocks = 0

with open(input_maf, "r") as mf:
    maf_part = open(f"{maf_split}maf_part_{batch_number}.maf", "w")
    for line in mf:
        maf_part.write(line)
        if line == "\n":
            current_blocks += 1
            if current_blocks == MAX_BLOCKS:
                batch_number += 1
                current_blocks = 0
                maf_part.close()
                maf_part = open(f"{maf_split}maf_part_{batch_number}.maf", "w")
                maf_part.write("##maf version=1 scoring=roast\n")
