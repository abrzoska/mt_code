from enum import Enum
from pathlib import Path
import pyranges as pr
import os
import logging
import sys

class Helper:

    def exit_w_error(self, e):
        logging.error(e)
        sys.exit(e)

    def read_bedfile(self, file):
        self.check_bfile(file)
        return pr.read_bed(file)

    def check_bfile(self, file):
        file_to_open=Path(file)
        print(file_to_open)
        if not os.path.exists(file_to_open):
            self.exit_w_error("Could not find file. Path does not exist: " +file)

    def cut_down_bed_file(self,old_file_name, new_file_name, size):
        file=open(old_file_name, 'r')
        new_file=open(new_file_name, 'w')
        for line in file:
            line = line.split()
            new_line = '\t'.join(line[0:size])+ "\n"
            new_file.write(new_line)
        file.close()
        new_file.close()

    def remove_quotation_marks(self, str):
        return str.replace("'", "")



