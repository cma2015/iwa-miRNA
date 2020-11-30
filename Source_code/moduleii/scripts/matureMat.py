import os
import re
import sys

seq_list = {}


with open('mature_miRNAs.fa') as input_flat:
    for eli in input_flat.readlines():
        eli = eli.strip().split('\t')
        if not eli[0].startswith('>') and eli[0] not in seq_list:
            eli[0] = re.sub("U", "T", eli[0])
            seq_list[eli[0]] = eli[0]

path = "{}/4rpmData/".format(sys.argv[1])
files= os.listdir(path)

with open('seq_mat.txt', 'w') as out_report:
    for file in files:
        with open( path + file) as file_input:
            for eli in file_input.readlines():
                eli = eli.strip().split('\t')
                if eli[1] in seq_list:
                    out_report.write('{}\t{}\t{}\n'.format(eli[0], eli[1], eli[2]) )
