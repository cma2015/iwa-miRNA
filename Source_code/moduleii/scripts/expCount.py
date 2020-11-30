import os
import sys

seq_list = {}
limit_seq_list = {}
limit_count_list = {}
limit_abundance_list = {}
count_list = {}
abundance_list = {}

with open('corr_pre-mirna_seq.txt') as input_flat:
    for eli in input_flat.readlines():
        eli = eli.strip().split('\t')
        if eli[1] not in seq_list:
            seq_list[eli[1]] = [eli[0]]
        else:
            seq_list[eli[1]].append(eli[0])
        if eli[2] == 'complex':
            if eli[1] not in limit_seq_list:
                limit_seq_list[eli[1]] = [eli[0]]
            else:
                limit_seq_list[eli[1]].append(eli[0])

path = "{}/4rpmData/".format(sys.argv[1])
files= os.listdir(path)
for file in files:
    with open( path + file) as file_input:
        for eli in file_input.readlines():
            eli = eli.strip().split('\t')
            ## limit and all
            if eli[0] not in limit_count_list:
                limit_count_list[eli[0]] = {}
            if eli[0] not in count_list:
                count_list[eli[0]] = {}
            # abundance
            if eli[0] not in limit_abundance_list:
                limit_abundance_list[eli[0]] = {}
            if eli[0] not in abundance_list:
                abundance_list[eli[0]] = {}
            # limit
            if eli[1] in limit_seq_list:
                for ii in limit_seq_list[eli[1]]:
                    if ii not in limit_count_list[eli[0]]:
                        limit_count_list[eli[0]][ii] = 1
                        limit_abundance_list[eli[0]][ii] = float(eli[2])
                    else:
                        limit_count_list[eli[0]][ii] += 1
                        limit_abundance_list[eli[0]][ii] += float(eli[2])
            # all
            if eli[1] in seq_list:
                for ii in seq_list[eli[1]]:
                    if ii not in count_list[eli[0]]:
                        count_list[eli[0]][ii] = 1
                        abundance_list[eli[0]][ii] = float(eli[2])
                    else:
                        count_list[eli[0]][ii] += 1
                        abundance_list[eli[0]][ii] += float(eli[2])

with open('corr_pre-mirna_seq_mat.txt', 'w') as out_report:
    for ii in limit_count_list:
        for jj in limit_count_list[ii]:
            tmp_value = limit_count_list[ii][jj]
            out_report.write('count\tlimit\t{}\t{}\t{}\n'.format(ii, jj, tmp_value))
    for ii in count_list:
        for jj in count_list[ii]:
            tmp_value = count_list[ii][jj]
            out_report.write('count\tall\t{}\t{}\t{}\n'.format(ii, jj, tmp_value))
    for ii in limit_abundance_list:
        for jj in limit_abundance_list[ii]:
            tmp_value = limit_abundance_list[ii][jj]
            out_report.write('abundance\tlimit\t{}\t{}\t{}\n'.format(ii, jj, tmp_value))
    for ii in abundance_list:
        for jj in abundance_list[ii]:
            tmp_value = abundance_list[ii][jj]
            out_report.write('abundance\tall\t{}\t{}\t{}\n'.format(ii, jj, tmp_value))
