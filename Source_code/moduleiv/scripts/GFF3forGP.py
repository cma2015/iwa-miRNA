import re
import sys
from tqdm import tqdm

readfile = sys.argv[1]
outGPs = sys.argv[2]
outPCG = sys.argv[3]

relationDict = {}
reDict = {}
tmpList = []

with open( outPCG, "w") as bed_file:
    with open( readfile ) as gff_file:
        for eli in tqdm(gff_file.readlines()):
            if not eli.startswith("#"):
                eli = eli.strip().split("\t")
                tmp_feat={}
                for x in eli[8].split(';'):
                    tmp_feat.update({x.split('=')[0]:x.split('=')[1]})
                if eli[2] == 'gene':
                    bed_file.write('{}\t{}\t{}\t{}\n'.format(
                        eli[0], re.sub("gene:", "", tmp_feat['ID']), eli[3], eli[4]))
                    relationDict[tmp_feat['ID']] = tmp_feat['ID']
                elif eli[2] == 'mRNA':
                    relationDict[tmp_feat['ID']] = tmp_feat['Parent']
                elif eli[2] == 'CDS':
                    relationDict[tmp_feat['protein_id']] = tmp_feat['Parent']
                    reDict[tmp_feat['protein_id']] = 1

with open( outGPs, "w") as co_file:
    for ii in reDict:
        co_file.write('{}\t{}\n'.format(
            ii, re.sub("gene:", "", relationDict[relationDict[ii]])))
