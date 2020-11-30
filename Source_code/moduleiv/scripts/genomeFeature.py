import re
import os
import sys
features = ["chromosome", "contig", "supercontig", "five_prime_UTR", "exon", "CDS", 
"mRNA", "three_prime_UTR", "lnc_RNA", "ncRNA", "miRNA", "tRNA", "snRNA", "snoRNA", "rRNA",
"SRP_RNA", "pre_miRNA", "RNase_MRP_RNA", "pseudogenic_transcript", "sense_intronic"]

count_list = {}
gene_list = []

readfile = sys.argv[1]
outfile = sys.argv[2]
inputType = sys.argv[3]

if inputType == "GENE":
    with open( outfile, "w") as bed_file:
        with open( readfile ) as gff_file:
            for eli in gff_file.readlines():
                if not eli.startswith("#"):
                    eli = eli.strip().split("\t")
                    if eli[2] not in features:
                        tmp_feat={}
                        for x in eli[8].split(';'):
                            tmp_feat.update({x.split('=')[0]:x.split('=')[1]})
                        if 'ID' in tmp_feat:
                            tmp_id = re.sub('gene:', '',tmp_feat['ID'])
                            tmp_type = "PCG"
                            top_type = "PCG"
                            bed_start = str(int(eli[3])-1)
                            if 'biotype' in tmp_feat:
                                if tmp_feat['biotype'] != "protein_coding":
                                    if 'Parent' in tmp_feat:
                                        tmp_id = re.sub('gene:', '',tmp_feat['Parent'])
                                    else:
                                        tmp_id = re.sub('gene:', '',tmp_feat['ID'])
                                    tmp_type = tmp_feat['biotype']
                                    top_type = eli[2]
                            if tmp_type not in count_list:
                                count_list[tmp_type] = 1
                            else:
                                count_list[tmp_type] += 1
                            if tmp_id not in gene_list:
                                gene_list.append(tmp_id)
                            else:
                                print(tmp_id, tmp_type)
                            if tmp_type != 'miRNA' or tmp_type != 'pre_miRNA':
                                bed_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                                    eli[0], bed_start, eli[4], tmp_id, '.', eli[6], top_type, tmp_type))
elif inputType == "TE":
    with open( outfile, "w") as bed_file:
        with open( readfile ) as gff_file:
            for eli in gff_file.readlines():
                if not eli.startswith("#"):
                    eli = eli.strip().split("\t")
                    tmp_feat={}
                    for x in eli[8].split(';'):
                        tmp_feat.update({x.split('=')[0]:x.split('=')[1]})
                    bed_start = str(int(eli[3])-1)
                    bed_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                        eli[0], bed_start, eli[4], tmp_feat['ID'], '.', eli[6], 'TE', eli[2]))
else:
    pass

# print(count_list, sum(count_list.values()))