import re
import os
import sys

features = ["chromosome", "contig", "supercontig", "five_prime_UTR", "exon", "CDS",
"mRNA", "three_prime_UTR", "lnc_RNA", "ncRNA", "miRNA", "tRNA", "snRNA", "snoRNA", "rRNA",
"SRP_RNA", "pre_miRNA", "RNase_MRP_RNA", "pseudogenic_transcript", "sense_intronic"]

count_list = {}
gene_list = []
outfile = '{0}/GFF3/Annotation.bed'.format(sys.argv[1])
if not os.path.isfile(outfile):
    with open(outfile, "w") as bed_file:
        with open('{0}/GFF3/Annotation.gff3'.format(sys.argv[1])) as gff_file:
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
                            bed_start = str(int(eli[3])-1)
                            if 'biotype' in tmp_feat:
                                if tmp_feat['biotype'] != "protein_coding":
                                    if 'Parent' in tmp_feat:
                                        tmp_id = re.sub('gene:', '',tmp_feat['Parent'])
                                    else:
                                        tmp_id = re.sub('gene:', '',tmp_feat['ID'])
                                    tmp_type = tmp_feat['biotype']
                            if tmp_type not in count_list:
                                count_list[tmp_type] = 1
                            else:
                                count_list[tmp_type] += 1
                            if tmp_id not in gene_list:
                                gene_list.append(tmp_id)
                            else:
                                print(tmp_id, tmp_type)
                            bed_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(eli[0], bed_start, eli[4], tmp_id, tmp_type, eli[6]))
