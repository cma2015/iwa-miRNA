#!/usr/bin/python
#-*-coding:utf-8 -*-

import sys
from Bio import SeqIO

beforeSeq={}
beforeModSeq={}
afterSeq={}
seqStrand={}
bed_file = open(sys.argv[1], "r")
for i in bed_file.readlines():
    i = i.strip()
    lineStrand = i.split("\t")
    if lineStrand[5]=="-":
        seqStrand[lineStrand[3]]="-"
    
bed_file.close()

for name_infor in SeqIO.parse(sys.argv[2], 'fasta'):
    extract_name = name_infor.id
    beforeSeq[extract_name]=str(name_infor.seq)
    if extract_name in seqStrand:
        name_infor.seq=name_infor.seq.complement()
        beforeModSeq[extract_name]=name_infor
    else:
        beforeModSeq[extract_name]=name_infor

SeqIO.write(beforeModSeq.values(), 'finalPremiRNAsMod.fasta', "fasta")
