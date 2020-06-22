#!/usr/bin/python
#-*-coding:utf-8 -*-

from reportlab.pdfgen import canvas
import os
import re
import sys
import copy
import commands
from PIL import Image
from Bio import SeqIO

def format_values_ps(value):
    val_tab = [" %s 0 0.8 0 cfmark" % v for i, v in enumerate(value)]
    return "".join(val_tab)

def text_location(text, font, size, color, width, height):
    t = c.beginText()
    t.setFont(font, size)
    t.setCharSpace(0)
    t.setFillColor(color)
    t.setTextOrigin(width, height)
    t.textLine(text)
    return t

order = 1
page_wid = 600
page_hei = 300
## file path
beforeSeq={}
afterSeq={}
seqStrand={}

bed_file = open('premirLoc.bed', "r")
for i in bed_file.readlines():
    i = i.strip()
    lineStrand = i.split("\t")
    if lineStrand[5]=="-":
        seqStrand[lineStrand[3]]="-"
    
bed_file.close()

for name_infor in SeqIO.parse('finalPremiRNAsMod.fasta', 'fasta'):
    extract_name = name_infor.id
    if extract_name in seqStrand:
        beforeSeq[extract_name]=str(name_infor.seq.complement())
    else:
        beforeSeq[extract_name]=str(name_infor.seq)

for name_mod in SeqIO.parse('finalPremiRNAsMod_SNPs.fasta', 'fasta'):
    extract_mod = name_mod.id
    if extract_mod in seqStrand:
        afterSeq[extract_mod]=str(name_mod.seq.complement())
    else:
        afterSeq[extract_mod]=str(name_mod.seq)

vcfLoc={}
vcfFile = open('SNP.vcf', "r")
for i in vcfFile.readlines():
    i = i.strip()
    if not i.startswith("#"):
        lineSNP = i.split("\t")
        if lineSNP[0] not in vcfLoc:
            vcfLoc[lineSNP[0]]=[int(lineSNP[1])]
        else:
            vcfLoc[lineSNP[0]].append(int(lineSNP[1]))

location = {}
location_file = open('mature_location.txt', "r")
for i in location_file.readlines():
    i = i.strip()
    line_loc = i.split("\t")
    location[line_loc[0]] = [int(line_loc[1]), int(line_loc[2]), int(line_loc[3]), int(line_loc[4])]
location_file.close()

## each per-miRNAs
for extract_name in afterSeq.keys():
    for seqType in ['before', 'after']:
        file_name = extract_name
        if seqType == 'before':
            sequence = beforeSeq[extract_name]
        else:
            sequence = afterSeq[extract_name]
        # RNAfold struture
        other, structure = commands.getstatusoutput("echo " + sequence + "| RNAfold")
        structure = re.split("\s", structure)
        ## Write the structure
        inputFile = file_name+'.str'
        fas_file = open(inputFile, "w")
        fas_file.write(">"+extract_name+"\n")
        fas_file.write(sequence+"\n")
        fas_file.write(structure[1] + " " + structure[2])
        fas_file.close()
        ## Calculate the color
        mature_location = sorted(location[extract_name])
        values = range(mature_location[0]+1, mature_location[1]+1)
        values.extend(range(mature_location[2]+1, mature_location[3]+1))
        color_tab = []
        sign1_rgd = ['0.8', '0', '0']
        sign2_rgd = ['0', '0.8', '0']
        pre_location = copy.deepcopy(location[extract_name])
        pre_location = sorted(pre_location)
        print vcfLoc[extract_name]
        for i in [i+1 for i in range(len(sequence))]:
            if i not in vcfLoc[extract_name]:
                if pre_location[0] < i < (pre_location[1]+1):
                    color_tab.append(str(i))
                    color_tab.extend(sign1_rgd)
                    color_tab.append('cfmark')
                elif pre_location[2] < i < (pre_location[3]+1):
                    color_tab.append(str(i))
                    color_tab.extend(sign2_rgd)
                    color_tab.append('cfmark')
                elif pre_location[1] < i < (pre_location[2]+1):
                    color_tab.extend([str(i), '1', '0.6', '0', 'cfmark'])
                else:
                    color_tab.extend([str(i), '0.8', '0.8', '0.8', 'cfmark'])
            else:
                color_tab.extend([str(i), '0.4', '0.6', '0.9', 'cfmark'])
        #print " ".join(color_tab)
        extraMacro = "/cfmark {setrgbcolor newpath 1 sub coor exch get aload pop fsize 2 div 0 360 arc fill} bind def"
        os.system("RNAplot --pre \" %s %s\" < %s" % (extraMacro, " ".join(color_tab), inputFile))
        if seqType == 'before':
            os.system('convert -trim -quality 500  -density 500 ' + extract_name + '_ss.ps ' + file_name + '_before.png')
            img = Image.open( file_name + '_before.png')
            (x, y) = img.size
            scaleBefore = float(max(x, y))/min(x, y)        
        else:
            os.system('convert -trim -quality 500  -density 500 ' + extract_name + '_ss.ps ' + file_name + '_after.png')
            img = Image.open( file_name + '_after.png')
            (x, y) = img.size
            scaleAfter = float(max(x, y))/min(x, y)  
    """The page
    """
    c = canvas.Canvas(extract_name + ".pdf", pagesize=(page_wid, page_hei))
    # pre-miRNAs and structures
    title = text_location(extract_name, 'Times-Bold', 15, "black", 290, 270)
    c.drawText(title)
    beforeText = text_location('Ref', 'Times-Bold', 15, "black", 200-1.5*(200/scaleBefore), 10)
    c.drawText(beforeText)
    AfterText = text_location('Alt', 'Times-Bold', 15, "black", 400, 10)
    c.drawText(AfterText)
    c.drawImage(file_name + '_before.png', 200-1.5*(200/scaleBefore), 30, width=(200/scaleBefore), height=200)
    c.drawImage(file_name + '_after.png', 400, 30, width=(200/scaleAfter), height=200)
    c.showPage()
    c.save()

# if extract_name in show_value:
#     loc_info = tab_info[1].split(":")
#     title = text_location(tab_info[21], 'Courier', 15, "black", word1_site, word1_height)
#     chr = text_location("Chromosome:    chr" + loc_info[0], 'Courier', 10, "black", word1_site, word1_height-40)
#     sta = text_location("Start:         " + loc_info[1], 'Courier', 10, "black", word1_site, word1_height-50)
#     end = text_location("End:           " + loc_info[2], 'Courier', 10, "black", word1_site, word1_height-60)
#     strand = text_location("Strand:        " + loc_info[3], 'Courier', 10, "black", word1_site, word1_height-70)
#     sbias = text_location("Strand bias:   " + str(round(float(tab_info[5]), 3)),
#                           'Courier', 10, "black", word1_site, word1_height-80)
#     abias = text_location("Abundance bias:" + str(round(float(tab_info[6]), 3)),
#                           'Courier', 10, "black", word1_site, word1_height-90)
#     mfe = text_location("MFE:           " + mfe_val, 'Courier', 10, "black", word1_site, word1_height-100)
#     c.drawText(title)
#     c.drawText(chr)
#     c.drawText(sta)
#     c.drawText(end)
#     c.drawText(strand)
#     c.drawText(sbias)
#     c.drawText(abias)
#     c.drawText(mfe)
