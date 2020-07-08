import re
import os
import sys
import time
import subprocess

# Common RNA function
def isParenthese(ch):
    if( ch=='(' or ch==')' ):
        return True
    else:
        return False

def matchRNAStructure( RNAStructure ):
    RNALen = len(RNAStructure)
    matchedPosList = [-1]*RNALen
    stack = []
    stackPos = []
    for i in range(RNALen):
        a = RNAStructure[i]
        if( isParenthese(a) == False ):
            continue
        if( len(stack) == 0 ):
            #empty stack, record item
            stack.append(a)
            stackPos.append(i)
        else:
            #pop last item in stack and stackPos
            stack_lastItem = stack.pop()
            stackPos_lastItem = stackPos.pop()
            if( stack_lastItem == '(' and a == ')' ):
                #meet a match, record matched position
                matchedPosList[i] = stackPos_lastItem
                matchedPosList[stackPos_lastItem] = i
                continue
            else:
                #have to record
                stack.append( stack_lastItem )
                stackPos.append( stackPos_lastItem )
                stack.append(a)
                stackPos.append(i)
    return matchedPosList

def getMatchedPositions( startPos, endPos, matchedStructList ):
    RNALen = len(matchedStructList)
    if( startPos < 0 ):
        startPos = 0
    if( endPos < 0 ):
        endPos = endPos
    if( startPos > RNALen ):
        startPos = RNALen - 1
    if( endPos > RNALen ):
        endPos = RNALen - 1
    if( startPos == 0 and endPos == 0 ):
        return [0,0]
    matchedPosList = []
    for i in range(startPos, endPos):
        curPos = matchedStructList[i]
        matchedPosList.append( curPos )
    if( matchedPosList[0] == -1 ):
        idx = 0
        for i in range(len(matchedPosList)):
            if( matchedPosList[i] != -1 ):
                idx = i
                break
        matchedPos = matchedPosList[idx]
        idxList = range(idx)
        for i in idxList[::-1]:
            matchedPos = matchedPos + 1
            matchedPosList[i] = matchedPos
    if( matchedPosList[-1] == -1 ):
        idxList = range(len(matchedPosList))
        idx = 0
        for i in idxList[::-1]:
            if( matchedPosList[i] != -1 ):
                idx = i
                break
        #end for i
        matchedPos =  matchedPosList[idx]
        idx = idx + 1
        for i in range(idx, len(matchedPosList)):
            matchedPos = matchedPos - 1
            matchedPosList[i] = matchedPos
    #end if
    if( matchedPosList[0] < 0 ):
        matchedPosList[0] = 0
    if( matchedPosList[0] > RNALen ):
        matchedPosList[0] = RNALen - 1
    if( matchedPosList[-1] < 0 ):
        matchedPosList[-1] = 0
    if( matchedPosList[-1] > RNALen ):
        matchedPosList[-1] = RNALen - 1
    return [matchedPosList[0], matchedPosList[-1]]

def getMiRNALen (miRNASeq):
    return len(miRNASeq)

def getMiRNAPosition (miRNASeq, RNASequence):
    startPos = RNASequence.find(miRNASeq)
    endPos = startPos + getMiRNALen(miRNASeq)
    return [startPos, endPos]

def getMiRNAStar(RNASequence, RNAStructure, miRNASeq, matchedStructList = [] ):
    if( len(matchedStructList) == 0 ):
        matchedStructList = matchRNAStructure( RNAStructure )
    [miRNAPosStart, miRNAPosEnd] = getMiRNAPosition( miRNASeq, RNASequence )
    [matchedPosStart, matchedPosEnd] = getMatchedPositions( miRNAPosStart,
        miRNAPosEnd-2, matchedStructList )
    matchedPosStart = matchedPosStart + 2
    matchedPosEnd = matchedPosEnd
    RNALen = len(RNASequence)
    if( matchedPosStart < 0 ):
        matchedPosStart = 0
    if( matchedPosEnd < 0 ):
        matchedPosEnd = 0
    if( matchedPosStart > RNALen ):
        matchedPosStart = RNALen - 1
    if( matchedPosEnd > RNALen ):
        matchedPosEnd = RNALen - 1
    starSeq = ""
    if( matchedPosStart < matchedPosEnd ):
        starSeq = RNASequence[matchedPosStart: (matchedPosEnd+1)]
        starSeq = starSeq[::-1]
    else:
        starSeq = RNASequence[matchedPosEnd: (matchedPosStart+1)]
    return starSeq

def getMiRNAStructure (miRNASeq, RNASequence, RNAStructure ):
    [startPos, endPos] = getMiRNAPosition( miRNASeq, RNASequence )
    mirnaStruct = RNAStructure[ startPos:endPos ]
    return mirnaStruct

def getMirnaIncludedInLoop( RNASequence, RNAStructure, miRNASeq ):
    flag = False
    mirnaStruct = getMiRNAStructure( miRNASeq, RNASequence, RNAStructure )
    if( ('(' in mirnaStruct) and (')' in mirnaStruct) ):
        flag = True
    return flag

def checkArm( RNASequence, RNAStructure, miRNASeq ):
    armDetailedType = ""
    if getMirnaIncludedInLoop(RNASequence, RNAStructure, miRNASeq):
        armDetailedType = "loop"
    #check arms
    mirnaStruct = getMiRNAStructure(miRNASeq, RNASequence, RNAStructure)
    armDetailedType = "unmatchedRegion"
    if list(mirnaStruct).count("(") > len(mirnaStruct)/2:  # '(' in mirnaStruct:
        armDetailedType = "5p"
    if list(mirnaStruct).count(")") > len(mirnaStruct)/2:  # ')' in mirnaStruct:
        armDetailedType = "3p"
    return armDetailedType

def star_info(input_list):
    '''miRNA* location, sequence, and abundance'''
    star_seq = getMiRNAStar(input_list[1], input_list[2], input_list[4])
    star_pos = getMiRNAPosition(star_seq, input_list[1])
    pre_info = re.split(':', input_list[0])
    if pre_info[3] == '+':
        star_location = ':'.join([pre_info[0],
        str(int(pre_info[1])+star_pos[0]),
        str(int(pre_info[1])+star_pos[1]-1), pre_info[3]])
    else:
        star_location = ':'.join([pre_info[0],
        str(int(pre_info[2])-star_pos[1]+1),
        str(int(pre_info[2])-star_pos[0]), pre_info[3]])
    if star_seq in tagId_count:
        input_list.extend([star_location, star_seq, str(tagId_count[star_seq])])
    else:
        input_list.extend([star_location, star_seq, '0'])
    return(input_list)

def parse_mircat(input_file, fai_file):
    fai_value = {}
    with open(fai_file) as fai_input:
        for eli in fai_input.readlines():
            eli = eli.strip().split('\t')
            fai_value[eli[0]] = int(eli[1])
    '''Read mircat2 result'''
    out_dic = {}
    with open(input_file + '.bed', 'w') as bed_out:
        with open(input_file) as infile:
            for eachline in infile.readlines()[1:]:
                el = eachline.strip().split(',')
                el[1] = el[1].split(" ")[0]
                pre_start = el[9]
                pre_end = el[10]
                mature_start = el[4]
                mature_end = el[5]
                pre_loc = ':'.join([el[1], pre_start, pre_end, el[6]])
                a,process_res = subprocess.getstatusoutput('echo {}| RNAfold --noPS '.format(el[8]))
                rnafold_res = re.split("\n|\s\(", process_res)
                mature_loc = ':'.join([el[1], mature_start, mature_end, el[6]])
                mature_abu = re.split("x", el[2])
                # print(rnafold_res)
                info_list = [pre_loc, el[8].upper(), rnafold_res[1], mature_loc, el[2].upper(), tagId_count[el[2].upper()]]
                if el[15] == "N/A":
                    info_list = star_info(info_list)
                else:
                    star_loc = ':'.join([el[1], el[17], el[18], el[6]])
                    info_list.extend([star_loc, el[15].upper(), tagId_count[el[15].upper()]])
                info_list.append(checkArm(el[8].upper(), rnafold_res[1], el[2].upper()))
                mat_tmp = info_list[3].split(":")
                star_tmp = info_list[6].split(":")
                mat_star = [int(mat_tmp[1]), int(mat_tmp[2]), int(star_tmp[1]), int(star_tmp[2])]
                mat_site = getMiRNAPosition(info_list[4], info_list[1])
                star_site = getMiRNAPosition(info_list[7], info_list[1])
                mat_star_site = [int(mat_site[0]), int(mat_site[1]), int(star_site[0]), int(star_site[1])]
                info_list[0] = ":".join([el[1], str(min(mat_star)), str(max(mat_star)), el[6]])
                info_list[1] = el[8].upper()[min(mat_star_site):max(mat_star_site)]
                info_list[2] = "-"
                ## expand
                if min(mat_star)-flank-1 < 0:
                    extend_left = str(0)
                else:
                    extend_left = str(min(mat_star)-flank)
                if max(mat_star)+flank > fai_value[el[1]]:
                    extend_right = str(fai_value[el[1]])
                else:
                    extend_right = str(max(mat_star)+flank)
                bed_out.write('\t'.join([el[1], extend_left, extend_right,
                             info_list[0], '.', el[6]]) + '\n')
                out_dic.update({info_list[0]:info_list})
        return(out_dic)

from optparse import OptionParser
import os

parser = OptionParser()
parser.add_option('-g','--genome',action = 'store',type = "string" ,dest = 'genome')
parser.add_option('-l','--faifile',action = 'store',type = "string" ,dest = 'faifile')
parser.add_option('-f','--fasta',action = 'store',type = "string" ,dest = 'fasta')
parser.add_option('-t','--txtfile',action = 'store',type = "string" ,dest = 'txtfile')
parser.add_option('-o','--output',action = 'store',type = "string" ,dest = 'output')
parser.add_option('-v','--version', action="store_false", dest="verbose", default='',help="version [default]")
(options,args)=parser.parse_args()

global flank
global tagId_seq
global tagId_count
flank = 20
tagId_seq, tagId_count = {}, {}
with open( options.fasta ) as ffread:
    for eli in ffread.readlines():
        if eli.startswith('>'):
            eli = eli.strip().split('(')
            eli[0] = eli[0].strip('>')
            eli[1] = eli[1].split(')')[0]
            seq_id = str(eli[0])
            # print(eli)
            seq_count = round(float(eli[1]),3)
        else:
            tagId_seq[seq_id] = eli.strip()
            tagId_count[eli.strip()] = seq_count

parse_dict = parse_mircat( options.txtfile, options.faifile )

print('bedtools getfasta -name -s -tab -fi {0} -bed {1}.bed | sed -r "s/::.+?\\t/\\t/" > {1}.fa'.format(options.genome, options.txtfile))

os.system('bedtools getfasta -name -s -tab -fi {0} -bed {1}.bed | sed -r "s/::.+?\\t/\\t/" > {1}.fa'.format(options.genome, options.txtfile))

new_seq = {}
with open( '{}.fa'.format(options.txtfile) ) as fa_input:
    for eli in fa_input.readlines():
        eli = eli.strip().split('\t')
        loc_tmp = re.sub("\(.+", "", eli[0])
        # print(loc_tmp)
        new_seq[loc_tmp] = eli[1].upper()

with open( options.output , 'w') as out_info:
    out_info.write('Precursors\tpLocation\tpSequence\tpLength\t\
Location5p\tSequence5p\tLength5p\tAbundance5p\tLocation3p\tSequence3p\tLength3p\tAbundance3p\tMature_arm\n')
    for i in new_seq:
        tmp_list = i.split(':')
        pre_miRNA_loc = ':'.join([tmp_list[0], str(int(tmp_list[1])-flank),
                                  str(int(tmp_list[2])+flank), tmp_list[3]])
        pre_miRNA_seq = new_seq[i]
        pre_miRNA_len = len(new_seq[i])
        if parse_dict[i][9] == "5p" or parse_dict[i][9] == "3p":
            out_info.write('{}\n'.format('\t'.join([pre_miRNA_loc, pre_miRNA_loc,
                                                    pre_miRNA_seq, str(pre_miRNA_len),
                                                   parse_dict[i][3], parse_dict[i][4],
                                                    str(len(parse_dict[i][4])), str(parse_dict[i][5]),
                                                   parse_dict[i][6], parse_dict[i][7],
                                                    str(len(parse_dict[i][7])), str(parse_dict[i][8]),
                                                   parse_dict[i][9]])))
