import re
import os
import sys
import time
import subprocess

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

def parse_mirdp(input_file, fai_file):
    fai_value = {}
    with open(fai_file) as fai_input:
        for eli in fai_input.readlines():
            eli = eli.strip().split('\t')
            fai_value[eli[0]] = int(eli[1])
    '''Read miRDP2 result'''
    out_dic = {}
    with open(input_file + '.bed', 'w') as bed_out:
        with open(input_file) as infile:
            for eachline in infile.readlines():
                el = eachline.strip().split('\t')
                pre_loc = el[5].split('.')
                mature_loc = el[4].split('.')
                if el[1] == '+':
                    pre_start = str(int(pre_loc[0])-1)
                    pre_end = str(int(pre_loc[2])-1)
                    mature_start = str(int(mature_loc[0])-1)
                    mature_end = str(int(mature_loc[2])-1)
                else:
                    pre_start = str(int(pre_loc[0])+2)
                    pre_end = str(int(pre_loc[2])+2)
                    mature_start = str(int(mature_loc[0])+2)
                    mature_end = str(int(mature_loc[2])+2)
                if int(pre_start)-flank-1 < 0:
                    extend_left = str(0)
                else:
                    extend_left = str(int(pre_start)-flank-1)
                if int(pre_end)+flank > fai_value[el[0]]:
                    extend_right = str(fai_value[el[0]])
                else:
                    extend_right = str(int(pre_end)+flank)
                pre_loc = ':'.join([el[0], pre_start, pre_end, el[1]])
                bed_out.write('\t'.join([el[0], extend_left, extend_right, 
                pre_loc, '.', el[1]]) + '\n')
                a,process_res = subprocess.getstatusoutput('echo {}| RNAfold --noPS '.format(el[7]))
                rnafold_res = re.split("\n|\s\(", process_res)
                mature_loc = ':'.join([el[0], mature_start, mature_end, el[1]])
                mature_abu = re.split("x", el[2])
                # print(rnafold_res)
                # print(pre_loc, el[7].upper(), rnafold_res[1], mature_loc, el[6].upper(), mature_abu[1])
                info_list = [pre_loc, el[7].upper(), rnafold_res[1], mature_loc, el[6].upper(), mature_abu[1]]
                info_list = star_info(info_list)
                mfe = re.sub(r'\)|\s', '', rnafold_res[2])
                info_list.append(mfe)
                info_list.append(float(mfe)/len(el[7])*100)
                info_list.append(checkArm(el[7].upper(), rnafold_res[1], el[6].upper()))
                info_list.append(checkArm(el[7].upper(), rnafold_res[1], info_list[7]))
                out_dic.update({info_list[0]:info_list})
        return(out_dic)

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
            eli = eli.strip().split(' ')
            eli[0] = eli[0].strip('>')
            seq_id = str(eli[0])
            # print(eli)
            seq_count = round(float(eli[1]),3)
        else:
            tagId_seq[seq_id] = eli.strip()
            tagId_count[eli.strip()] = seq_count

parse_dict = parse_mirdp( options.txtfile, options.faifile )

print('bedtools getfasta -name -s -tab -fi {0} -bed {1}.bed | sed -r "s/::.+?\\t/\\t/" > {1}.fa'.format(options.genome, options.txtfile))

os.system('bedtools getfasta -name -s -tab -fi {0} -bed {1}.bed | sed -r "s/::.+?\\t/\\t/" > {1}.fa'.format(options.genome, options.txtfile))

new_seq = {}
with open( '{}.fa'.format(options.txtfile) ) as fa_input:
    for eli in fa_input.readlines():
        eli = eli.strip().split('\t')
        loc_tmp = re.sub("\(.+", "", eli[0])
        # print(loc_tmp)
        new_seq[loc_tmp] = eli[1].upper()

with open( options.output , 'w') as outinfo:
    outinfo.write('## mirGFF3\n## Genome: \n## Data:{}\n'.format(
    time.strftime("%Y-%m-%d", time.localtime())))
    for i in new_seq:
        parse_dict[i][1] = new_seq[i]
        tmp_list = i.split(':')
        pre_miRNA_loc = ':'.join([tmp_list[0], str(int(tmp_list[1])-flank), str(int(tmp_list[2])+flank), tmp_list[3]])
        con1 = parse_dict[i][11] == '3p' and parse_dict[i][12] == '5p'
        con2 = parse_dict[i][11] == '5p' and parse_dict[i][12] == '3p'
        if con1 or con2:
            ## pre-miRNA features
            pre_miRNA = pre_miRNA_loc.split(':')
            outinfo.write('{}\t{}\tpre-miRNA\t{}\t{}\t.\t{}\t.\tID={};Seq={};Structure={};MFE={};AMFE={}\n'.format(
            pre_miRNA[0],"miRDeeP-P2", pre_miRNA[1], pre_miRNA[2],pre_miRNA[3], pre_miRNA_loc,
            parse_dict[i][1], parse_dict[i][2], parse_dict[i][9], parse_dict[i][10]))
            ## mature miRNA feature
            mature_mir = parse_dict[i][3].split(':')
            outinfo.write('{}\t{}\tref_miRNA\t{}\t{}\t.\t{}\t.\tID={};Parent={};Seq={};Arm={}\n'.format(
            mature_mir[0],"miRDeeP-P2", mature_mir[1], mature_mir[2],mature_mir[3],parse_dict[i][3],
            pre_miRNA_loc, parse_dict[i][4], parse_dict[i][11]))
            ## star features
            star_mir = parse_dict[i][6].split(':')
            outinfo.write('{}\t{}\tstar_miRNA\t{}\t{}\t.\t{}\t.\tID={};Parent={};Seq={};Arm={}\n'.format(
            star_mir[0],"miRDeeP-P2", star_mir[1], star_mir[2],star_mir[3],parse_dict[i][6],
            pre_miRNA_loc, parse_dict[i][7], parse_dict[i][12]))
        else:
            print(parse_dict[i])
