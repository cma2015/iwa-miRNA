import re
import os
import subprocess
import multiprocessing

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-o','--outpath',action = 'store',type = "string" ,dest = 'outputpath')
parser.add_option('-v','--version', action="store_false", dest="verbose", default='',help="version [default]")
(options,args)=parser.parse_args()

outpath = options.outputpath

def centroid_fold(tmp_dict):
    ## Centroidfold structure
    global centroid_str
    centroid_str = {}
    with open("{}/tmp.fasta".format(outpath), "w") as fas_file:
        for rna_site in tmp_dict:
            fas_file.write('>%s\n%s\n' % (rna_site, tmp_dict[rna_site][2]))
    fas_str = 'centroid_fold -e CONTRAfold -g 4 {}/tmp.fasta > {}/tmp.structure'.format(outpath, outpath)
    other, res = subprocess.getstatusoutput(fas_str)
    with open("{}/tmp.structure".format(outpath), "r") as str_file:
        for i in str_file.readlines():
            if i[:1] == '>':
                name=i.replace('>','').split()[0]
                centroid_str[name]=''
            elif i[:1] in [".", "("]:
                centroid_str[name]+=i.split(" ")[0]
    other, res = subprocess.getstatusoutput('rm {}/tmp.fasta {}/tmp.structure'.format(outpath, outpath))


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


# structure filtering
def struct_filter(test_seq, mature_mi_rna):
    mi_rna_pos = getMiRNAPosition(mature_mi_rna, test_seq[1])
    mi_rna_bracket = len(set(test_seq[2][mi_rna_pos[0]:mi_rna_pos[1]]))
    mi_star = getMiRNAStar(test_seq[1], test_seq[2], mature_mi_rna)
    mi_star_pos = getMiRNAPosition(mi_star, test_seq[1])
    mi_star_bracket = len(set(test_seq[2][mi_star_pos[0]:mi_star_pos[1]]))
    if mi_rna_bracket == 3 or mi_star_bracket == 3:
        return(['False', '-', '-', 'nonor'])
    else:
        pairStr = matchRNAStructure(test_seq[2])
        duplexS = mi_rna_pos[0]
        duplexE = mi_rna_pos[1]-2
        location=pairStr[duplexS:duplexE]
        location_bc=pairStr[duplexS:duplexE]
        while -1 in location_bc:
            location_bc.remove(-1)
        if all([location_bc[i] > location_bc[i+1] for i in range(len(location_bc)-1)]):
            if mi_rna_pos[0] > mi_star_pos[1]:
                loop_len = mi_rna_pos[0]-mi_star_pos[1]
                loop_ture = pairStr[(mi_star_pos[1]):(mi_star_pos[1]+3)]
            else:
                loop_len = mi_star_pos[0]-mi_rna_pos[1]
                loop_ture = pairStr[(mi_rna_pos[1]):(mi_rna_pos[1]+3)]
            if loop_len <= MIN_SPACE and loop_ture.count(-1) >0:
                return(['False', str(loop_len), str(loop_ture.count(-1)), str('short_loop')])
            elif(location.count(-1) >5):
                return(['False', '-', '-', '-'])
            else:
                mismatch=0
                while location[0] == -1:
                    mismatch+=1
                    del location[0]
                bugleOne=0
                bugleTwo=0
                nomatch=0
                for num in range(len(location)-1):
                    if(location.count(-1)>5):
                        mismatch = 99
                        break
                    if(location[num]-location[num+1]==1):
                        nomatch=0
                    elif(location[num]-location[num+1]==2):
                        bugleOne+=1
                    elif(location[num]-location[num+1]>2 and location[num] != -1 and location[num+1] != -1):
                        bugleLen = location[num]-location[num+1]-1
                        bugleTwo+=bugleLen
                        # break
                    elif(location[num+1] == -1):
                        nomatch+=1
                        if(location[num] != -1):
                            leftnum = location[num]
                        elif(num+1 == len(location)):
                            comstrand = leftnum-location[num+1]-1
                            mismatch = 99
                            break
                    else:
                        comstrand = leftnum-location[num+1]-1
                        if(abs(nomatch-comstrand)==1):
                            bugleOne+=1
                        elif(abs(nomatch-comstrand)>1):
                            bugleTwo+=abs(nomatch-comstrand)
                        else:
                            mismatch+=min(nomatch, comstrand)
                            nomatch=0
                if mismatch+bugleOne+bugleTwo <= 5 and bugleOne+bugleTwo <= 3:
                    return(['True', str(mismatch+bugleOne), str(bugleOne), str(bugleTwo)])
                elif(mismatch==99):
                    return(['False', '-', '-', '-'])
                else:
                    return(['False', str(mismatch+bugleOne), str(bugleOne), str(bugleTwo)])
        else:
            return(['False', '-', '-', 'nonor'])


def filter_feature(results_dict):
    rnafold_pass = {}
    for pre_site in results_dict:
        rnafold_res = struct_filter(results_dict[pre_site][0:3], results_dict[pre_site][4])
        if rnafold_res[0] != 'False' :
            rnafold_pass[pre_site] = results_dict[pre_site][0:11]
            rnafold_pass[pre_site].extend(rnafold_res)
    centroidfold_pass = {}
    for pre_site in rnafold_pass:
        all_info = rnafold_pass[pre_site]
        pre_centroid = centroid_str[pre_site]
        centroid_res = struct_filter([all_info[0], all_info[1], pre_centroid], all_info[4])
        if centroid_res[0] != 'False' :
            centroidfold_pass[pre_site] = all_info[0:15]
            centroidfold_pass[pre_site].extend(centroid_res)
    mature_len_pass = {}
    for pre_site in centroidfold_pass:
        info_pre = centroidfold_pass[pre_site]
        if len(info_pre[4]) > 22 and info_pre[8] == 0:
            pass
        else:
            # if int(info_pre[12]) >0 and int(info_pre[16]) >0:
            mature_len_pass[pre_site] = info_pre
    overhang_pass = {}
    for pre_loc in mature_len_pass:
        pre_info = mature_len_pass[pre_loc]
        arm_info = checkArm(pre_info[1], pre_info[2], pre_info[4])
        if arm_info == '5p':
            [mir_start, mir_end] = getMiRNAPosition( pre_info[4], pre_info[1])
            match_struct = matchRNAStructure(pre_info[2])
            [match_start, match_end] = getMatchedPositions( mir_start, mir_end-2, match_struct)
            match_start = match_start + 2
            match_end = match_end
            if match_start+1 <= len(pre_info[1]):
                overhang_pass[pre_loc] = pre_info
        elif arm_info == '3p':
            overhang_pass[pre_loc] = pre_info
    # print('rnafold_pass', len(rnafold_pass))
    # print('centroidfold_pass', len(centroidfold_pass))
    # print('mature_len_pass', len(mature_len_pass))
    # print('overhang_pass', len(overhang_pass))
    return(overhang_pass)


def mat_star(pre_info):
    if pre_info[10] == '5p':
        mat_seq = pre_info[5]
        star_seq = pre_info[8]
    elif pre_info[10] == '3p':
        mat_seq = pre_info[8]
        star_seq = pre_info[5]
    else:
        mat_seq = pre_info[5]
        star_seq = pre_info[8]
    return([mat_seq, star_seq])


def create_report(dict_key):
    info_list = input_dict[dict_key]
    pre_name = re.sub(':','_',info_list[1])
    rnafoldstr = '{}/png/{}_r'.format(outpath, pre_name)
    centroidstr = '{}/png/{}_c'.format(outpath, pre_name)
    a,process_res = subprocess.getstatusoutput('echo {}| RNAfold --noPS '.format(info_list[2]))
    rnafold_res = re.split("\n|\s\(", process_res)
    # print(rnafold_res[1])
    with open(rnafoldstr+'.str', 'w') as str_f:
        str_f.write('>{}_r\n{}\n{}'.format(pre_name,
                                        info_list[2],
                                        rnafold_res[1]))
    with open(centroidstr+'.str', 'w') as str_f:
        str_f.write('>{}_c\n{}\n{}'.format(pre_name,
                                        info_list[2],
                                        centroid_str[dict_key]))
    mature_se = getMiRNAPosition(mat_star(info_list)[0],info_list[2])
    star_se = getMiRNAPosition(mat_star(info_list)[1],info_list[2])
    site_loc = sorted([0, mature_se[0], mature_se[1], star_se[0],
                star_se[1], len(info_list[2])])
    if info_list[10] == "5p":
        sign1_name = "Mature"
        sign2_name = "Star"
        sign1_col = ['0.8', '0', '0']
        sign2_col = ['0', '0.8', '0']
    elif info_list[10] == "3p":
        sign1_name = "Star"
        sign2_name = "Mature"
        sign1_col = ['0', '0.8', '0']
        sign2_col = ['0.8', '0', '0']
    else:
        sign1_col = ['0', '0.8', '0']
        sign2_col = ['0', '0.8', '0']
    color_list = [['0.8', '0.8', '0.8'], sign1_col, ['1', '0.6', '0'], sign2_col, ['0.8', '0.8', '0.8']]
    color_tab = []
    for i in range(len(site_loc)-1):
        for j in range(len(info_list[2])):
            if site_loc[i] <= j < site_loc[i+1]:
                color_tab.append(str(j+1))
                color_tab.extend(color_list[i])
                color_tab.append('cfmark')
    extraMacro = "/cfmark {setrgbcolor newpath 1 sub coor exch get aload pop fsize 2 div 0 360 arc fill} bind def"
    os.system("RNAplot --pre \" {} {}\" < {}.str".format(extraMacro, " ".join(color_tab), rnafoldstr))
    os.system("RNAplot --pre \" {} {}\" < {}.str".format(extraMacro, " ".join(color_tab), centroidstr))
    os.system("rm {}.str {}.str".format(rnafoldstr, centroidstr))
    os.system("mv {0}_r_ss.ps {1}.ps &&mv {0}_c_ss.ps {2}.ps".format(
        pre_name, rnafoldstr, centroidstr))
    os.system('convert -trim -quality 300  -density 300 ' + rnafoldstr + '.ps ' + rnafoldstr + '.png')
    os.system('convert -trim -quality 300  -density 300 ' + centroidstr + '.ps ' + centroidstr + '.png')
    os.system("rm {0}.ps".format(rnafoldstr))
    os.system("rm {0}.ps".format(centroidstr))

# os.system("{1}gmt psconvert -Tf {0}.ps && {1}pdf2svg {0}.pdf {0}.svg".format(rnafoldstr, '~/sRNAbox/miniconda3/bin/'))
# os.system("{1}gmt psconvert -Tf {0}.ps && {1}pdf2svg {0}.pdf {0}.svg".format(centroidstr, '~/sRNAbox/miniconda3/bin/'))
# os.system("rm {0}.ps {0}.pdf ".format(rnafoldstr))
# os.system("rm {0}.ps {0}.pdf ".format(centroidstr))

if __name__ == '__main__':
    global input_dict
    input_dict = {}
    with open('{}/Translate_result.txt'.format(outpath), 'r') as mergefile:
        for eli in mergefile.readlines()[1:]:
            eli = eli.strip().split('\t')
            input_dict[eli[1]] = eli
    centroid_fold(input_dict)
    pool = multiprocessing.Pool(processes=2)
    pool.map(create_report, list(input_dict.keys()))
