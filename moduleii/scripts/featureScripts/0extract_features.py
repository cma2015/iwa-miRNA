import re
import os
import subprocess
import math
from tqdm import tqdm

def bracket_or_not(single_pos):
    if single_pos == '(' or single_pos == ')':
        return(True)
    else:
        return(False)
    
def parse_rna_str(rna_str):
    rna_length = len(rna_str)
    matched_pos_list = [-1]*rna_length
    stack = []
    stackPos = []
    for i in range(rna_length):
        a = rna_str[i]
        if bracket_or_not(a) == False:
            continue
        if len(stack) == 0:
            stack.append(a)
            stackPos.append(i)
        else:
            stack_lastItem = stack.pop()
            stackPos_lastItem = stackPos.pop()
            if stack_lastItem == '(' and a == ')':
                matched_pos_list[i] = stackPos_lastItem
                matched_pos_list[stackPos_lastItem] = i
                continue
            else:
                stack.append(stack_lastItem)
                stackPos.append(stackPos_lastItem)
                stack.append(a)
                stackPos.append(i)
    return(matched_pos_list)

def gc_content (rna_seq):
    gc_val = rna_seq.count('G')
    gc_val += rna_seq.count('g')
    gc_val += rna_seq.count('C')
    gc_val += rna_seq.count('c')
    return(round(1.0*gc_val/len(rna_seq), 4))

def get_mononuc_ratio(rna_seq, mnn):
    mnv = rna_seq.count(mnn)
    return round(1.0*mnv/len(rna_seq), 4)

def mononucleotide_ratio(rna_seq):
    ratio_A = get_mononuc_ratio(rna_seq, "A")
    ratio_C = get_mononuc_ratio(rna_seq, "C")
    ratio_G = get_mononuc_ratio(rna_seq, "G")
    ratio_U = get_mononuc_ratio(rna_seq, "U")
    return({'sub_A': ratio_A, 'sub_C': ratio_C, 
            'sub_G': ratio_G, 'sub_U': ratio_U})

##get dinucleotide content
def get_dinucleotide_ratio(rna_seq, mnn1, mnn2):
    sub_str = mnn1 + mnn2
    mnv = rna_seq.count(sub_str)
    return(round(1.0*mnv/len(rna_seq), 4))

def dinucleotide_ratio(rna_seq):
    nus = ["A", "C", "G", "U"]
    tmp_dict = {}
    for i in nus:
        for j in nus:
            tmp_dict['sub_' + i + j] = get_dinucleotide_ratio(rna_seq, i, j)
    return(tmp_dict)

def pos_in_pre(mature_seq, rna_seq):
    start_pos = rna_seq.find(mature_seq)
    end_pos = start_pos + len(mature_seq)
    return([start_pos, end_pos])


def get_matched_pos(start_pos, end_pos, matched_str):
    rna_length = len(matched_str)
    if start_pos < 0:
        start_pos = 0
    if end_pos < 0:
        end_pos = end_pos
    if start_pos > rna_length:
        start_pos = rna_length - 1
    if end_pos > rna_length:
        end_pos = rna_length - 1
    if start_pos == 0 and end_pos == 0:
        return([0,0])
    matched_pos_list = []
    for i in range(start_pos, end_pos):
        curPos = matched_str[i]
        matched_pos_list.append(curPos)
    if matched_pos_list[0] == -1:
        idx = 0
        for i in range(len(matched_pos_list)):
            if matched_pos_list[i] != -1:
                idx = i
                break 
        matched_pos = matched_pos_list[idx]
        idxList = range(idx)
        for i in idxList[::-1]:
            matched_pos = matched_pos + 1
            matched_pos_list[i] = matched_pos
    if matched_pos_list[-1] == -1:
        idxList = range(len(matched_pos_list))
        idx = 0
        for i in idxList[::-1]:
            if matched_pos_list[i] != -1:
                idx = i
                break
        matched_pos =  matched_pos_list[idx]
        idx = idx + 1
        for i in range(idx, len(matched_pos_list)):
            matched_pos = matched_pos - 1
            matched_pos_list[i] = matched_pos
    if matched_pos_list[0] < 0:
        matched_pos_list[0] = 0
    if matched_pos_list[0] > rna_length:
        matched_pos_list[0] = rna_length - 1
    if matched_pos_list[-1] < 0:
        matched_pos_list[-1] = 0
    if matched_pos_list[-1] > rna_length:
        matched_pos_list[-1] = rna_length - 1
    return([matched_pos_list[0], matched_pos_list[-1]])

def get_star(premir_seq, premir_str, mature_seq):
    matched_str = parse_rna_str(premir_str)
    [mature_pos_start, mature_pos_end] = pos_in_pre(mature_seq, premir_seq)
    [matched_start, matched_end] = get_matched_pos(mature_pos_start, mature_pos_end, matched_str)
    matched_start = matched_start + 2
    matched_end = matched_end + 2
    rna_length = len(premir_seq)
    if matched_start < 0:
        matched_start = 0
    if matched_end < 0:
        matched_end = 0
    if matched_start > rna_length:
        matched_start = rna_length - 1
    if matched_end > rna_length:
        matched_end = rna_length - 1
    star_seq = ""
    if matched_start < matched_end:
        star_seq = premir_seq[matched_start: (matched_end+1)]
        star_seq = star_seq[::-1]
    else:
        star_seq = premir_seq[matched_end: (matched_start+1)]
    return(star_seq)    

def get_sub_str(sub_seq, premir_seq, premir_str):
    [startPos, endPos] = pos_in_pre(sub_seq, premir_seq)
    mir_sub_str = premir_str[startPos:endPos]
    return mir_sub_str

def getMirnaIncludedInLoop(premir_seq, premir_str, sub_seq):
    flag = False
    mir_sub_str = get_sub_str(sub_seq, premir_seq, premir_str)
    if(('(' in mir_sub_str) and (')' in mir_sub_str)):
        flag = True
    return flag
#end fun


###check miRNA located in which arm
def check_arm(premir_seq, premir_str, sub_seq):
    armDetailedType = ""    
    if(getMirnaIncludedInLoop(premir_seq, premir_str, sub_seq)):
        armDetailedType = "loop"
    #check arms
    mir_sub_str = get_sub_str(sub_seq, premir_seq, premir_str)
    armDetailedType = "unmatchedRegion"
    if('(' in mir_sub_str):
        armDetailedType = "arm5"
    if(')' in mir_sub_str):
        armDetailedType = "arm3"
    return armDetailedType

def get_complementarity_seq(premir_seq, premir_str, sub_seq, matchedStructList):
    [mirnaStart, mirnaEnd] = pos_in_pre(sub_seq, premir_seq)
    [matchedPosStart, matchedPosEnd] = get_matched_pos(mirnaStart, mirnaEnd, matchedStructList)
    if(matchedPosStart == 0 and matchedPosEnd == 0):
        return 'N'*len(sub_seq)
    if(matchedPosStart < matchedPosEnd):
        return premir_seq[matchedPosStart: (matchedPosEnd+1)]
    else:
        subStr = premir_seq[matchedPosEnd: (matchedPosStart+1)]
        return subStr[::-1]
    size = getMiRNALen(sub_seq)
    [mirnaStart, mirnaEnd] = pos_in_pre(sub_seq, premir_seq)
    armType = check_arm(premir_seq, premir_str, sub_seq)
    miRNAStar = ""
    if(armType == "arm5"):
        startToMirna = premir_str[0:mirnaStart]
        openPar = startToMirna.count('(') - startToMirna.count(')')
        #analysing from the end of precursor if presence of loop
        openParEnd=0
        i=len(premir_str)-1
        while (openPar != openParEnd):
            if(premir_str[i]== ')'):
                openParEnd = openParEnd + 1
            if(premir_str[i]== '('):
                openParEnd = openParEnd - 1
            if(openPar != openParEnd):
                i = i - 1
        star=""
        posOnStar = i
        i = mirnaStart
        while(i >= mirnaStart and i < mirnaEnd):
            posOnStar = posOnStar - 1
            a = premir_str[i]
            b = premir_str[posOnStar]
            #for bulge before miRNA start position
            if(i == mirnaStart and bracket_or_not(a) == True and bracket_or_not(b) == False):
                continue
            if(bracket_or_not(a) == True and bracket_or_not(b) == True or a==b):
                star = star + premir_seq[posOnStar]
            elif(a =='.' and bracket_or_not(b) == True):
                posOnStar = posOnStar + 1
            elif(b =='.' and bracket_or_not(a) == True):
                star = star + premir_seq[posOnStar]
                i = i - 1
            i = i + 1
        miRNAStar = star[::-1]
    elif(armType == "arm3"):
        EndToMirna = premir_str[mirnaStart+size: len(premir_str)]
        #get open parenthesis minus close parenthesis between precursor end and mirna end
        closePar= EndToMirna.count(')') - EndToMirna.count('(');
        #analysing from the end of precursor if presence of loop
        openPar = 0
        i = 0
        while(openPar != closePar):
            if(premir_str[i] == '('):
                openPar = openPar + 1
            if(premir_str[i] == ')'):
                openPar = openPar - 1
            if(openPar != closePar):
                i = i + 1
        #get complement            
        star=""
        posOnStar = i
        #print premir_str[posOnStar + 1: posOnStar+ 31]
        idxVec = range(mirnaStart, mirnaEnd)
        i = mirnaEnd - 1
        while(i >= mirnaStart and i < mirnaEnd):
            posOnStar = posOnStar + 1
            a = premir_str[i]
            b = premir_str[posOnStar]
            #for a bulge after the miRNA end position
            if(bracket_or_not(a) and i == (mirnaEnd - 1) and bracket_or_not(b) == False):
                continue
            
            if(bracket_or_not(a) and bracket_or_not(b) or a==b):
                star = star + premir_seq[posOnStar]
            elif(a == '.' and bracket_or_not(b)):
                posOnStar = posOnStar - 1
            elif(b =='.' and bracket_or_not(a)):
                star = star + premir_seq[posOnStar]
                i = i + 1
            i = i - 1
        miRNAStar = star[::-1]
    else:
        miRNAStar = 'N'*len(sub_seq)
    return miRNAStar

def checkOSType():
    from sys import platform as _platform
    if _platform == "linux" or _platform == "linux2":
        return "linux"
    elif _platform == "darwin":
        return "mac"
    elif _platform == "win32":
        return "win"
    else:
        return "undefined OS"

def readlines_from_file(fileDir):
    fpr = open(fileDir, "r")
    lines = fpr.readlines()
    fpr.close()
    return lines

def removeFile(files):
    if(os.path.isfile(files) == True):
        try:
            os.remove(files)
        except:
            print("Warning: Failed to remove ", files)

def get_MFE(sub_seq, miRNAStarSeq, tempFileDir):
    with open('{}test.fa'.format(tempFileDir), 'w') as seq_out:
        seq_out.write('>miRNASeq\n{}\n>miRNAStarSeq\n{}\n'.format(sub_seq, miRNAStarSeq))
    resultFileDir = tempFileDir + "test.txt"
    osType = checkOSType()
    order = ""
    if(osType == "win"):
        order = "RNAduplex.exe < " + tempFileDir + "test.fa > " +  resultFileDir
    elif(osType == "linux"):
        order = "RNAduplex < " + tempFileDir + "test.fa > " +  resultFileDir
    else:
        sys.exit("Error: undefined OS")
    res = os.system(order)
    if(res != 0):
        print("Error to run RNAduplex for calculating MFE for miRNASeq:", sub_seq, " and miRNAStarSeq:", miRNAStarSeq) 
        sys.exit() 
    #get MFE
    lines = readlines_from_file(resultFileDir)
    if(len(lines) != 3):
        sys.exit("Error: line number is not THREE in RNAduplex output file:")
    lastLine = lines[-1]
    lastLine = lastLine.strip()
    items = lastLine.split("(")
    lastItem = items[-1]
    MFE = lastItem.replace(")", "")
    removeFile('{}test.fa'.format(tempFileDir))
    removeFile(resultFileDir)
    return({'sub_MFE': MFE})

def maxlen_without_bulges(mir_sub_str):
    Len = Max = 0
    for i in range(len(mir_sub_str)):
        ch = mir_sub_str[i]
        if ch == '(' or ch == ')':
            Len += 1
        else:
            Len = 0
        if Max < Len:
            Max = Len
    return({'len_without_bulges': Max})

def per_maxlen_without_bulges(mir_sub_str):
    b_val = list(maxlen_without_bulges(mir_sub_str).values())[0]
    pc_val = round(1.0*b_val/len(mir_sub_str), 4)
    return({'percentage_len_without_bulges': pc_val})

def basepair_in_duplex(mir_sub_str):
    p = mir_sub_str.count('(') + mir_sub_str.count(')')
    return({'basepair_induplex': p})

def dist_from_terminal_loop(sub_seq, premir_seq, premir_str):
    arm_check = check_arm(premir_seq, premir_str, sub_seq)
    [sub_pos_start, sub_pos_end] = pos_in_pre(sub_seq, premir_seq)
    if(arm_check == "arm5"):
        leftRNAStruct = premir_str[sub_pos_start:len(premir_str)]
        closestLoopEndPar = leftRNAStruct.find(')') + sub_pos_start
        leftRNAStruct = premir_str[sub_pos_start:closestLoopEndPar]
        loopStart = leftRNAStruct.rfind('(') + sub_pos_start
        res = loopStart - (sub_pos_start + len(sub_seq) - 1)
    elif(arm_check == "arm3"):
        leftRNAStruct = premir_str[0: sub_pos_start]
        closestLoopOpenPar = leftRNAStruct.rfind('(')
        leftRNAStruct = premir_str[closestLoopOpenPar: sub_pos_start]
        loopEnd = leftRNAStruct.find(')') + closestLoopOpenPar
        res = sub_pos_start - loopEnd
    else:
        res = -1
    return({'dist_from_terminal_loop': res})

def dist_from_hairpin(sub_seq, premir_seq, premir_str):
    arm_check = check_arm(premir_seq, premir_str, sub_seq)
    [sub_pos_start, sub_pos_end] = pos_in_pre(sub_seq, premir_seq)
    tmp_val = 0
    if arm_check == "arm5":
        tmp_val = sub_pos_start
    elif arm_check == "arm3":
        tmp_val = len(premir_str) - sub_pos_start - len(sub_seq)
    else:
        if ')' not in premir_str:
            tmp_val = sub_pos_start
        elif '(' not in premir_str:
            tmp_val = len(premir_str) - sub_pos_start - len(sub_seq)
        else:
            tmp_val = sub_pos_start
    return({'dist_from_hairpin': tmp_val})

def overlap_in_loop(sub_seq, premir_seq, premir_str):
    t1 = list(dist_from_terminal_loop(sub_seq, premir_seq, premir_str).values())[0]
    tmp_val = 0
    if(t1 < 0):
        tmp_val =  -1.0*t1
    return({'overlap_in_loop': tmp_val})

def count_bulges(sub_seq, premir_seq, premir_str):
    mir_sub_str = get_sub_str(sub_seq, premir_seq, premir_str)
    num = 0
    for i in range(len(mir_sub_str)-2):
        ch1 = mir_sub_str[i]
        ch2 = mir_sub_str[i+1]
        if(ch1 == '.' and bracket_or_not(ch2)):
            num = num + 1
    return({'bulge number': num})

def count_perfect_20mer(mir_sub_str):
    if mir_sub_str.find("((((((((((((((((((((") >= 0 or mir_sub_str.find(")))))))))))))))))))))") >= 0:
        return({'perfect_20mer': 1})
    else:
        return({'perfect_20mer': 0})

##start of 20mer base paired
def start_perfect_20mer(mir_sub_str):
    tmp_val = list(count_perfect_20mer(mir_sub_str).values())[0]
    if tmp_val == 1:
        s = mir_sub_str.find("((((((((((((((((((((")
        if s == -1:
            s = mir_sub_str.find(")))))))))))))))))))))")
        return({'start_perfect_20mer': s+1})
    else:
        return({'start_perfect_20mer': -1})

##10mer base paired
def count_perfect_10mer(mir_sub_str):
    if mir_sub_str.find("((((((((((") >= 0 or mir_sub_str.find("))))))))))") >= 0:
        return({'perfect_10mer': 1})
    else:
        return({'perfect_10mer': 0})

##start of 10mer base paired
def start_perfect_10mer(mir_sub_str):
    tmp_val = list(count_perfect_10mer(mir_sub_str).values())[0]
    if tmp_val == 1:
        s = mir_sub_str.find("((((((((((")
        if s == -1:
            s = mir_sub_str.find("))))))))))")
        return({'start_perfect_10mer': s+1})
    else:
        return({'start_perfect_10mer': -1})

##5mer base paired
def count_perfect_5mer(mir_sub_str):
    if mir_sub_str.find("(((((") >= 0 or mir_sub_str.find(")))))") >= 0:
        return({'perfect_5mer': 1})
    else:
        return({'perfect_5mer': 0})

##start of 5mer base paired
def start_perfect_5mer(mir_sub_str):
    tmp_val = list(count_perfect_5mer(mir_sub_str).values())[0]
    if tmp_val == 1:
        s = mir_sub_str.find("(((((")
        if s == -1:
            s = mir_sub_str.find(")))))")
        return({'start_perfect_5mer': s+1})
    else:
        return({'start_perfect_5mer': -1})

def mean_bp_in_window(mir_sub_str, w):
    numList = []
    total = 0
    for i in range(len(mir_sub_str) - w):
        subStruct = mir_sub_str[i:i+w+1]
        m1 = subStruct.count("(")
        m2 = subStruct.count(")")
        numList.append(m1 + m2)
        total = total + m1
        total = total + m2
    tmp_val = round(1.0*total/len(numList), 4)    
    return({'mean_bp_in_{}window'.format(w): tmp_val})

#Bulge at specific position around the start of miRNA, 0: False; 1:True
def start_pos_is_bulge(miRNA, premir_seq, premir_str, offset):
    [sub_pos_start, sub_pos_end] = pos_in_pre(miRNA, premir_seq)
    spPos = sub_pos_start + offset
    tmp_val = 0
    if(spPos < 0 or spPos > len(premir_str)-1):
        tmp_val = 0
    else:
        ch = premir_str[spPos]
        if(ch == '.'):
            tmp_val = 1
        else:
            tmp_val = 0
    return({'start_pos{}_is_bulge'.format(offset): tmp_val})

#Bulge at end position around the end of miRNA. 0:False; 1:True
def end_pos_is_bulge(sub_seq, premir_seq, premir_str, offset):
    [sub_pos_start, sub_pos_end]  = pos_in_pre(sub_seq, premir_seq)
    spPos = sub_pos_end + offset
    tmp_val = 0
    if(spPos < 0 or spPos >= len(premir_str)-1):
        tmp_val = 0
    else:
        ch = premir_str[spPos]
        if(ch == '.'):
            tmp_val = 1
        else:
            tmp_val = 0
    return({'end_pos{}_is_bulge'.format(offset): tmp_val})    

##get all mono-Nucleotide&Structure
def get_monoseq_monostr(sub_seq, premir_seq, premir_str):
    mir_sub_str = get_sub_str(sub_seq, premir_seq, premir_str)
    miRNALen = len(sub_seq)
    tmp_dict = {"A.":0, "A(":0, "A)":0,
                "C.":0, "C(":0, "C)":0, 
                "G.":0, "G(":0, "G)":0, 
                "U.":0, "U(":0, "U)":0, 
                "mono_other":0}
    for i in range(miRNALen):
        cur_key = sub_seq[i] + mir_sub_str[i]
        if cur_key in tmp_dict:
            tmp_dict[cur_key] += 1
        else:
            tmp_dict["mono_other"] += 1
    res_dict = {}        
    for cur_key in sorted(tmp_dict.keys()):
        res_dict['sub_' + cur_key] = round(1.0*tmp_dict[cur_key]/miRNALen, 4)
    return(res_dict)

def get_diseq_distr(sub_seq, premir_seq, premir_str):
    mir_sub_str = get_sub_str (sub_seq, premir_seq, premir_str)
    chars = ["A", "C", "G", "U"]
    dots = ["(", ".", ")"]
    tmp_dict = {}
    for i in range(len(chars)):
        for j in range(len(chars)):
            for x in range(len(dots)):
                for y in range(len(dots)):
                    cur_key = chars[i] + chars[j] + dots[x] + dots[y]
                    tmp_dict[cur_key] = 0
    tmp_dict["di_other"] = 0
    for i in range(1, len(sub_seq)):
        ch1 = sub_seq[i-1]
        ch2 = sub_seq[i]
        dot1 = mir_sub_str[i-1]
        dot2 = mir_sub_str[i]
        cur_key = ch1 + ch2 + dot1 + dot2
        if(cur_key in tmp_dict.keys()):
            tmp_dict[cur_key] += 1
        else:
            tmp_dict["di_other"] += 1
    res_dict = {}         
    for cur_key in sorted(tmp_dict.keys()):
        res_dict['sub_' + cur_key] = round(1.0*tmp_dict[cur_key]/len(sub_seq), 4)
    return(res_dict)            

##get all triplets
def get_monoseq_tristr(sub_seq, premir_seq, premir_str):
    mir_sub_str = get_sub_str(sub_seq, premir_seq, premir_str)
    tmp_dict = {"A...":0,"C...":0,"G...":0,"U...":0,
              "A(..":0,"C(..":0,"G(..":0,"U(..":0,
              "A.(.":0,"C.(.":0,"G.(.":0,"U.(.":0,
              "A..(":0,"C..(":0,"G..(":0,"U..(":0,
              "A((.":0,"C((.":0,"G((.":0,"U((.":0,
              "A(.(":0,"C(.(":0,"G(.(":0,"U(.(":0,
              "A.((":0,"C.((":0,"G.((":0,"U.((":0,
              "A(((":0,"C(((":0,"G(((":0,"U(((":0,
              "tri_other":0}
    for i in range(1,len(sub_seq)-2):
        cur_key = sub_seq[i] + mir_sub_str[i-1:i+2]
        if cur_key in tmp_dict:
            tmp_dict[cur_key] += 1
        else:
            tmp_dict["tri_other"] += 1
    res_dict = {}         
    for cur_key in sorted(tmp_dict.keys()):
        res_dict['sub_' + cur_key] = round(1.0*tmp_dict[cur_key]/len(sub_seq), 4)
    return(res_dict)

def getDPPSFileDir (dpSSFileDic, curPreMiRNAID):
    dpPSFile = dpSSFileDic + curPreMiRNAID + "_dp.ps"
    if(os.path.exists(dpPSFile)):
        return dpPSFile
    dpPSFile = dpSSFileDic + curPreMiRNAID + "-5p_dp.ps"
    if(os.path.exists(dpPSFile)):
        return dpPSFile
    dpPSFile = dpSSFileDic + curPreMiRNAID + "-3p_dp.ps"
    if(os.path.exists(dpPSFile)):
        return dpPSFile 
    idx = curPreMiRNAID.rfind("-")
    dpPSFile = dpSSFileDic + curPreMiRNAID[:idx] + "_dp.ps"
    if(os.path.exists(dpPSFile)):
        return dpPSFile
    return "NotFind"

##decoding *_dp.ps file for get positional Entropy
def get_seq_entropy(dpPSFile, premir_seq):
    length = len(premir_seq)
    mm = {}
    mp = {}
    pp = {}
    sp = {}
    p0 = {}
    for i in range(0, length+1):
        mm[i] = 0.0
        mp[i] = 0.0
        pp[i] = 0.0
        sp[i] = 0.0
        p0[i] = 0.0
    max = 0.0
    lines = readlines_from_file(dpPSFile)
    flag = 0
    for i in range(len(lines)):
        curLine = lines[i]
        if("%start of base pair probability data" in curLine):
            flag = 1
            continue
        if("showpage" in curLine):
            flag = 0
        if(flag != 1):
            continue
        curLine = curLine.strip()
        [posi, posj, prob, boxType] = curLine.split(" ")
        #print posi, posj, prob, boxType, length
        posi = int(posi)
        posj = int(posj)
        prob = float(prob)
        if(boxType == "ubox"):
            #square prob to probability
            prob = prob*prob
            ss = 0
            if(prob > 0):
                ss = prob*math.log(prob)
            mp[posi + 1] = mp[posi+1] + prob
            mp[posj] = mp[posj] - prob
            sp[posi] = sp[posi] + ss
            sp[posj] = sp[posj] + ss
            pp[posi] = pp[posi] + prob
            pp[posj] = pp[posj] + prob
        elif(boxType == "lbox"):
            mm[posi+1] = mm[posi+1] + 1
            mm[posj] = mm[posj] - 1
    mp[0] = 0
    mm[0] = 0
    max = 0
    posProb = []
    posMFE = []
    posEntropy = []
    for i in range(1,length+1):
        mp[i] = mp[i] + mp[i-1]
        if(mp[i] > max):
           max = mp[i]
        mm[i] = mm[i] + mm[i-1]
        sp[i] = sp[i] + (1-pp[i])*math.log(1-pp[i])
        tmp = -1.0*sp[i]/math.log(2)
        posEntropy.append(tmp)
        posProb.append(mp[i])
        posMFE.append(mm[i])
    removeFile(dpPSFile)
    return [posProb, posMFE, posEntropy]

#get positional entropies around a specific positions.
def entropy_in_pos(posEntropy, position, w, l, tmp_tp):
    upList = []
    downList = []
    p1 = position - w
    p2 = position + l + 1
    if p1 < 0:
        upList = [posEntropy[0]]*abs(p1)
        p1 = 0
    if p2 > len(posEntropy):
        downList = [posEntropy[-1]]*abs(p2 - len(posEntropy))
        p2 = len(posEntropy)
    res = upList + posEntropy[p1:p2] + downList
    ## output
    res_dict = {}
    tmp_list = list(range(-w, l+1))
    for ii in range(len(res)):
        res_dict['sub_{}_{}_entropy'.format(tmp_tp, tmp_list[ii])] = round(res[ii], 6)
    return(res_dict)

#get nucleotide type at specific position
def codeDiNucleotideForPosition(premir_seq, Position1, Position2):
    RNASeqLen = len(premir_seq)
    pos1 = Position1
    pos2 = Position2
    if(pos1 < 0):
        pos1 = 0
    if(pos2 < 0):
        pos2 = 1
    if(pos1 > RNASeqLen - 1):
        pos1 = RNASeqLen - 2
    if(pos2 > RNASeqLen - 1):
        pos2 = RNASeqLen - 1
    return (codeDiNucleotide(premir_seq[pos1], premir_seq[pos2]))

#code nucleotide type
def codeNucleotide(ch):
    dist = {"A":[0,0,0,1], "C":[0,0,1,0], "G":[0,1,0,0], "U":[1,0,0,0], 
            "a":[0,0,0,1], "c":[0,0,1,0], "g":[0,1,0,0], "u":[1,0,0,0]}
    if ch in "ACGUacgu":
        return dist[ch]
    else:
        return [0,0,0,0]

##get nucleotide type at specific position
def codeNucleotideForPosition(premir_seq, Position):
    pos = Position
    if(pos < 0):
        pos = 0
    if(pos > len(premir_seq) - 1):
        pos = len(premir_seq) - 1
    ch = premir_seq[pos]
    return(codeNucleotide(ch))

def decode_pos_up_down(premir_seq, Position, up, down, tmp_tp):
    offsetList = range(-1*up, down+1)
    res_dict = {}
    for curOffset in offsetList:
        res_dict['sub_{}_{}_decode'.format(tmp_tp, curOffset)] = codeNucleotideForPosition(premir_seq, Position + curOffset)
    return(res_dict)

def featureExtract(mature_seq, premir_seq, premir_str, curr_dir, premir_id, 
                   upOffLen = 5, downOffLen = 5):
    mature_seq = re.sub('T', 'U', mature_seq)
    premir_seq = re.sub('T', 'U', premir_seq)
    matched_str = parse_rna_str(premir_str)
    feature_dict = {}
    #1st feature:miRNA length
    mature_len = len(mature_seq)
    #2nd feature: miRNA GCContent
    mature_gc = gc_content(mature_seq)
    feature_dict.update({'sub_Len': mature_len, 'sub_GC':mature_gc})
    #3-6 feature: content of mono nucleotide 
    feature_dict.update(mononucleotide_ratio(mature_seq))
    #7-22 features: dinucleotide content 
    feature_dict.update(dinucleotide_ratio(mature_seq))
    #miRNA star 
    star_seq = get_star(premir_seq, premir_str, mature_seq)
    [star_start_pos, star_end_pos] = pos_in_pre(star_seq, premir_seq)
    cs_seq = get_complementarity_seq(premir_seq, premir_str, mature_seq, matched_str)    
    #23rd feature: MFE
    feature_dict.update(get_MFE(mature_seq, cs_seq, curr_dir))
    #get miRNA structure
    mir_sub_str = get_sub_str(mature_seq, premir_seq, premir_str)
    if len(mir_sub_str) < 16 or len(mir_sub_str) > 32:
        print(premir_id)
    #feature: maximal length of miRNA without bulge
    feature_dict.update(maxlen_without_bulges(mir_sub_str))
    #feature: percenate of maximal length of miRNA without bulge
    feature_dict.update(per_maxlen_without_bulges(mir_sub_str))
    #feature: number of base pairs in miRNA and miRNA* duplex
    feature_dict.update(basepair_in_duplex(mir_sub_str))
    #feautre:distance to the terminal loop 
    feature_dict.update(dist_from_terminal_loop(mature_seq, premir_seq, premir_str))
    #feature: distance to the start of helix
    feature_dict.update(dist_from_hairpin(mature_seq, premir_seq, premir_str))
    #feature: number of bases in loop 
    feature_dict.update(overlap_in_loop(mature_seq, premir_seq, premir_str))
    #1 feature: number of bulges in miRNA sequence
    feature_dict.update(count_bulges(mature_seq, premir_seq, premir_str))
    #feature: perfect base pair
    feature_dict.update(count_perfect_20mer(mir_sub_str))
    feature_dict.update(start_perfect_20mer(mir_sub_str))
    feature_dict.update(count_perfect_10mer(mir_sub_str))
    feature_dict.update(start_perfect_10mer(mir_sub_str))
    feature_dict.update(count_perfect_5mer(mir_sub_str))
    feature_dict.update(start_perfect_5mer(mir_sub_str))
    #feature: base pair
    feature_dict.update(mean_bp_in_window(mir_sub_str, 7))
    feature_dict.update(mean_bp_in_window(mir_sub_str, 5))
    feature_dict.update(mean_bp_in_window(mir_sub_str, 3))
    #feature: bulge at specific position around the start of miRNA sequence 
    feature_dict.update(start_pos_is_bulge(mature_seq, premir_seq, premir_str, -3))
    feature_dict.update(start_pos_is_bulge(mature_seq, premir_seq, premir_str, -2))
    feature_dict.update(start_pos_is_bulge(mature_seq, premir_seq, premir_str, -1))
    feature_dict.update(start_pos_is_bulge(mature_seq, premir_seq, premir_str, 0))
    feature_dict.update(start_pos_is_bulge(mature_seq, premir_seq, premir_str, 1))
    feature_dict.update(start_pos_is_bulge(mature_seq, premir_seq, premir_str, 2))
    feature_dict.update(start_pos_is_bulge(mature_seq, premir_seq, premir_str, 3))
    #feature: bulge at specific position around the end of miRNA sequence 
    feature_dict.update(end_pos_is_bulge(mature_seq, premir_seq, premir_str, -3))
    feature_dict.update(end_pos_is_bulge(mature_seq, premir_seq, premir_str, -2))
    feature_dict.update(end_pos_is_bulge(mature_seq, premir_seq, premir_str, -1))
    feature_dict.update(end_pos_is_bulge(mature_seq, premir_seq, premir_str, 0))
    feature_dict.update(end_pos_is_bulge(mature_seq, premir_seq, premir_str, 1))
    feature_dict.update(end_pos_is_bulge(mature_seq, premir_seq, premir_str, 2))
    feature_dict.update(end_pos_is_bulge(mature_seq, premir_seq, premir_str, 3))
    #12 features: mono-nucleotide & strucutre features
    feature_dict.update(get_monoseq_monostr(mature_seq, premir_seq, premir_str))
    # features: di-nucleotide and structure features
    feature_dict.update(get_diseq_distr(mature_seq, premir_seq, premir_str))
    #33 features: triplets sequence & strucutre features
    feature_dict.update(get_monoseq_tristr(mature_seq, premir_seq, premir_str))
    #19*2 features: positional entropy features
    [mature_pos_start, mature_pos_end] = pos_in_pre(mature_seq, premir_seq)
    # previous use 5
    dpPSFile = getDPPSFileDir(curr_dir,  premir_id)
    [curPosProb, curPosMFE, RNAPosEntropy] = get_seq_entropy(dpPSFile, premir_seq) 
    feature_dict.update(entropy_in_pos(RNAPosEntropy, mature_pos_start, upOffLen, downOffLen, 'ms'))
    feature_dict.update(entropy_in_pos(RNAPosEntropy, mature_pos_end, upOffLen, downOffLen, 'me'))
    #nucleotide type at position around the start and end of miRNA sequence
    feature_dict.update(decode_pos_up_down(premir_seq, mature_pos_start, upOffLen, downOffLen, 'ms'))
    feature_dict.update(decode_pos_up_down(premir_seq, mature_pos_end, upOffLen, downOffLen, 'me'))
    #for miRNA start, region around the starts and ends.
    feature_dict.update(decode_pos_up_down(premir_seq, star_start_pos, upOffLen, downOffLen, 'ss'))
    feature_dict.update(decode_pos_up_down(premir_seq, star_end_pos, upOffLen, downOffLen, 'se'))
    # return all features
    return(feature_dict)

## input data
input_file = '../out_pool_merge.txt'
curr_dir = './'
ii_num = 1
mirplant = ['ID', 'cg', '%AA', '%AC', '%AG', '%AU', '%CA', '%CC', '%CG', '%CU', '%GA', 
'%GC', '%GG', '%GU', '%UA', '%UC', '%UG', '%UU', 'mfe1', 'mfe2', 'dG',
'dP', 'dQ', 'dD', 'dF', 'zG', 'zP', 'zQ', 'zD', 'zF', 'mfe3',
'mfe4', 'nefe', 'freq', 'div', 'diff', 'dH', 'dHL', 'dS', 'dSL', 'Tm',
'TmL', 'auL', 'gcL', 'guL', 'bpStems', 'auStems', 'gcStems', 'guStems',
'mfe5', 'mfe6', 'avg_mis_num', 'mfe7', 'mfe8', 'mfe9', 'mis_num_begin',
'mis_num_end', ['triple_begin']*32, ['triple_end']*32, ['triple']*32]
header_name = []
for tmp_name in mirplant:
    if type(tmp_name) is list:
        for tmp_name_each in tmp_name:
            header_name.append(tmp_name_each)
    else:
        header_name.append(tmp_name)

ii = 1
with open('seq.fasta', 'w') as tmp_write:
    with open(input_file) as file_in:
        for eli in file_in.readlines()[1:]:
            eli = eli.strip().split('\t')
            res_seq = re.sub('N', '', eli[2])
            tmp_write.write('>pre_{}\n{}\n'.format(ii, res_seq))
            ii += 1
# perl script
os.system('perl extract_features.pl')

# perl results
ii_num = 1
feature_dict = {}
with open("../data/seq.svm.data") as sub_feat:
    for eli in sub_feat.readlines():
        eli = eli.strip().split(' ')
        feature_dict[ii_num] = [x.split(':')[1] for x in eli[1:]]
        ii_num += 1

print(len(feature_dict))
ii_num = 1 
with open('../00feature_out.txt', 'w') as feat_out:
    with open( input_file ) as input_rna:
        for eli in tqdm(input_rna.readlines()[1:]):
            eli = eli.strip().split('\t')
            a,b = subprocess.getstatusoutput('echo {}| RNAfold --id-prefix={} -p --noLP'.format(
            eli[2], ii_num))
            a,process_res = subprocess.getstatusoutput('echo {}| RNAfold --noPS'.format(eli[2]))
            rnafold_res = re.split("\n|\s\(", process_res)
            if '(' in rnafold_res[1]:
                fe_out = featureExtract(eli[7], eli[2], rnafold_res[1], curr_dir, str(ii_num)+'_0001')
                for tmp_key in fe_out:
                    if type(fe_out[tmp_key]) is list:
                        for tmp_val in fe_out[tmp_key]:
                            feature_dict[ii_num].append(str(tmp_val))
                            if ii_num == 1:
                                header_name.append(tmp_key)
                    else:
                        feature_dict[ii_num].append(str(fe_out[tmp_key]))
                        if ii_num == 1:
                            header_name.append(tmp_key)
                if ii_num == 1:
                    feat_out.write('\t'.join(header_name) + '\n')
                feat_out.write('{}\t{}\n'.format(eli[1], '\t'.join(feature_dict[ii_num])))
                ii_num += 1
            else:
                print(eli, process_res, rnafold_res)  

os.system('rm *_ss.ps')
