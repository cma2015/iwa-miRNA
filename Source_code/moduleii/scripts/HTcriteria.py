# Load the package
import re
import os
import sys
import time
import random
import subprocess
import pandas as pd
import multiprocessing

def tag_id_abu(fa_file):
    '''Extract sequence and abundance from fasta file'''
    global tagId_seq
    global tagId_count
    tagId_seq, tagId_count = {}, {}
    with open(fa_file) as ffread:
        for eli in ffread.readlines():
            if eli.startswith('>'):
                eli = eli.strip().split(' ')
                eli[0] = eli[0].strip('>')
                seq_id = str(eli[0])
                seq_count = round(float(eli[1]),3)
            else:
                tagId_seq[seq_id] = eli.strip()
                tagId_count[eli.strip()] = seq_count

def extract_map_read(read_map_file ,input_dict):
    df_sran = pd.read_csv(read_map_file, sep='\t', header=None,
    names=['ID','chr','start','end','strand','TPM','SNP'])
    separate_key_df = {}
    for ikey in input_dict:
        pre_location = ikey.split(':')
        pre_start = int(pre_location[1])
        pre_end = int(pre_location[2])
        sub_frame_index = ((df_sran['chr'] == pre_location[0]) &
        (df_sran['start'] >= pre_start-3) &
        (df_sran['end'] <= pre_end+3))
        ## extract reads from file
        sub_frame = df_sran[sub_frame_index]
        sub_frame = sub_frame[['ID','start','end','strand','TPM','SNP']]
        sub_frame = sub_frame.sort_values(by=['start','end'],axis=0, ascending=True)
        separate_key_df[ikey] = sub_frame
    return separate_key_df

def centroid_fold(input_dict):
    get_centroid_dict = {}
    with open("{}/tmp.fasta".format(out_path), "w") as fas_file:
        for rna_site in input_dict:
            fas_file.write('>%s\n%s\n' % (rna_site, input_dict[rna_site][2]))
    fas_str = 'centroid_fold -e CONTRAfold -g 4 {0}/tmp.fasta > {0}/tmp.structure'.format(out_path)
    other, res = subprocess.getstatusoutput(fas_str)
    with open("{}/tmp.structure".format(out_path), "r") as str_file:
        for i in str_file.readlines():
            if i[:1] == '>':
                name=i.replace('>','').split()[0]
                get_centroid_dict[name]=[]
            elif i[:1] in [".", "("]:
                get_centroid_dict[name]=[i.split(" ")[0]]
    os.system('rm {}/tmp.fasta {}/tmp.structure'.format(out_path, out_path))
    return get_centroid_dict

def rnafold_fold(input_dict):
    rnafold_str = {}
    for rna_site in input_dict:
        a,process_res = subprocess.getstatusoutput('echo {}| RNAfold --noPS '.format(
        input_dict[rna_site][2]))
        rnafold_res = re.split("\n|\s\(", process_res)
        mfe = re.sub(r'\)|\s', '', rnafold_res[2])
        amfe = float(mfe)/len(input_dict[rna_site][2])*100
        rnafold_str[rna_site] = [rnafold_res[1], mfe, amfe]
    return(rnafold_str)

def iso_type(ref_location,query_location):
    type_out = '-'
    if -3 <= query_location[0]-ref_location[0] < 0:
        type_out = 'add5'
        if query_location[1] > ref_location[1]:
            type_out = 'add5_add3'
        elif query_location[1] < ref_location[1]:
            type_out = 'add5_sub3'
    elif 0 < query_location[0]-ref_location[0] <= 3:
        type_out = 'sub5'
        if query_location[1] > ref_location[1]:
            type_out = 'sub5_add3'
        elif query_location[1] < ref_location[1]:
            type_out = 'sub5_sub3'
    elif 0 < query_location[1]-ref_location[1] <= 3:
        type_out = 'add3'
        if query_location[0] < ref_location[0]:
            type_out = 'add3_add5'
        elif query_location[0] > ref_location[0]:
            type_out = 'add3_sub5'
    elif -3 <= query_location[1]-ref_location[1] < 0:
        type_out = 'sub3'
        if query_location[0] < ref_location[0]:
            type_out = 'sub3_add5'
        elif query_location[0] > ref_location[0]:
            type_out = 'sub3_sub5'
    elif query_location[0] == ref_location[0] and query_location[1] == ref_location[1]:
        type_out = 'ref'
    else:
        type_out = '-'
    return type_out

def iso_snp_type(ref_location,query_location,tmp_mutant):
    type_out = '-'
    seed_dist = int(tmp_mutant) + query_location[0] - ref_location[0]
    if query_location[1] - ref_location[1] == 1 and seed_dist == ref_location[1]-ref_location[0]+2:
        type_out = 'nt_add'
    if query_location[0] == ref_location[0] and query_location[1] == ref_location[1]:
        type_out = 'seed_snp'
        if 9 <= seed_dist <= (ref_location[1]-ref_location[0]+1):
            type_out = 'tail_snp'
    return type_out

def sub_in_pre(pre_start, pre_end, mature_start, mature_end, pre_strand):
    if pre_strand == '+':
        return [mature_start-pre_start+1, mature_end-pre_start+1]
    else:
        return [pre_end-mature_end+1, pre_end-mature_start+1]

def bias_get(input_dict, iso_df):
    with open( '{}/{}'.format(out_path, "corr_pre-mirna_seq.txt"), 'w') as out_tmp:
        bias_str = {}
        for rna_site in input_dict:
            info_list = input_dict[rna_site]
            pre_location = rna_site.split(':')
            pre_start = int(pre_location[1])
            pre_end = int(pre_location[2])
            pre_strand = pre_location[3]
            sub_frame = iso_df[rna_site]
            match_num = sub_frame.shape[0]
            if not sub_frame.empty:
                with open('{}/{}.out'.format(out_data, rna_site), 'w') as out_iso:
                    strand_index = (sub_frame['strand'] == pre_strand)
                    all_sequence = [tagId_seq[i] for i in sub_frame[strand_index]['ID']]
                    all_rpm = sum([tagId_count[i] for i in all_sequence])
                    mature = info_list[4].split(':')
                    mature_start = int(mature[1])
                    mature_end = int(mature[2])
                    m_mod_add = mature_end-mature_start+2
                    star = info_list[7].split(':')
                    star_start = int(star[1])
                    star_end = int(star[2])
                    s_mod_add = star_end-star_start+2
                    # mature
                    mature_iso_index = (abs(sub_frame['start'] - mature_start) <=3) & \
                    (abs(sub_frame['end'] - mature_end) <=3) & \
                    (sub_frame['strand'] == pre_strand)
                    mature_sub_df = sub_frame[mature_iso_index]
                    mature_iso_seq = [tagId_seq[i] for i in mature_sub_df['ID']]
                    mature_iso_rpm = sum([tagId_count[i] for i in mature_iso_seq])
                    ##
                    if len(mature_sub_df) > 0:
                        for index in mature_sub_df.index:
                            isomiR_site = [int(mature_sub_df.loc[index, 'start']), int(mature_sub_df.loc[index, 'end'])]
                            isomiR_site = sub_in_pre(pre_start, pre_end, isomiR_site[0], isomiR_site[1], pre_strand)
                            mir_site = sub_in_pre(pre_start, pre_end, mature_start, mature_end, pre_strand)
                            if mature_sub_df.loc[index, 'SNP'] == '-':
                                tmp_type = iso_type(mir_site, isomiR_site)
                            else:
                                snp_info = mature_sub_df.loc[index, 'SNP'].split(':')
                                tmp_type = iso_snp_type(mir_site, isomiR_site, snp_info[0])
                            df_index_list = [str(x) for x in list(mature_sub_df.loc[index, ])]
                            out_iso.write('{}\t5p\t{}\t{}\t{}\n'.format(rna_site,tagId_seq[df_index_list[0]], '\t'.join(df_index_list), tmp_type))
                    # star
                    star_iso_index = (abs(sub_frame['start'] - star_start) <=3) & \
                    (abs(sub_frame['end'] - star_end) <=3) & \
                    (sub_frame['strand'] == pre_strand)
                    star_sub_df = sub_frame[star_iso_index]
                    star_iso_seq = [tagId_seq[i] for i in sub_frame[star_iso_index]['ID']]
                    star_iso_rpm = sum([tagId_count[i] for i in star_iso_seq])
                    ##
                    if len(star_sub_df) > 0:
                        for index in star_sub_df.index:
                            isomiR_site = [int(star_sub_df.loc[index, 'start']), int(star_sub_df.loc[index, 'end'])]
                            isomiR_site = sub_in_pre(pre_start, pre_end, isomiR_site[0], isomiR_site[1], pre_strand)
                            mir_site = sub_in_pre(pre_start, pre_end, star_start, star_end, pre_strand)
                            if star_sub_df.loc[index, 'SNP'] == '-':
                                tmp_type = iso_type(mir_site, isomiR_site)
                            else:
                                snp_info = star_sub_df.loc[index, 'SNP'].split(':')
                                tmp_type = iso_snp_type(mir_site, isomiR_site, snp_info[0])
                            df_index_list = [str(x) for x in list(star_sub_df.loc[index, ])]
                            out_iso.write('{}\t3p\t{}\t{}\t{}\n'.format(rna_site,tagId_seq[df_index_list[0]],
                            '\t'.join(df_index_list), tmp_type))
                    for ii in all_sequence:
                        if ii in mature_iso_seq or ii in star_iso_seq:
                            out_tmp.write('{}\t{}\t{}\t{}\n'.format(rna_site, ii, 'complex', tagId_count[ii]))
                        else:
                            out_tmp.write('{}\t{}\t{}\t{}\n'.format(rna_site, ii, 'all', tagId_count[ii]))
                    ## mature reverse
                    mature_reverse = (abs(sub_frame['start'] - mature_start) <=3) & \
                    (abs(sub_frame['end'] - mature_end) <=3) & \
                    (sub_frame['strand'] != pre_strand)
                    reverse_sequence = [tagId_seq[i] for i in sub_frame[mature_reverse]['ID']]
                    reverse_mature = sum([tagId_count[i] for i in reverse_sequence])
                    ## star reverse
                    mature_reverse = (abs(sub_frame['start'] - star_start) <=3) & \
                    (abs(sub_frame['end'] - star_end) <=3) & \
                    (sub_frame['strand'] != pre_strand)
                    reverse_sequence = [tagId_seq[i] for i in sub_frame[mature_reverse]['ID']]
                    reverse_star = sum([tagId_count[i] for i in reverse_sequence])
                    ## calculate bias
                    arround_rpm = mature_iso_rpm + star_iso_rpm
                    strand_rpm = arround_rpm + reverse_mature + reverse_star
                    if all_rpm == 0:
                        abundance_bias = 0
                    else:
                        abundance_bias = round(arround_rpm/all_rpm, 4)
                    if strand_rpm == 0:
                        strand_bias = 0
                    else:
                        strand_bias = round(arround_rpm/strand_rpm, 4)
                    map_read_count = len(mature_iso_seq) + len(star_iso_seq)
                    if abundance_bias == 0 and strand_bias == 0:
                        bias_str[rna_site] = [round(arround_rpm,2), map_read_count, len(all_sequence),0,0]
                    else:
                        bias_str[rna_site] = [round(arround_rpm,2), map_read_count, len(all_sequence),
                                              abundance_bias, strand_bias]
            else:
                bias_str[rna_site] = [0]*5
        return(bias_str)


def isParenthese(ch):
    if ch=='(' or ch==')':
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
        if len(location) == 0:
            return(['False', '-', '-', 'nonor'])
        else:
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

def mat_star(pre_info):
    if pre_info[10] == '5p':
        mat_seq = pre_info[5]
        star_seq = pre_info[8]
    elif pre_info[10] == '3p':
        mat_seq = pre_info[8]
        star_seq = pre_info[5]
    else:
        if pre_info[5] in tagId_count:
            fp_rpm = tagId_count[pre_info[5]]
        else:
            fp_rpm = 0
        if pre_info[8] in tagId_count:
            tp_rpm = tagId_count[pre_info[8]]
        else:
            tp_rpm = 0
        if fp_rpm >= tp_rpm:
            mat_seq = pre_info[5]
            star_seq = pre_info[8]
        elif fp_rpm < tp_rpm:
            mat_seq = pre_info[8]
            star_seq = pre_info[5]
    return([mat_seq, star_seq])


def format_values_ps(value):
    val_tab = [" %s 0 0.8 0 cfmark" % v for i, v in enumerate(value)]
    return "".join(val_tab)


def structure_plot(info_list):
    pre_file = re.sub(':','_',info_list[1])
    rnafoldstr = '{}/{}_r'.format(out_data, pre_file)
    centroidstr = '{}/{}_c'.format(out_data, pre_file)
    with open(rnafoldstr+'.str', 'w') as str_f:
        str_f.write('>{}_r\n{}\n{}'.format(pre_file,
                                           info_list[2],
                                           rnafold_dict[info_list[1]][0]))
    with open(centroidstr+'.str', 'w') as str_f:
        str_f.write('>{}_c\n{}\n{}'.format(pre_file,
                                           info_list[2],
                                           centroid_dict[info_list[1]][0]))
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
        pre_file, rnafoldstr, centroidstr))
    os.system("gmt psconvert -Tf {0}.ps && pdf2svg {0}.pdf {0}.svg".format(rnafoldstr))
    os.system("gmt psconvert -Tf {0}.ps && pdf2svg {0}.pdf {0}.svg".format(centroidstr))
    os.system("rm {0}.ps {0}.pdf ".format(rnafoldstr))
    os.system("rm {0}.ps {0}.pdf ".format(centroidstr))


def seq_table(seq_info, mature_se, star_se):
    out_raw = '<tr>'
    tmp_seq = seq_info[0]
    if mature_se[0] < star_se[0]:
        seq_td = '<td>{}<span style="background-color: #DE3025">{}</span>{}<span style="background-color: #257ADE">{}</span>{}</td>'.format(
        tmp_seq[0:mature_se[0]], tmp_seq[mature_se[0]:mature_se[1]],tmp_seq[mature_se[1]:star_se[0]],
        tmp_seq[star_se[0]:star_se[1]],tmp_seq[star_se[1]:])
    else:
        seq_td = '<td>{}<span style="background-color: #257ADE">{}</span>{}<span style="background-color: #DE3025">{}</span>{}</td>'.format(
        tmp_seq[0:star_se[0]], tmp_seq[star_se[0]:star_se[1]],tmp_seq[star_se[1]:mature_se[0]],
        tmp_seq[mature_se[0]:mature_se[1]],tmp_seq[mature_se[1]:])
    out_raw += seq_td
    for kk in seq_info[1:]:
        out_raw += '<td >{}</td>'.format(kk)
    out_raw += '</tr>'
    return out_raw

def site_tolow(sequence, site):
    seq_site = sequence[site].lower()
    if site == 0:
        seq_out = seq_site + sequence[(site+1)::]
    else:
        seq_out = sequence[:site] + seq_site + sequence[(site+1):]
    return seq_out

def miR_iso_add(info_list):
    with open('{0}/{1}_map.html'.format(out_data, info_list[1]), 'w') as seq_map:
        iso_df = iso_dict[info_list[1]]
        pre_strand = info_list[1].split(':')[3]
        iso_index = (iso_df['strand'] == pre_strand)
        iso_df = iso_df[iso_index]
        pre_loc = info_list[1].split(':')
        ##
        mature_se = getMiRNAPosition(mat_star(info_list)[0],info_list[2])
        star_se = getMiRNAPosition(mat_star(info_list)[1],info_list[2])
        ##
        seq_map.write('<table style="font-family:monospace">\n')
        seq_map.write('<tr><td></td><td>Type</td><td>TPM</td></tr>\n')
        seq_map.write(seq_table([info_list[2], 'Sequence', '-'], mature_se, star_se)  + '\n')
        seq_map.write(seq_table([rnafold_dict[info_list[1]][0], 'RNAfold', '-'], mature_se, star_se)  + '\n')
        seq_map.write(seq_table([centroid_dict[info_list[1]][0], 'Centroidfold', '-'], mature_se, star_se)  + '\n')
        ##
        one_iso_info = {}
        seq_file = '{0}/{1}.out'.format(out_data, info_list[1])
        if os.path.isfile(seq_file):
            with open(seq_file) as iso_info:
                for iso_line in iso_info.readlines():
                    iso_line = iso_line.strip().split('\t')
                    one_iso_info[iso_line[2]] = iso_line[8:10]
        ##
        for ii in range(len(iso_df)):
            text_seq = tagId_seq[iso_df.iloc[ii, 0]]
            info_tpm = tagId_count[text_seq]
            if text_seq in one_iso_info:
                text_seq_info = one_iso_info[text_seq]
                if text_seq_info[0] != '-':
                    mutant_site = int(text_seq_info[0].split(':')[0])
                    text_seq = site_tolow(text_seq, mutant_site)
                info_type = text_seq_info[1]
            else:
                info_type = '-'
            if pre_loc[3] == '+':
                text_left = iso_df.iloc[ii, 1]-int(pre_loc[1])
                text_right = int(pre_loc[2])-iso_df.iloc[ii, 2]
            else:
                text_left = int(pre_loc[2])-iso_df.iloc[ii, 2]
                text_right = iso_df.iloc[ii, 1]-int(pre_loc[1])
            only_seq = "."*text_left + text_seq + "."*text_right
            seq_map.write(seq_table([only_seq, info_type,
                                     info_tpm], mature_se, star_se)  + '\n')
        seq_map.write('</table>')


def randomString(length):
    raw = ""
    range1 = range(58, 65) # between 0~9 and A~Z
    range2 = range(91, 97) # between A~Z and a~z
    i = 0
    while i < length:
        seed = random.randint(48, 122)
        if ((seed in range1) or (seed in range2)):
            continue;
        raw += chr(seed);
        i += 1
    return raw

from optparse import OptionParser
parser = OptionParser()
parser.add_option('-c','--outpath',action = 'store',type = "string" ,dest = 'outpath')
parser.add_option('-f','--fasta',action = 'store',type = "string" ,dest = 'fasta')
parser.add_option('-t','--txtfile',action = 'store',type = "string" ,dest = 'txtfile')
parser.add_option('-m','--mergefile',action = 'store',type = "string" ,dest = 'mergefile')
parser.add_option('-v','--version', action="store_false", dest="verbose", default='',help="version [default]")

# fakeArgs = ['-c','.', '-f', 'aradata/reads_18_26.fa', '-t', 'aradata/reads_18_26.txt', '-m', 'aradata/tmp_all_miRNAs.txt',
#             '-o', 'aradata/out_pool_merge.txt', '-r', 'aradata/corr_pre-mirna_seq.txt']
(options,args)=parser.parse_args()

## pdf and html report
out_path = options.outpath
os.makedirs(out_path, exist_ok=True)
out_data = out_path + '/data'
os.makedirs(out_data, exist_ok=True)

print(out_path)

global MIN_SPACE
MIN_SPACE = 5
tag_id_abu( options.fasta )

parse_dict = {}

with open( options.mergefile ) as file_input:
    for eli in file_input.readlines()[1:]:
        eli = eli.strip().split('\t')
        parse_dict[eli[1]] = eli

rnafold_dict = rnafold_fold(parse_dict)
centroid_dict = centroid_fold(parse_dict)
print(options.txtfile)
iso_dict = extract_map_read( options.txtfile, parse_dict)

print(len(parse_dict))
bias_dict = bias_get(parse_dict, iso_dict)

print(len(parse_dict), "======\n")

def create_ps_map(tmp_info):
    structure_plot(tmp_info)
    miR_iso_add(tmp_info)

print('Parent process {0} is Running'.format(os.getpid()))

p = multiprocessing.Pool(processes=2)
for poolsite in parse_dict:
    p.apply_async(create_ps_map, args=(parse_dict[poolsite],))
print('process start')
p.close()
p.join()

print('All processes done!')

# pool = multiprocessing.Pool(processes=2)
# pool.map(create_ps_map, list(parse_dict.keys()))

with open( '{}/{}'.format(out_path, "out_pool_merge.txt"), 'w') as out_file:
    inline_len = len(parse_dict[list(parse_dict.keys())[1]])
    if inline_len >12:
        tmp_name = ['Precursors', 'Extended_stem_loop_loc', 'Extended_stem_loop_seq', 'Extended_stem_loop_len',
                    'Stem_loop_loc', 'Stem_loop_seq', 'Loc5p', 'Seq5p', 'Len5p', 'Loc3p', 'Seq3p', 'Len3p', 'Mature_arm',
                    ['-']*(inline_len-12), 'Source', 'TPM5p', 'TPM3p', 'Stem_loop_len', 'Stem_loop_MFE',
                    'Stem_loop_AMFE', 'The_total_abundance', 'The_number_of_sequences_in_miRNA/miRNA*_and_3nt_variant_region',
                    'The_number_of_sequences_in_pre-miRNAs', 'Abundance_bias', 'Strand_bias',
                    'RNAfold', 'Centroidfold']
        tmp_tname = []
        for i in tmp_name:
            if type(i) == list:
                for j in i:
                    tmp_tname.append(j)
            else:
                tmp_tname.append(i)
        out_file.write('\t'.join(tmp_tname) + '\n')
    else:
        out_file.write('\t'.join(['Precursors', 'Extended_stem_loop_loc', 'Extended_stem_loop_seq',
                                  'Extended_stem_loop_len', 'Stem_loop_loc', 'Stem_loop_seq',
                                  'Loc5p', 'Seq5p', 'Len5p', 'Loc3p', 'Seq3p', 'Len3p', 'Mature_arm',
                                  'Source', 'TPM5p', 'TPM3p', 'Stem_loop_len', 'Stem_loop_MFE',
                                  'Stem_loop_AMFE', 'The_total_abundance', 'The_number_of_sequences_in_miRNA/miRNA*_and_3nt_variant_region',
                                  'The_number_of_sequences_in_pre-miRNAs', 'Abundance_bias', 'Strand_bias',
                                  'RNAfold', 'Centroidfold']) + '\n')
    ## write to text
    for pre_site in parse_dict:
        basic_info = parse_dict[pre_site]
        with open('{0}/{1}_fold.txt'.format(out_data, basic_info[1]), 'w') as fold_summary:
            ## expression value
            if basic_info[5] in tagId_count:
                rmpl = tagId_count[basic_info[5]]
            else:
                rmpl = 0
            if basic_info[8] in tagId_count:
                rmpf = tagId_count[basic_info[8]]
            else:
                rmpf = 0
            ## stem_loop
            loc_f = basic_info[4].split(':')
            loc_t = basic_info[7].split(':')
            if basic_info[5] in basic_info[2] and basic_info[8] in basic_info[2]:
                pos_f = getMiRNAPosition(basic_info[5], basic_info[2])
                pos_t = getMiRNAPosition(basic_info[8], basic_info[2])
                loc_list = [loc_f[1], loc_f[2], loc_t[1], loc_t[2]]
                loc_list = [int(i) for i in loc_list]
                stem_loop_loc = ':'.join([loc_f[0],str(min(loc_list)),str(max(loc_list)),loc_f[3]])
                if pos_f[0] > pos_t[1]:
                    pass
                    # print(basic_info[1])
                else:
                    stem_loop_seq = basic_info[2][pos_f[0]:pos_t[1]]
                    stem_loop_seq = stem_loop_seq.replace('T', 'U')
                    ## MFE and AMFE
                    # print(basic_info, stem_loop_seq, "\n")
                    other, stem_loop_structure = subprocess.getstatusoutput("echo " + stem_loop_seq + "| RNAfold --noPS")
                    rnafold_tmp = re.split("\n|\s\(", stem_loop_structure)
                    rmfe = re.sub(r'\)|\s', '', rnafold_tmp[2])
                    ramfe = float(rmfe)/len(stem_loop_seq)*100
                    ramfe = round(ramfe, 2)
                    ## structure
                    rnafold_res_o = struct_filter([basic_info[1], basic_info[2], rnafold_dict[pre_site][0]], basic_info[5])
                    rnafold_res_t = struct_filter([basic_info[1], basic_info[2], rnafold_dict[pre_site][0]], basic_info[8])
                    if rnafold_res_o[0] == 'False' and rnafold_res_t[0] == 'False':
                        rnafold_out = 'False'
                    else:
                        rnafold_out = 'True'
                    fold_summary.write('{}\t{}\n'.format(basic_info[5], '\t'.join(rnafold_res_o)))
                    centroid_res_o = struct_filter([basic_info[1], basic_info[2], centroid_dict[pre_site][0]], basic_info[5])
                    centroid_res_t = struct_filter([basic_info[1], basic_info[2], centroid_dict[pre_site][0]], basic_info[8])
                    if centroid_res_o[0] == 'False' and centroid_res_t[0] == 'False':
                        centroid_out = 'False'
                    else:
                        centroid_out = 'True'
                    fold_summary.write('{}\t{}\n'.format(basic_info[5], '\t'.join(centroid_res_o)))
                    ## output list
                    basic_info[2] = basic_info[2].replace('T', 'U')
                    basic_info[5] = basic_info[5].replace('T', 'U')
                    basic_info[8] = basic_info[8].replace('T', 'U')
                    basic_info[3] = len(basic_info[2])
                    infor_list = [basic_info[0:4], stem_loop_loc, stem_loop_seq, basic_info[4:], rmpl, rmpf,
                                    len(stem_loop_seq), rmfe, ramfe, bias_dict[pre_site], rnafold_out, centroid_out]
                    infor_out = []
                    for i in infor_list:
                        if type(i) == list:
                            for j in i:
                                infor_out.append(str(j))
                        else:
                            infor_out.append(str(i))
                    ## write to table
                    out_file.write('\t'.join(infor_out) +'\n')
