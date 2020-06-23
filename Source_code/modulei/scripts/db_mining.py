import yaml
import subprocess
import os
import re
import random
import urllib.request
from Bio import SeqIO

## function
def isParenthese(ch):
    if( ch=='(' or ch==')' ):
        return True
    else:
        return False

def getMiRNALen (miRNASeq):
    return len(miRNASeq)

def getMiRNAPosition (miRNASeq, RNASequence):
    startPos = RNASequence.find(miRNASeq)
    endPos = startPos + getMiRNALen(miRNASeq)
    return [startPos, endPos]

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

def getMiRNAStar(RNASequence, RNAStructure, miRNASeq, matchedStructList = [] ):
    if( len(matchedStructList) == 0 ):
        matchedStructList = matchRNAStructure( RNAStructure )
    [miRNAPosStart, miRNAPosEnd] = getMiRNAPosition( miRNASeq, RNASequence )
    [matchedPosStart, matchedPosEnd] = getMatchedPositions( miRNAPosStart, miRNAPosEnd-2, matchedStructList )
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

def seqTostr( RNASequence):
    import re
    import subprocess
    other, structure = subprocess.getstatusoutput("echo " + RNASequence + "| RNAfold --noPS")
    structure = re.split("\s", structure)[1]
    return structure

def getLocation( RNALoc, RNASequence, miRNASeq):
    [startPos, endPos] = getMiRNAPosition(miRNASeq, RNASequence)
    [chromo, start, end, strand] = RNALoc.split(':')
    if strand == '+':
        miRLoc = ':'.join([chromo, str(int(start)+startPos),
                         str(int(start)+endPos-1), strand])
    else:
        miRLoc = ':'.join([chromo, str(int(end)-endPos+1),
                         str(int(end)-startPos), strand])
    return miRLoc

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

## get options
from optparse import OptionParser
parser = OptionParser()
# parser.add_option('-t','--type',action = 'store',type = "string" ,dest = 'type')
parser.add_option('-i','--curpath',action = 'store',type = "string" ,dest = 'curpath')
parser.add_option('-d','--database',action = 'store',type = "string" ,dest = 'database')
parser.add_option('-s','--species',action = 'store',type = "string" ,dest = 'species')
parser.add_option('-o','--outfile',action = 'store',type = "string" ,dest = 'outfile')
parser.add_option('-f','--outformat',action = 'store',type = "string" ,dest = 'outformat')
parser.add_option('-v','--version', action="store_false", dest="verbose", default='',help="version [default]")

# fakeArgs = ['-s','Zea_mays', '-d', 'miRBase,PmiREN,sRNAanno,PlantsmallRNAgenes',
#             '-i', '~/sRNAbox/tools/srnaDatabase', '-f', 'Rmarkdown']

(options,args)=parser.parse_args()

query_type = 'miRNA'
input_species = options.species
input_database = options.database
cur_dir = options.curpath
out_dir = '{}/tmp/{}'.format(cur_dir, 'vvvvv') #randomString(10)
os.makedirs(out_dir, exist_ok=True)
os.makedirs('{}/source'.format(cur_dir), exist_ok=True)
out_file = options.outfile
outformat = options.outformat
input_database = input_database.split(',')

## input data
with open("{}/table_data/data_info.yml".format(cur_dir)) as conf_file:
    conf_data = yaml.load(conf_file, Loader=yaml.FullLoader)
species_list = [x for k in conf_data for x in conf_data[k]]

species_database = [k for k in conf_data if input_species in conf_data[k]
                    and k in input_database]
print(species_database)

## database file
for tmp_name in species_database:
    os.makedirs('{}/source/{}'.format(cur_dir, tmp_name), exist_ok=True)

## import data source
url_dict = {}
with open("{}/table_data/datasource.txt".format(cur_dir), 'r') as ds:
    for eli in ds.readlines()[1:]:
        eli = eli.strip().split('\t')
        if eli[0] not in url_dict:
            url_dict[eli[0]] = {}
            url_dict[eli[0]][eli[1]] = [eli[2], eli[3]]
        else:
            url_dict[eli[0]][eli[1]] = [eli[2], eli[3]]

## download miRBase data
if 'miRBase' in species_database:
    print("download miRBase data")
    abbr_name = conf_data['miRBase'][input_species]
    tmp_file = '{}/source/miRBase/{}.gff3'.format(cur_dir, abbr_name)
    hc_file = '{}/source/miRBase/hairpin_high_conf.fa'.format(cur_dir)
    hairpin_file = '{}/source/miRBase/hairpin.fa'.format(cur_dir)
    mature_file = '{}/source/miRBase/mature.fa'.format(cur_dir)
    if not os.path.exists(hc_file):
        url_mirbase = 'ftp://mirbase.org/pub/mirbase/CURRENT/hairpin_high_conf.fa.gz'.format(abbr_name)
        os.system('wget -O {}.gz {} && gunzip {}.gz'.format(hc_file, url_mirbase, hc_file))
    if not os.path.exists(hairpin_file):
        url_mirbase = 'ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz'.format(abbr_name)
        os.system('wget -O {}.gz {} && gunzip {}.gz'.format(hairpin_file, url_mirbase, hairpin_file))
    if not os.path.exists(mature_file):
        url_mirbase = 'ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz'.format(abbr_name)
        os.system('wget -O {}.gz {} && gunzip {}.gz'.format(mature_file, url_mirbase, mature_file))
    if not os.path.exists(tmp_file):
        url_mirbase = 'ftp://mirbase.org/pub/mirbase/CURRENT/genomes/{}.gff3'.format(abbr_name)
        try:
            # print('wget -O {} {}'.format(tmp_file, url_mirbase))
            os.system('wget -O {} {}'.format(tmp_file, url_mirbase))
        except:
            pass

## download PmiREN data
if 'PmiREN' in species_database:
    print("Extract data from PmiREN")
    abbr_pmir = '{}{}'.format(input_species.split('_')[0][0:1], input_species.split('_')[1][0:2])
    pmir_path = '{}/source/PmiREN/{}'.format(cur_dir, input_species)
    if not os.path.exists(pmir_path):
        url_mirbase = 'http://www.pmiren.com/ftp-download/{}_{}/'.format(input_species, abbr_pmir)
        os.system('wget -r -np -nH {} && mv ftp-download/{}_{} {}/source/PmiREN/{}'.format(
            url_mirbase, input_species, abbr_pmir, cur_dir, input_species))

## download sRNAanno data
if 'sRNAanno' in species_database:
    print("download sRNAanno data")
    if query_type in ['miRNA', 'PHAS21', 'PHAS24', 'hc-siRNA']:
        tmp_file = '{}/source/sRNAanno/{}_{}.gff3'.format(cur_dir, input_species, query_type)
        if not os.path.exists(tmp_file):
            url_sRNAanno = 'http://www.plantsrnas.org/{}/{}.{}.gff3'.format(query_type, input_species, query_type)
            # print('wget -O {} {}'.format(tmp_file, url_sRNAanno))
            os.system('wget -O {} {}'.format(tmp_file, url_sRNAanno))

## download PlantsmallRNAgenes data
if 'PlantsmallRNAgenes' in species_database:
    print("download PlantsmallRNAgenes data")
    out_gff = '{}/source/PlantsmallRNAgenes/{}.gff3'.format(cur_dir, input_species)
    out_sRNA_fa = '{}/source/PlantsmallRNAgenes/{}_sRNA.fasta'.format(cur_dir, input_species)
    out_bed = '{}/source/PlantsmallRNAgenes/{}.bed'.format(cur_dir, input_species)
    out_genome = '{}/source/PlantsmallRNAgenes/{}_genome.fasta'.format(cur_dir, input_species)
    out_sequence = '{}/source/PlantsmallRNAgenes/{}_sequence.fasta'.format(cur_dir, input_species)
    if not os.path.exists(out_sequence):
        abbr_number = conf_data['PlantsmallRNAgenes'][input_species]
        gff_url = 'https://plantsmallrnagenes.science.psu.edu/annotations_dumper.php?genomes_id={}&dump_type=gff3'.format(abbr_number)
        # print(abbr_number, out_gff, gff_url)
        urllib.request.urlretrieve(gff_url, out_gff)
        fa_url = 'https://plantsmallrnagenes.science.psu.edu/annotations_dumper.php?genomes_id={}&dump_type=fasta_sRNA'.format(abbr_number)
        urllib.request.urlretrieve(fa_url, out_sRNA_fa)
        with open(out_bed, 'w') as output_bed:
            with open(out_gff) as input_gff:
                for eli in input_gff.readlines():
                    if not eli.startswith('#'):
                        eli = eli.strip().split('\t')
                        tmp_feat={}
                        for x in eli[8].split(';'):
                            tmp_feat.update({x.split('=')[0]:x.split('=')[1]})
                        output_bed.write('{}\t{}\t{}\t{}\t.\t{}\n'.format(eli[0], int(eli[3])-1, eli[4],
                                                                         tmp_feat['Name'], eli[6]))
        ## extract sequence
        genome_url = url_dict['PlantsmallRNAgenes'][input_species][0]
        urllib.request.urlretrieve(genome_url, out_genome)
        # os.system('wget {} {}'.format(out_genome, genome_url))
        os.system('bedtools getfasta -s -name -fi {} -bed {} | sed -r "s/\(.\)//" > {}'.format(
            out_genome, out_bed, out_sequence))

## Extract data from miRBase:
if 'miRBase' in species_database:
    print("Extract data from miRBase")
    abbr_name = conf_data['miRBase'][input_species]
    high_conf = []
    tmp_dir = '{}/source/miRBase'.format(cur_dir)
    for rec in SeqIO.parse('{}/hairpin_high_conf.fa'.format(tmp_dir), "fasta"):
        high_conf.append(rec.id)
    pre_dict = {}
    for rec in SeqIO.parse('{}/hairpin.fa'.format(tmp_dir), "fasta"):
        pre_dict[rec.id] = str(rec.seq)
    mature_dict = {}
    for rec in SeqIO.parse('{}/mature.fa'.format(tmp_dir), "fasta"):
        mature_dict[rec.id] = str(rec.seq)
    out_mirbase = {}
    count_pre = {}
    with open('{}/{}.gff3'.format(tmp_dir, abbr_name), 'r') as data_mirbase:
        for eli in data_mirbase.readlines():
            if not eli.startswith('#'):
                eli = eli.strip().split('\t')
                tmp_feat={}
                for x in eli[8].split(';'):
                    tmp_feat.update({x.split('=')[0]:x.split('=')[1]})
                eli[0] = re.sub(r'[Cc]hr', '', eli[0])
                eli[0] = re.sub(r'\D+-', '', eli[0])
                paste_loc = '{}:{}:{}:{}'.format(eli[0], eli[3], eli[4], eli[6])
                if tmp_feat['Name'] in pre_dict:
                    seq_out = pre_dict[tmp_feat['Name']]
                elif re.sub('\.\d', '', tmp_feat['Name']) in pre_dict:
                    seq_out = pre_dict[re.sub('\.\d', '', tmp_feat['Name'])]
                else:
                    seq_out = mature_dict[tmp_feat['Name']]
                table_info = [tmp_feat['Name'], paste_loc, seq_out, len(seq_out)]
                if 'Derives_from' not in tmp_feat: # and tmp_feat['Alias'] not in out_mirbase:
                    count_pre[tmp_feat['Alias']] = tmp_feat['ID']
                    if table_info[0] not in count_pre:
                        count_pre[table_info[0]] = 0
                    else:
                        count_pre[table_info[0]] += 1
                        table_info[0] = table_info[0]  + "_" + str(count_pre[table_info[0]])
                    out_mirbase[tmp_feat['ID']] = 12*['-']
                    out_mirbase[tmp_feat['ID']][0:4] = table_info
                    if tmp_feat['Name'] in high_conf:
                        out_mirbase[tmp_feat['Alias']][11] = 'High-conf'
                    else:
                        out_mirbase[tmp_feat['Alias']][11] = 'Non-high-conf'
                elif 'Derives_from' in tmp_feat and tmp_feat['Derives_from'] in out_mirbase:
                    pre_seq = out_mirbase[count_pre[tmp_feat['Derives_from']]][2]
                    pre_str = seqTostr(pre_seq)
                    arm_out = checkArm(pre_seq, pre_str, seq_out)
                    mat_vec = tmp_feat['Name'].split('-')
                    if arm_out == '5p':
                        out_mirbase[count_pre[tmp_feat['Derives_from']]][5:8] = [paste_loc, seq_out, len(seq_out)]
                    elif arm_out == '3p':
                        out_mirbase[count_pre[tmp_feat['Derives_from']]][8:11] = [paste_loc, seq_out, len(seq_out)]
                    else:
                        pass
    ## Complete the corresponding sequence
    for eli in out_mirbase:
        infor_list = out_mirbase[eli]
        pre_seq = infor_list[2]
        pre_str = seqTostr(infor_list[2])
        if infor_list[5] == '-' and infor_list[8] != '-':
            mirseq = infor_list[9]
            # print(pre_seq, pre_str, mirseq)
            seq_out = getMiRNAStar(pre_seq, pre_str, mirseq)
            out_mirbase[eli][5] = getLocation(infor_list[1], infor_list[2], seq_out)
            out_mirbase[eli][6] = seq_out
            out_mirbase[eli][7] = len(seq_out)
        elif infor_list[8] == '-' and infor_list[5] != '-':
            mirseq = infor_list[6]
            seq_out = getMiRNAStar(pre_seq, pre_str, mirseq)
            out_mirbase[eli][8] = getLocation(infor_list[1], infor_list[2], seq_out)
            out_mirbase[eli][9] = seq_out
            out_mirbase[eli][10] = len(seq_out)
    print('pre_miRNA_number:{}, pre_in_gff:{}'.format(str(len(pre_dict)), str(len(out_mirbase))))
    ## Output results
    with open('{}/miRBase.txt'.format(out_dir), 'w') as out_tmp:
        out_tmp.write('Pre-miRNAs\tpLocation\tpSequence\tpLength\tMature_arm\t\
Location5p\tSequence5p\tLength5p\tLocation3p\tSequence3p\tLength3p\tConf\n')
        for i in sorted(list(out_mirbase.keys())):
            out_list = [str(x) for x in out_mirbase[i]]
            out_tmp.write('{}\n'.format('\t'.join(out_list)))

## Extract data from PmiREN:
if 'PmiREN' in species_database:
    print("Extract data from PmiREN")
    pmir_path = '{}/source/PmiREN/{}'.format(cur_dir, input_species)
    out_pmir = {}
    with open('{}/{}_basicInfo.txt'.format(pmir_path, input_species)) as data_pmir:
        for eli in data_pmir.readlines()[1:]:
            if len(eli)>0:
                eli = eli.split('\t')
                eli[4] = re.sub(r'[Cc]hr', '', eli[4])
                eli[4] = re.sub(r'\D+-', '', eli[4])                
                pre_loc = '{}:{}:{}:{}'.format(eli[4], eli[5], eli[6], eli[7])
                mat_loc = '{}:{}:{}:{}'.format(eli[4], eli[15], eli[16], eli[7])
                star_loc = '{}:{}:{}:{}'.format(eli[4], eli[19], eli[20], eli[7])
                mat_pos = checkArm(eli[8], eli[9], eli[17])
                if mat_pos == '5p':
                    out_pmir[eli[0]] = [eli[0], pre_loc, eli[8], len(eli[8]), mat_pos,
                                       mat_loc, eli[17], len(eli[17]),
                                       star_loc, eli[21], len(eli[21]),]
                else:
                    out_pmir[eli[0]] = [eli[0], pre_loc, eli[8], len(eli[8]), mat_pos,
                                       star_loc, eli[21], len(eli[21]),
                                       mat_loc, eli[17], len(eli[17]),]
    ## Output results
    with open('{}/PmiREN.txt'.format(out_dir), 'w') as out_tmp:
        out_tmp.write('Pre-miRNAs\tpLocation\tpSequence\tpLength\tMature_ter\t\
Location5p\tSequence5p\tLength5p\tLocation3p\tSequence3p\tLength3p\n')
        for i in sorted(list(out_pmir.keys())):
            out_list = [str(x) for x in out_pmir[i]]
            out_tmp.write('{}\n'.format('\t'.join(out_list)))
    ## Copy file
    import shutil
    shutil.copy('{}/{}_mature.tissue.exp.txt'.format(pmir_path, input_species),
                '{}/PmiREN_tissue.txt'.format(out_dir))


## Extract data from sRNAanno:
if 'sRNAanno' in species_database:
    print("Extract data from sRNAanno")
    out_srna = {}
    with open('{}/source/sRNAanno/{}_{}.gff3'.format(cur_dir, input_species, query_type)) as data_srna:
        for eli in data_srna.readlines():
            if not eli.startswith('#') and len(eli)>8:
                eli = eli.strip().split('\t')
                tmp_feat={}
                for x in eli[8].split(';'):
                    tmp_feat.update({x.split('=')[0]:x.split('=')[1]})
                eli[0] = re.sub(r'[Cc]hr', '', eli[0])
                eli[0] = re.sub(r'\D+-', '', eli[0])                    
                paste_loc = '{}:{}:{}:{}'.format(eli[0], eli[3], eli[4], eli[6])
                table_info = [tmp_feat['ID'], paste_loc, tmp_feat['seq'], len(tmp_feat['seq'])]
                if eli[1] not in out_srna:
                    out_srna[eli[1]] = 11*['']
                    if 'Derives_from' not in tmp_feat:
                        out_srna[eli[1]][0:4] = table_info
                elif 'Derives_from' in tmp_feat and tmp_feat['Derives_from'] in out_srna:
                    mat_vec = tmp_feat['ID'].split('-')
                    if mat_vec[-2] == '5p':
                        out_srna[eli[1]][5:8] = [paste_loc, tmp_feat['seq'], len(tmp_feat['seq'])]
                    else:
                        out_srna[eli[1]][8:11] = [paste_loc, tmp_feat['seq'], len(tmp_feat['seq'])]
                    if mat_vec[-1] == 'mature':
                        out_srna[eli[1]][4] = mat_vec[-2]
    ## Output results
    with open('{}/sRNAanno.txt'.format(out_dir), 'w') as out_tmp:
        out_tmp.write('Pre-miRNAs\tpLocation\tpSequence\tpLength\tMature_ter\t\
Location5p\tSequence5p\tLength5p\tLocation3p\tSequence3p\tLength3p\n')
        for i in sorted(list(out_srna.keys())):
            out_list = [str(x) for x in out_srna[i]]
            out_tmp.write('{}\n'.format('\t'.join(out_list)))

## Extract data from PlantsmallRNAgenes:
if query_type == 'miRNA':
    psrg = ['MIRNA'] #, 'nearMIRNA'
elif query_type == 'siRNA':
    psrg = ['siRNA21', 'siRNA22', 'siRNA23', 'siRNA24']

if 'PlantsmallRNAgenes' in species_database:
    print("Extract data from PlantsmallRNAgenes")
    pre_dict = {}
    tmp_dir = '{}/source/PlantsmallRNAgenes'.format(cur_dir)
    for rec in SeqIO.parse('{}/{}_sequence.fasta'.format(tmp_dir, input_species), "fasta"):
        pre_dict[rec.id] = re.sub('T', 'U', str(rec.seq))
    mature_dict = {}
    for rec in SeqIO.parse('{}/{}_sRNA.fasta'.format(tmp_dir, input_species), "fasta"):
        tmp_read={}
        for x in rec.description.split(' ')[1:]:
            tmp_read.update({x.split('=')[0]:x.split('=')[1]})
        mature_dict[rec.id] = [str(rec.seq), tmp_read['reads']]
    out_anno = {}
    with open('{}/{}.gff3'.format(tmp_dir, input_species), 'r') as data_anno:
        for eli in data_anno.readlines()[1:]:
            eli = eli.split('\t')
            if eli[2] in psrg:
                tmp_feat={}
                for x in eli[8].split(';'):
                    tmp_feat.update({x.split('=')[0]:x.split('=')[1]})
                eli[0] = re.sub(r'[Cc]hr', '', eli[0])
                eli[0] = re.sub(r'\D+-', '', eli[0])                    
                pre_loc = '{}:{}:{}:{}'.format(eli[0], eli[3], eli[4], eli[6])
                pre_seq = pre_dict[tmp_feat['Name']]
                pre_str = seqTostr(pre_seq)
                mat_seq = mature_dict[tmp_feat['Name']+'_majorRNA'][0]
                mat_abundance = mature_dict[tmp_feat['Name']+'_majorRNA'][1]
                if mat_seq in pre_seq:
                    star_seq = getMiRNAStar(pre_seq, pre_str, mat_seq)
                    mat_loc = getLocation(pre_loc, pre_seq, mat_seq)
                    star_loc = getLocation(pre_loc, pre_seq, star_seq)
                    arm_out = checkArm(pre_seq, pre_str, mat_seq)
                    # print(tmp_feat['Name'], pre_loc, pre_seq, len(pre_seq), arm_out, mat_abundance)
                    out_anno[tmp_feat['Name']] = [tmp_feat['Name'], pre_loc, pre_seq, len(pre_seq), arm_out, mat_abundance]
                    if arm_out == '5p':
                        out_anno[tmp_feat['Name']].extend([mat_loc, mat_seq, len(mat_seq),
                                                        star_loc, star_seq, len(star_seq), eli[1], eli[2]])
                    else:
                        out_anno[tmp_feat['Name']].extend([star_loc, star_seq, len(star_seq),
                                                        mat_loc, mat_seq, len(mat_seq), eli[1], eli[2]])
    ## Output results
    with open('{}/PlantsmallRNAgenes.txt'.format(out_dir), 'w') as out_tmp:
        out_tmp.write('Pre-miRNAs\tpLocation\tpSequence\tpLength\tMature_ter\treadnumber\t\
Location5p\tSequence5p\tLength5p\tLocation3p\tSequence3p\tLength3p\tSource\tType\n')
        for i in sorted(list(out_anno.keys())):
            out_list = [str(x) for x in out_anno[i]]
            out_tmp.write('{}\n'.format('\t'.join(out_list)))

## Rmarkdown file
def report_out():
    with open("{}/output.Rmd".format(out_dir), 'w') as out_report:
        out_report.write(side_bar)
        ## miRBase
        if 'miRBase' in species_database:
            out_report.write(miRBase_report)
        ## PmiREN
        if 'PmiREN' in species_database:
            out_report.write(PmiREN_report)
        ## sRNAanno
        if 'sRNAanno' in species_database:
            out_report.write(sRNAanno_report)
        ## Plant small RNA
        if 'PlantsmallRNAgenes' in species_database:
            out_report.write(PlantsmallRNAgenes_report)
        ### output
        # print('Rscript -e "rmarkdown::render(\'{0}/output.Rmd\',output_dir=\'{0}\')"'.format(out_dir))
        os.system(
        'Rscript -e "rmarkdown::render(\'{0}/output.Rmd\',output_dir=\'{0}\')"'.format(out_dir))
        os.system('cp {0}/output.html {1}'.format(out_dir, out_file))

side_bar = """---
title: "The available miRNAs in several public databases"
output:
    flexdashboard::flex_dashboard:
        orientation: rows
        vertical_layout: scroll
        social: menu
        source_code: embed
        theme: cosmo
---

```{r setup, include=FALSE}
library(DT)
library(purrr)
library(ggplot2)
library(dplyr)
library(highcharter)
options(stringsAsFactors = F)
knitr::opts_chunk$set(echo = FALSE,
                      message = FALSE,
                      warning = FALSE)
```

```{css, echo=FALSE}
.limitrow .chart-stage {
width:1500px;
overflow-x:scroll;
}

p.seq {
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    width: 200px;
    min-width: 0;
    margin: 0;
    font-size: 15px;
}

p.seq:hover {
    overflow: none;
    width: auto;
}
```

Sidebar {.sidebar}
=====================================

There are four public databases containing miRNAs in plants.

1. [miRBase](http://www.mirbase.org)

2. [PmiREN](http://pmiren.com/)

3. [sRNAanno](http://www.plantsrnas.org/)

4. [Plant small RNA gene](https://plantsmallrnagenes.science.psu.edu/index.php)

In this section, we summarize an overview of miRNAs in related databases according to the input species.
"""

miRBase_report = """
miRBase {data-icon=fa-bar-chart}
=====================================

### Introduction

The miRBase database is a searchable database of published miRNA sequences and annotation. Each entry in the miRBase Sequence database represents a predicted hairpin portion of a miRNA transcript (termed mir in the database), with information on the location and sequence of the mature miRNA sequence (termed miR). Both hairpin and mature sequences are available for searching and browsing, and entries can also be retrieved by name, keyword, references and annotation. All sequence and annotation data are also available for download.

### Result

Row
-------------------------------------

### Table: Detailed information of pre-miRNAs and miRNAs {.limitrow}

```{r }
sRNA_data <- read.table("miRBase.txt", sep = "\t", header = T)
sRNA_data <- sRNA_data[sRNA_data[,8]<30&sRNA_data[,11]<30, ]
sRNA_data[,3] <- paste0("<p class=\\"seq\\">", sRNA_data[,3],"</p>")
datatable(data = sRNA_data, 
          extensions = 'Buttons',
          options = list(dom = "Blfrtip",
                         buttons = list("copy",list(extend = "collection",
                                                     buttons = c("csv", "excel"),
                                                     text = "Download")),
                         lengthMenu = list( c(10, 20, -1), c(10, 20, "All")),
                         pageLength = 10, autoWidth = TRUE),
          rownames = FALSE, escape = FALSE, class = 'compact nowrap stripe hover')%>%
formatStyle(TRUE, `text-align` = 'center')
```

Row
-------------------------------------

### Categories of miRNA identification

```{r}
name_col_two <- as.data.frame(table(sRNA_data[,12]))
name_col_two[,1] <- paste0(name_col_two[,1], " (", name_col_two[,2], ")")
highchart() %>%
    hc_title(text = "miRNAs" ,align = "center",verticalAlign = "middle") %>%
    hc_tooltip(headerFormat ="", pointFormat = "{series.y} <b>{point.percentage:.1f}%</b>")%>%
    hc_plotOptions(pie = list(dataLabels = list(enabled = TRUE,distance = -50
                                                ,style = list(fontWeight = "bold",color = "white",
                                                              fontSize=16, textOutline = "")),
                              center = c('50%','50%'))) %>%
    hc_add_series(name_col_two, type = "pie", hcaes(name = Var1, y = Freq),
            innerSize = "50%") %>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Pie_miRNA_identification")

```

### miRNA family sizes and total number of miRNA family

```{r}
name_family <- gsub("\\\\D$", "", sRNA_data[,1])
family_table <- table(name_family)
name_col_pre <- table(family_table)
name_col_num <- as.numeric(names(name_col_pre))
tmp_name <- as.character(min(name_col_num):max(name_col_num))
name_col_pre <- name_col_pre[tmp_name]
names(name_col_pre) <- tmp_name
name_col_pre[is.na(name_col_pre)] <- 0
name_col_one <- as.data.frame(name_col_pre)
name_col_one[,1] <- as.numeric(name_col_one[,1])
highchart() %>%
    hc_title(text = "miRNA family") %>%
    hc_add_series(name="miRNA family", name_col_one, type = "column",hcaes(name = Var1, y = Freq))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 12)))%>%
    hc_exporting(enabled = TRUE, filename = "total number of miRNA family")
```

### miRNA family sizes

```{r}
family_index <- name_family%in%names(family_table[family_table>1])
family_show <- data.frame('mirna' = name_family[family_index])
family_show %>%
    count(mirna) %>%
    hchart('treemap', hcaes(x = 'mirna', value = 'n', color = 'n'))%>%
    hc_exporting(enabled = TRUE, filename = "miRNA family sizes")
```

Row
-------------------------------------

### Length distribution

```{r}
ds <- map(unique(sRNA_data[,12]), function(x){
    dt <- density(sRNA_data[,4][sRNA_data[,12] == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>%
    hc_add_series_list(ds)%>%
    hc_add_theme(hc_theme_google())
```

### Length and distribution of all miRNAs

```{r}
ter_len <- data.frame('Ter'=rep(c("5p","3p"), each = nrow(sRNA_data)),
                      "Len" = c(sRNA_data[,8], sRNA_data[,11]))
ter_len <- ter_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(ter_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution")
```

### Composition of the first base

```{r}
mat_char <- data.frame('Type'=rep(c("5p","3p"), each = nrow(sRNA_data)),
                      "Letter" = c(substr(sRNA_data[,7], 1, 1), substr(sRNA_data[,10], 1, 1)))
mat_char <- mat_char %>% group_by(Letter, Type) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal", borderColor = "",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_char, type = "column", hcaes(x = Type, y = n , group = Letter))%>%
    hc_add_theme(hc_theme_google()) %>%
    hc_tooltip(headerFormat ="", pointFormat = "<b>{point.percentage:.1f}%</b>")%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 16)))%>%
    hc_exporting(enabled = TRUE, filename = "Composition of the first base")
```
"""

PmiREN_report = """
PmiREN {data-icon=fa-bar-chart}
=====================================

### Introduction

PmiREN (Plant miRNA ENcyclopedia) is a comprehensive functional plant miRNA database. In current version, PmiREN contains 20388 miRNA loci (MIRs) belonging to 5757 families, 1365 clusters, 1668 syntenic blocks and 141327 predicted miRNA-target pairs in 88 species phylogenetically ranging from chlorophytes to angiosperms. In addition, 1537 deeply sequenced small RNA libraries were used in quantification of miRNA expression patterns, and 116 PARE-Seq libraries were employed to validate predicted miRNA-target pairs. PmiREN has many unique features such as annotating miRNAs with an uniform method, solid and extensive information on miRNA synteny, conservation, expression, references, and providing many convenient data accesses, powerful search and download engines. We believe PmiREN is a truly functional knowledgebase and will provide novel insights into plant miRNA research and benefit the whole research community.

### Result

Row
-------------------------------------

### Table: Detailed information of pre-miRNAs and miRNAs {.limitrow}

```{r }
sRNA_data2 <- read.table("PmiREN.txt", sep = "\t", header = T)
sRNA_data2[,3] <- paste0("<p class=\\"seq\\">", sRNA_data2[,3],"</p>")
datatable(data = sRNA_data2, 
          extensions = 'Buttons',
          options = list(dom = "Blfrtip",
                         buttons = list("copy",list(extend = "collection",
                                                     buttons = c("csv", "excel"),
                                                     text = "Download")),
                         lengthMenu = list( c(10, 20, -1), c(10, 20, "All")),
                         pageLength = 10, autoWidth = TRUE),
          rownames = FALSE, escape = FALSE, class = 'compact nowrap stripe hover')%>%
formatStyle(TRUE, `text-align` = 'center')
```

Row
-------------------------------------

### miRNA family sizes and total number of miRNA family

```{r}
name_family <- gsub("\\\\D$", "", sRNA_data2[,1])
family_table <- table(name_family)
name_col_pre <- table(family_table)
name_col_num <- as.numeric(names(name_col_pre))
tmp_name <- as.character(min(name_col_num):max(name_col_num))
name_col_pre <- name_col_pre[tmp_name]
names(name_col_pre) <- tmp_name
name_col_pre[is.na(name_col_pre)] <- 0
name_col_one <- as.data.frame(name_col_pre)
name_col_one[,1] <- as.numeric(name_col_one[,1])
highchart() %>%
    hc_title(text = "miRNA family") %>%
    hc_add_series(name="miRNA family", name_col_one, type = "column",hcaes(name = Var1, y = Freq))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 12)))%>%
    hc_exporting(enabled = TRUE, filename = "total number of miRNA family")
```

### miRNA family sizes

```{r}
family_index <- name_family%in%names(family_table[family_table>1])
family_show <- data.frame('mirna' = name_family[family_index])
family_show %>%
    count(mirna) %>%
    hchart('treemap', hcaes(x = 'mirna', value = 'n', color = 'n'))%>%
    hc_exporting(enabled = TRUE, filename = "miRNA family sizes")
```

### Categories of miRNA identification

```{r}
def_col_two <- rep("Known", nrow(sRNA_data2))
def_col_two[grepl("MIRN", sRNA_data2[,1])] <- "Novel"
name_col_two <- as.data.frame(table(def_col_two))
name_col_two[,1] <- paste0(name_col_two[,1], " (", name_col_two[,2], ")")
highchart() %>%
    hc_title(text = "miRNAs" ,align = "center",verticalAlign = "middle") %>%
    hc_tooltip(headerFormat ="", pointFormat = "{series.y} <b>{point.percentage:.1f}%</b>")%>%
    hc_plotOptions(pie = list(dataLabels = list(enabled = TRUE,distance = -50
                                                ,style = list(fontWeight = "bold",color = "white",
                                                              fontSize=16, textOutline = "")),
                              center = c('50%','50%'))) %>%
    hc_add_series(name_col_two, type = "pie", hcaes(name = def_col_two, y = Freq),
            innerSize = "50%") %>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Pie_miRNA_identification")
```

Row
-------------------------------------

### Length distribution

```{r}
ds <- map(unique(def_col_two), function(x){
    dt <- density(sRNA_data2[,4][def_col_two == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>%
    hc_add_series_list(ds)%>%
    hc_add_theme(hc_theme_google())
```

### Length and distribution of all miRNAs

```{r}
ter_len <- data.frame('Ter'=rep(c("5p","3p"), each = nrow(sRNA_data2)),
                      "Len" = c(sRNA_data2[,8], sRNA_data2[,11]))
ter_len <- ter_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(ter_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution")
```

### Length and distribution of mature miRNAs

```{r}
mat_len <- data.frame('Ter'=sRNA_data2[,5],
                      "Len" = sRNA_data2[,8])
mat_len <- mat_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Length and distribution of mature miRNA")
```
"""

sRNAanno_report = """
sRNAanno {data-icon=fa-bar-chart}
=====================================

### Introduction

Small RNAs (sRNAs) are essential regulatory molecules in plants. The sRNAanno database hosts a large collection of miRNA, phasiRNA- and hc-siRNA-generating loci annotated from ~140 plants using consistent and high-confident criteria. All the annotations are made freely available to the scientific community via various services and tools.

### Result

Row
-------------------------------------

### Table: Detailed information of pre-miRNAs and miRNAs {.limitrow}

```{r }
sRNA_data3 <- read.table("sRNAanno.txt", sep = "\t", header = T)
sRNA_data3[,3] <- paste0("<p class=\\"seq\\">", sRNA_data3[,3],"</p>")
datatable(data = sRNA_data3, 
          extensions = 'Buttons',
          options = list(dom = "Blfrtip",
                         buttons = list("copy",list(extend = "collection",
                                                     buttons = c("csv", "excel"),
                                                     text = "Download")),
                         lengthMenu = list( c(10, 20, -1), c(10, 20, "All")),
                         pageLength = 10, autoWidth = TRUE),
          rownames = FALSE, escape = FALSE, class = 'compact nowrap stripe hover')%>%
formatStyle(TRUE, `text-align` = 'center')
```

Row
-------------------------------------

### miRNA family sizes and total number of miRNA family

```{r}
name_col <- t(do.call("cbind", strsplit(sRNA_data3[,1], "-")))
name_col[,1] <- gsub("\\\\D$", "", name_col[,1] )
colnames(name_col) <- c("mirna", 'Type')

name_col_pre <- table(table(name_col[,1]))
name_col_num <- as.numeric(names(name_col_pre))
tmp_name <- as.character(min(name_col_num):max(name_col_num))
name_col_pre <- name_col_pre[tmp_name]
names(name_col_pre) <- tmp_name
name_col_pre[is.na(name_col_pre)] <- 0
name_col_one <- as.data.frame(name_col_pre)
name_col_one[,1] <- as.numeric(name_col_one[,1])

highchart() %>%
    hc_title(text = "miRNA family") %>%
    hc_add_series(name="miRNA family", name_col_one, type = "column",hcaes(name = Var1, y = Freq))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 12)))%>%
    hc_exporting(enabled = TRUE, filename = "total number of miRNA family")
```

### miRNA family sizes

```{r}
family_index <- name_col[,1]%in%names(table(name_col[,1])[table(name_col[,1])>1])
family_show <- as.data.frame(name_col[family_index, ])
family_show %>%
    count(mirna) %>%
    hchart('treemap', hcaes(x = 'mirna', value = 'n', color = 'n'))%>%
    hc_exporting(enabled = TRUE, filename = "miRNA family sizes")
```

### Categories of miRNA identification

```{r}
name_col_two <- as.data.frame(table(name_col[,2]))
name_col_two[,1] <- paste0(name_col_two[,1], " (", name_col_two[,2], ")")
highchart() %>%
    hc_title(text = "miRNAs" ,align = "center",verticalAlign = "middle") %>%
    hc_tooltip(headerFormat ="", pointFormat = "{series.y} <b>{point.percentage:.1f}%</b>")%>%
    hc_plotOptions(pie = list(dataLabels = list(enabled = TRUE,distance = -50
                                                ,style = list(fontWeight = "bold",color = "white",
                                                              fontSize=16, textOutline = "")),
                              center = c('50%','50%'))) %>%
    hc_add_series(name_col_two, type = "pie", hcaes(name = Var1, y = Freq),
            innerSize = "50%") %>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Pie_miRNA_identification")
```

Row
-------------------------------------

### Length distribution

```{r}
ds <- map(unique(name_col[,2]), function(x){
    dt <- density(sRNA_data3[,4][name_col[,2] == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>%
    hc_add_series_list(ds)%>%
    hc_add_theme(hc_theme_google())

ter_len <- data.frame('Ter'=rep(c("5p","3p"), each = nrow(sRNA_data3)),
                      "Len" = c(sRNA_data3[,8], sRNA_data3[,11]))
ter_len <- ter_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(ter_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution")
```

### Length and distribution of mature miRNA

```{r}
mat_len <- data.frame('Ter'=sRNA_data3[,5],
                      "Len" = sRNA_data3[,8])
mat_len <- mat_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Length and distribution of mature miRNA")
```

### Composition of the first base

```{r}
mat_char <- data.frame('Type'=name_col[,2],
                      "Letter" = substr(sRNA_data3[,7], 1, 1))
mat_char <- mat_char %>% group_by(Letter, Type) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal", borderColor = "",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_char, type = "column", hcaes(x = Type, y = n , group = Letter))%>%
    hc_add_theme(hc_theme_google()) %>%
    hc_tooltip(headerFormat ="", pointFormat = "<b>{point.percentage:.1f}%</b>")%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 16)))%>%
    hc_exporting(enabled = TRUE, filename = "Composition of the first base")
```
"""

PlantsmallRNAgenes_report = """
Plant small RNA genes {data-icon=fa-bar-chart}
=====================================

### Introduction

Small non-coding RNAs between 20 and 24 nucleotides in length function to regulate gene expression, as well as suppress accumulation of viruses. In plants, there are several distinct types of endogenous small RNAs. MicroRNAs (miRNAs) are produced from single-stranded stem-loop RNA precursors, and function to target mRNAs and long non-coding RNAs for post-transcriptional regulation. In plants most miRNAs are 21 or 22 nucleotides in length. Endogenous short interfering RNAs (siRNAs) arise from double-stranded RNA precursors that are frequently produced by endogenous RNA-dependent RNA polymerases. Some siRNAs are 21 or 22 nucleotides long and can function like miRNAs to post-transcriptionally regulate mRNAs. Many more siRNAs are 24 nts in length, and function to target chromatin-associated non-coding RNAs and reinforce repressive chromatin modifcations. MicroRNAs have been well-annotated in many species, but siRNA locus annotations have lagged.

### Result

Row
-------------------------------------

### Table: Detailed information of pre-miRNAs and miRNAs {.limitrow}

```{r }
sRNA_data4 <- read.table("PlantsmallRNAgenes.txt", sep = "\t", header = T)
data_index <- sRNA_data4$Mature_ter!='unmatchedRegion'&sRNA_data4$Length5p<30&sRNA_data4$Length3p<30
sRNA_data4 <- sRNA_data4[data_index, ]
sRNA_data4[,3] <- paste0("<p class=\\"seq\\">", sRNA_data4[,3],"</p>")
datatable(data = sRNA_data4, 
          extensions = 'Buttons',
          options = list(dom = "Blfrtip",
                         buttons = list("copy",list(extend = "collection",
                                                     buttons = c("csv", "excel"),
                                                     text = "Download")),
                         lengthMenu = list( c(10, 20, -1), c(10, 20, "All")),
                         pageLength = 10, autoWidth = TRUE),
          rownames = FALSE, escape = FALSE, class = 'compact nowrap stripe hover')%>%
formatStyle(TRUE, `text-align` = 'center')
```

Row
-------------------------------------

### Categories of miRNA identification

```{r}
name_col_two <- as.data.frame(table(sRNA_data4$Source))
name_col_two[,1] <- paste0(name_col_two[,1], " (", name_col_two[,2], ")")
highchart() %>%
    hc_title(text = "miRNAs" ,align = "center",verticalAlign = "middle") %>%
    hc_tooltip(headerFormat ="", pointFormat = "{series.y} <b>{point.percentage:.1f}%</b>")%>%
    hc_plotOptions(pie = list(dataLabels = list(enabled = TRUE,distance = -50
                                                ,style = list(fontWeight = "bold",color = "white",
                                                              fontSize=16, textOutline = "")),
                              center = c('50%','50%'))) %>%
    hc_add_series(name_col_two, type = "pie", hcaes(name = Var1, y = Freq),
            innerSize = "50%") %>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Pie_miRNA_identification")
```

### Categories of miRNA identification

```{r}
name_col_two <- as.data.frame(table(sRNA_data4$Type))
name_col_two[,1] <- paste0(name_col_two[,1], " (", name_col_two[,2], ")")
highchart() %>%
    hc_title(text = "miRNAs" ,align = "center",verticalAlign = "middle") %>%
    hc_tooltip(headerFormat ="", pointFormat = "{series.y} <b>{point.percentage:.1f}%</b>")%>%
    hc_plotOptions(pie = list(dataLabels = list(enabled = TRUE,distance = -50
                                                ,style = list(fontWeight = "bold",color = "white",
                                                              fontSize=16, textOutline = "")),
                              center = c('50%','50%'))) %>%
    hc_add_series(name_col_two, type = "pie", hcaes(name = Var1, y = Freq),
            innerSize = "50%") %>%
    hc_add_theme(hc_theme_google())%>%
    hc_exporting(enabled = TRUE, filename = "Pie_miRNA_identification")
```


### miRNA family sizes and total number of miRNA family

```{r}
name_family <- gsub("\\\\D$", "", sRNA_data4$Pre.miRNAs)
family_table <- table(name_family)
name_col_pre <- table(family_table)
name_col_num <- as.numeric(names(name_col_pre))
tmp_name <- as.character(min(name_col_num):max(name_col_num))
name_col_pre <- name_col_pre[tmp_name]
names(name_col_pre) <- tmp_name
name_col_pre[is.na(name_col_pre)] <- 0
name_col_one <- as.data.frame(name_col_pre)
name_col_one[,1] <- as.numeric(name_col_one[,1])
highchart() %>%
    hc_title(text = "miRNA family") %>%
    hc_add_series(name="miRNA family", name_col_one, type = "column",hcaes(name = Var1, y = Freq))%>%
    hc_add_theme(hc_theme_google())%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 12)))%>%
    hc_exporting(enabled = TRUE, filename = "total number of miRNA family")
```

### miRNA family sizes

```{r}
family_index <- name_family%in%names(family_table[family_table>1])
family_show <- data.frame('mirna' = name_family[family_index])
family_show %>%
    count(mirna) %>%
    hchart('treemap', hcaes(x = 'mirna', value = 'n', color = 'n'))%>%
    hc_exporting(enabled = TRUE, filename = "miRNA family sizes")
```

Row
-------------------------------------

### Length distribution

```{r}
ds <- map(unique(sRNA_data4$Source), function(x){
    dt <- density(sRNA_data4$pLength[sRNA_data4$Source == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

db <- map(unique(sRNA_data4$Type), function(x){
    dt <- density(sRNA_data4$pLength[sRNA_data4$Type == x])[1:2]
    dt <- list_parse2(as.data.frame(dt))
    list(data = dt, name = x)
})

highchart() %>%
    hc_add_series_list(ds)%>%
    hc_add_series_list(db)%>%
    hc_add_theme(hc_theme_google())
```

### Length and distribution of all miRNAs

```{r}
ter_len <- data.frame('Ter'=rep(c("5p","3p"), each = nrow(sRNA_data4)),
                      "Len" = c(sRNA_data4$Length5p, sRNA_data4$Length3p))
ter_len <- ter_len %>% group_by(Len, Ter) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(ter_len, type = "column", hcaes(x = Len, y = n , group = Ter))%>%
    hc_exporting(enabled = TRUE, filename = "Length distribution")
```

### Composition of the first base

```{r}
firstbase <- apply(sRNA_data4, 1, function(x){
  if(x[5] == "5p"){
  return(substr(x[8], 1, 1))
  }else{
  return(substr(x[11], 1, 1))
  }
})
mat_char <- data.frame('Type'=rep(c(sRNA_data4$Source, sRNA_data4$Type)),
                      "Letter" = c(firstbase, firstbase))
mat_char <- mat_char %>% group_by(Letter, Type) %>%  summarise(n = n())
highchart() %>%
    hc_plotOptions(column = list(
        dataLabels = list(enabled = FALSE),
        stacking = "normal", borderColor = "",
        enableMouseTracking = TRUE)
    ) %>%
    hc_add_series(mat_char, type = "column", hcaes(x = Type, y = n , group = Letter))%>%
    hc_add_theme(hc_theme_google()) %>%
    hc_tooltip(headerFormat ="", pointFormat = "<b>{point.percentage:.1f}%</b>")%>%
    hc_xAxis(type = "category", labels = list(style = list(fontSize = 16)))%>%
    hc_exporting(enabled = TRUE, filename = "Composition of the first base")
```
"""

if outformat == "Rmarkdown":
    report_out();
else:
    print(out_dir)
