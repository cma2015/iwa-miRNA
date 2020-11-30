import os
import urllib
import urllib.request
import re
import time
import sys
from Bio import SeqIO

miRNAinput = ''
with open('mature_miRNAs.fa', 'r') as sRNAfile:
    for srna in sRNAfile.readlines():
        srna = srna.strip()
        miRNAinput = miRNAinput + '\n' + srna

miRNAinput = miRNAinput.replace('\n', '', 1)
print( 'miRNAs have been uploaded!')

data = {}
with open('{}/Transcripts/transcripts.fa'.format(sys.argv[1])) as f:
    for seq in SeqIO.parse(f, "fasta"):
        data[">" + seq.id] = str(seq.seq)

targetinput = '\n'.join(['\n'.join(list(i)) for i in data.items()])

print(len(targetinput))

# targetinput= ''
# for tarseq in data.keys():
#     targetinput = targetinput + '\n' + ">" + tarseq + "\n" + data[tarseq]
# targetinput = targetinput.replace('\n', '', 1)

print ('Trnascripts have been uploaded!\n')

print( "Start connect the psRNATarget website!")
url = 'http://plantgrn.noble.org/psRNATarget/analysis?function=3'
values = {'function' : 3,
    'srna_content' : miRNAinput,
    'target_content' : targetinput,
    'curschema':'s2',
    'top':200,
    'expect':3,
    'gup':0.5,
    'misp':1,
    'seedfactor':1.5,
    'seedpos1':2,
    'seedpos2':13,
    'maxnummismatchinseed':2,
    'hspsize':19,
    'allowbulge':'yes',
    'gapstartp':2,
    'gapextp':2,
    'cutpos1':10,
    'cutpos2':11
    }

data = bytes(urllib.parse.urlencode(values), encoding='utf8')
req = urllib.request.Request(url, data)
response = urllib.request.urlopen(req)
the_page = response.read().decode('utf-8')

back = re.findall(r'result\?sessionid=\d+', the_page)
print(back)
downpath = 'http://plantgrn.noble.org/psRNATarget/result?sessionid='+back[0].split('=')[1]+'\&contentype=text'

with open(sys.argv[2], 'w') as pspath:
    pspath.write('wget ' + downpath + ' -O psRNAtarget_MIT.out' + '\n')
