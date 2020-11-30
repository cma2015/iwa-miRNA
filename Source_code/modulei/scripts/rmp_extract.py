#! /bin/python
import glob,os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-i','--path',action = 'store',type = "string" ,dest = 'path')
parser.add_option('-m','--mean',action = 'store',type = "string" ,dest = 'mean')
parser.add_option('-t','--total',action = 'store',type = "string" ,dest = 'total')
parser.add_option('-l','--least',action = 'store',type = "string" ,dest = 'least')
parser.add_option('-s','--sample',action = 'store',type = "string" ,dest = 'sample')
parser.add_option('-o','--output',action = 'store',type = "string" ,dest = 'output')
parser.add_option('-g','--log',action = 'store',type = "string" ,dest = 'log')
parser.add_option('-v','--version', action="store_false", dest="verbose", default='',help="version [default]")
(options,args)=parser.parse_args()

#meanval = options.mean
leastval = options.least
samval = options.sample
path = options.path
file = glob.glob(os.path.join(path, "*.txt"))

seqindex = {}
for f in file:
    with open(f, 'r') as rmpfile:
        for eli in rmpfile.readlines():
            eli = eli.strip().split('\t')
            eli[2] = float(eli[2])
            if eli[1] not in seqindex:
                seqindex[eli[1]] = [eli[2], eli[2], 1]
            else:
                seqindex[eli[1]][0] += eli[2]
                seqindex[eli[1]][2] += 1
                if eli[2] > seqindex[eli[1]][1]:
                    seqindex[eli[1]][1] = eli[2]

## write to file
#print('sum value:{}\tleast value{}'.format(len(file)*float(meanval), float(leastval)))

with open( options.log,'w') as logfile:
    with open( options.output, 'w') as outfile:
        tmp = 1
        for ii in seqindex:
            seqindex[ii] = [round(seqindex[ii][0],2), 
            round(seqindex[ii][1],2), seqindex[ii][2]]
            logfile.write('{}\t{}\t{}\n'.format(seqindex[ii][0], seqindex[ii][1], seqindex[ii][2]))
            if seqindex[ii][2] >= int(samval):
                if  seqindex[ii][1] >= float(leastval) or seqindex[ii][0] >= int(options.total):
                    outfile.write('>Seq{} {}\n{}\n'.format(tmp,seqindex[ii][0],ii))
                    tmp += 1
