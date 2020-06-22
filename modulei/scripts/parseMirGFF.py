
from optparse import OptionParser
parser = OptionParser()
parser.add_option('-i','--input',action = 'store',type = "string" ,dest = 'input')
parser.add_option('-o','--output',action = 'store',type = "string" ,dest = 'output')
parser.add_option('-v','--version', action="store_false", dest="verbose", default='',help="version [default]")
(options,args)=parser.parse_args()

out_anno = {}
with open( options.input ) as mirna_input:
    for eli in mirna_input.readlines():
        if not eli.startswith("#"):
            eli = eli.strip().split("\t")
            tmp_feat={}
            for x in eli[8].split(';'):
                tmp_feat.update({x.split('=')[0]:x.split('=')[1]})  
            if eli[2] == "pre-miRNA":
                out_anno[tmp_feat['ID']] = ['-']*11
                out_anno[tmp_feat['ID']][0:4] = [tmp_feat['ID'],tmp_feat['ID'],
                         tmp_feat['Seq'],len(tmp_feat['Seq'])]
            else:
                if tmp_feat['Parent'] in out_anno:
                    if eli[2] == 'ref_miRNA':
                        out_anno[tmp_feat['Parent']][10] = tmp_feat['Arm']
                    if tmp_feat['Arm'] == '5p':
                        out_anno[tmp_feat['Parent']][4:7] = [tmp_feat['ID'],
                                 tmp_feat['Seq'],len(tmp_feat['Seq'])]
                    elif tmp_feat['Arm'] == '3p':
                        out_anno[tmp_feat['Parent']][7:10] = [tmp_feat['ID'],
                                 tmp_feat['Seq'],len(tmp_feat['Seq'])]
                    else:
                        del out_anno[tmp_feat['Parent']]
                                        
with open(options.output, 'w') as out_tmp:
    out_tmp.write('Pre-miRNAs\tpLocation\tpSequence\tpLength\t\
Location5p\tSequence5p\tLength5p\tLocation3p\tSequence3p\tLength3p\tMature_arm\n')
    for i in sorted(list(out_anno.keys())):
        out_list = [str(x) for x in out_anno[i]]
        out_tmp.write('{}\n'.format('\t'.join(out_list)))