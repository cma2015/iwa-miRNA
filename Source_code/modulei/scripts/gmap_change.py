
import re

## get options
from optparse import OptionParser

parser = OptionParser()
# parser.add_option('-t','--type',action = 'store',type = "string" ,dest = 'type')
parser.add_option('-i','--inputfile',action = 'store',type = "string" ,dest = 'inputfile')
parser.add_option('-p','--psl',action = 'store',type = "string" ,dest = 'psl')
parser.add_option('-o','--outfile',action = 'store',type = "string" ,dest = 'outfile')
parser.add_option('-v','--version', action="store_false", dest="verbose", default='',help="version [default]")
(options,args)=parser.parse_args()

out_dict = {}
tmp_dict = {}
loc_name = {}
with open( options.psl ) as input_psl:
    for eli in input_psl.readlines():
        eli = eli.strip().split('\t')
        loc_line = ':'.join([eli[13], eli[15], eli[16], eli[8]])
        loc_name[loc_line] = eli[9]
        if eli[9] not in out_dict and loc_line not in tmp_dict:
            out_dict[eli[9]] = loc_line
            tmp_dict[loc_line] = 'op'
        else:
            pass

print(len(out_dict))
## multiple names
remain_loc = [x for x in loc_name.keys() if x not in tmp_dict.keys()]
remain_loc_add = {}
for i in remain_loc:
    if loc_name[i] not in remain_loc_add:
        remain_loc_add[loc_name[i]] = 1
        tmp_name = loc_name[i] + '_' + str(remain_loc_add[loc_name[i]])
        out_dict[tmp_name] = i
    else:
        remain_loc_add[loc_name[i]] += 1
        tmp_name = loc_name[i] + '_' + str(remain_loc_add[loc_name[i]])
        out_dict[tmp_name] = i

print(len(out_dict))

def loc_change(pre_loc, mat_loc, new_pre):
    pre_loc = pre_loc.split(':')
    mat_loc = mat_loc.split(':')
    new_pre = new_pre.split(':')
    print(pre_loc, mat_loc, new_pre)
    if pre_loc[3] == '+':
        loc1 = int(mat_loc[1]) - int(pre_loc[1])
        loc2 = int(mat_loc[2]) - int(pre_loc[1])
    else:
        loc1 = int(pre_loc[2]) - int(mat_loc[2])
        loc2 = int(pre_loc[2]) - int(mat_loc[1])
    if new_pre[3] == '+':
        new1 = int(new_pre[1]) + loc1
        new2 = int(new_pre[1]) + loc2
    else:
        new1 = int(new_pre[2]) - loc2
        new2 = int(new_pre[2]) - loc1
    return(':'.join([new_pre[0], str(new1), str(new2), new_pre[3]]))

## location file
tmp_dict = {}
with open( options.inputfile ) as input_txt:
    for eli in input_txt.readlines()[1:]:
        eli = eli.strip().split('\t')
        tmp_dict[eli[0]] = eli

## Conversion location
with open( options.outfile, 'w') as result_file:
    for i in out_dict:
        if i in tmp_dict:
            eli = tmp_dict[i]
            print(eli)
            mat_res = loc_change(eli[1], eli[5], out_dict[eli[0]])
            star_res = loc_change(eli[1], eli[8], out_dict[eli[0]])
            print(eli)
            result_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                i, out_dict[i], eli[2], eli[3], eli[4], mat_res, eli[6],
                eli[7], star_res, eli[9], eli[10]))
        elif re.sub("_\d$", "", i) in tmp_dict:
            eli = tmp_dict[re.sub("_\d$", "", i)]
            mat_res = loc_change(eli[1], eli[5], out_dict[i])
            star_res = loc_change(eli[1], eli[8], out_dict[i])
            print(eli)
            result_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                i, out_dict[i], eli[2], eli[3], eli[4], mat_res, eli[6],
                eli[7], star_res, eli[9], eli[10]))
