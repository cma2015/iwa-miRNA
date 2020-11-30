#! /usr/bin/env python

import re
import sys
from collections import Counter

with open(sys.argv[1]) as tmp_seq:
    adapter = []
    base_seq = "TGACAGAAGAGAGTGAGCAC"
    seq_list = tmp_seq.readlines()
    if base_seq in seq_list:
        adatper_seq = '-'
    else:
        seq_len=len(seq_list)
        for eli in seq_list:
            eli = eli.strip()
            adaseq = re.search(base_seq, eli).span()[1]
            adapter.append(eli[int(adaseq):(int(adaseq)+12)])
        adatper_count = Counter(adapter).most_common(1)[0][1]
        if float(adatper_count)/seq_len>0.5:
            adatper_seq = Counter(adapter).most_common(1)[0][0]
        else:
            adatper_seq = 'NA'

print(adatper_seq)
