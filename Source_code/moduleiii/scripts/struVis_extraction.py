import os
import sys

with open(sys.argv[4]) as miRNA_input:
    for eli in miRNA_input.readlines()[1:]:
        eli = eli.strip().split("\t")
        mir_location = eli[4].split(":")
        chr_star_end = "{}:{}-{}".format(mir_location[0], mir_location[1], mir_location[2])
        if mir_location[3]=="+":
            strand = "plus"
        else:
            strand = "minus"
        os.system("{}/strucVis -b {} -g {} -c {} -s {} -p {}.ps -n {}".format(
            sys.argv[1], sys.argv[2], sys.argv[3], chr_star_end, strand, eli[4], eli[1]
        ))
        os.system("gmt psconvert -Tf {0}.ps && pdf2svg {0}.pdf ../miRNASelection/data/{0}.svg".format(eli[4]))

os.system("rm *.ps *.pdf")