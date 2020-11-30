
def base_change(basei):
  if basei == "A":
    return("U")
  elif basei == "C":
    return("G")
  elif basei == "G":
    return("C")
  elif basei == "T":
    return("A")
  elif basei == "U":
    return("A")

def UT_cahnge(basei):
  if basei == "T":
    return("U")
  else:
    return(basei)

inter_precusors = []
with open("SNP.vcf", "w") as output:
  output.write("##fileformat=VCFv4.2\n")
  output.write("#CHROM\tPOS\tID\tREF\tALT\n")
  with open("intersect_miSNPs.txt") as interfile:
    for eli in interfile.readlines():
      eli = eli.strip().split("\t")
      inter_precusors.append(eli[3])
      if eli[12] in ["A", "T", "C", "G", "U"] and eli[13] in ["A", "T", "C", "G", "U"]:
        if eli[5]=="+":
          start = str(int(eli[8])-int(eli[1]))
          output.write('{}\t{}\t{}\t{}\t{}\n'.format(eli[3], start, eli[9], 
          UT_cahnge(eli[12]), UT_cahnge(eli[13])))
        else:
          start = str(int(eli[2])-int(eli[7]))
          output.write('{}\t{}\t{}\t{}\t{}\n'.format(eli[3], start, eli[9], 
          base_change(eli[12]), base_change(eli[13])))

with open("inter_pre_miRNA.bed", "w") as interbed:
  with open("pre_miRNA.bed") as prebed:
    for eli in prebed.readlines():
      eli = eli.strip().split("\t") 
      if eli[3] in inter_precusors:
        interbed.write('\t'.join(eli) + "\n")

