

awk 'NR>1{print ">"$1"#5p#"$9"\n"$8"\n>"$1"#3p#"$12"\n"$11}' out_pool_merge.txt >mature_seq.fa

ssearch36 -m 8 -E 1 mature_seq.fa mature_seq.fa >mature_seq.out

awk -F"\t" '{OFS="\t";split($1,c,"#");split($2,d,"#");a=c[3];b=int($3*$4/100+0.5-$6);if(a-b<=2&&c[1]!=d[1]) print c[1],d[1]}' mature_seq.out \
 | sort -k 1,1 -k 2,2 | uniq >Conserved.txt

cut -f 1 out_pool_merge.txt | awk -F"\t" 'NR==FNR{a[$1]=$1;a[$2]=$2;next}{if(!a[$1]) print $1}' Conserved.txt - >Non_conserved.txt

Rscript igraph_rename.R

awk -F"\t" 'NR==FNR{a[$1]=$2;next}{OFS="\t";{print $1,a[$1]}}' NameChange.txt out_pool_merge.txt  >NameChange_out.txt
