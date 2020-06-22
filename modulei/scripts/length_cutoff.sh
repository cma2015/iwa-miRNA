
set -x
awk -v RS=">" -v aa=$2 -v bb=$3 'NR>1{if(length($2)>=aa&&length($2)<=bb){print ">"$1"\n"$2}}' $1 >$4