
## all 152 features
################################
$fn="seq.fasta";

## miPred features
## miPred sequential and structural features
system("miPred/RNAfold -noPS < $fn > ../data/$fn.fold");
system("perl miPred/genRNAStats.pl < ../data/$fn.fold > ../data/$fn.data1");
system("miPred/RNAspectral.exe < ../data/$fn.fold > ../data/$fn.data2");

## z-features
system("perl miPred/genRandomRNA.pl -n 10 -m d < $fn > ../data/$fn.random.fasta");
system("miPred/RNAfold < ../data/$fn.random.fasta > ../data/$fn.random.fold");
system("perl miPred/genRNARandomStats.pl -n 10 -i ../data/$fn.random.fold -o ../data/$fn.zdata -m ../data/$fn.fold");

## MFEI1, MFEI2, MFEI3, MFEI4
system("java mfe14 $fn");
system("java mfe23 $fn");

## basepair-related features
system("java bpcount1 $fn");
system("java bpcount2 $fn");

## RNAfold-related features
system("RNAfold/RNAfold -p2 < $fn > ../data/$fn.RNAfold1");  
system("java RNAfoldfilter $fn");
system("rm *.ps");

## Mfold-related features
system("perl melt.pl $fn > ../data/$fn.mfold");
system("java Mfoldfilter $fn"); 

## filtering all features
system("java filter_all $fn");

## get the parameters besides microPred
system("miPred/RNAfold < $fn > ../data/$fn.fold");
system("perl fold_combine.pl ../data/$fn.fold");
system("perl parameters.pl ../data/$fn.fold1");

##combine all the 152 features
system("perl combine.pl");
