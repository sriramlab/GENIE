gen=test
phen=test.pheno 
covar=test.cov
annot=single.annot
env=test.env


../build/GENIE -g $gen -p $phen -c $covar -e $env -m G+GxE+NxE  -k 10 -jn 10  -o test.out.txt -annot $annot -t 6







