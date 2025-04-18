gen=test
phen=test.pheno 
covar=test.cov
annot=single.annot
env=test.env

# time ../build/GENIE -g $gen -p $phen -c $covar -e $env -m G+GxE+NxE  -k 10 -jn 10  -o test.out --annot $annot -t 6 -s 1
# ../build/GENIE_multi_pheno -g $gen -p ./test.2.pheno -c $covar -e $env -m G+GxE+NxE  -k 10 -jn 10  -o test.2.out --annot $annot -t 6 -tr 1
# ../build/GENIE_mem --config ./config.txt

../build/GENIE_trace1 -g $gen -p ./test.pheno -c $covar -e $env -m G+GxE+NxE  -k 10 -jn 10  -o test.temp.out --annot $annot -t 6 -tr 1
../build/GENIE_trace2 -g $gen -p ./test.2.pheno -c $covar -e $env -m G+GxE+NxE  -k 10 -jn 10  -o test.2.trace.out --annot $annot -t 6 -tr_input test.temp.out.trace
