head -n -6 ../model06/bestrun/model06_maxL.par | cat - ~/data_analysis/costataeEvolution/Result/08fastsimcoal/par_ending

################# par_ending file
#//Number of independent loci [chromosome]
#14 0
#//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
#1
#//per Block:data type, number of loci, per gen recomb and mut rates
#DNA 555000 0 9.5e-9 OUTEXP

~/biosoft/fsc27_linux64/fsc27093 -i model06_maxL.bootstrap.par -j -d -s0 -x -I -q -c 4 -B 4 -n100 -k50000000