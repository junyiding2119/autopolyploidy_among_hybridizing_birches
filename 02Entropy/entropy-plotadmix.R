## R script to plot admixture proportions from inds
## and output text files with Q-matrix coefficients


library(rhdf5)
library(coda)
set.seed(10041996)

for(k in 2:8){
    for (chain in 1:4){
        qch1<-h5read(paste0("mcmcoutchain", chain, ".k",k,".hdf5"),"q")

        qest1<-apply(qch1, c(2,3), mean)

        # creating a simple Q-matrix file for use with CLUMPAK, pong, etc
        write.table(format(t(qest1), digits=5), file=paste0("simpleQchain",chain,".k",k,".txt"),
                    row.names=F, col.names=F, quote=F)
	}
}

