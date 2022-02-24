library(sequenza)
args = commandArgs(T)
# args = c('small.seqz.gz','RA14-gDNA')

test <- sequenza.extract(args[1], verbose = T, chromosome.list=c(paste0("chr", c(1:22, "X", "Y"))))
CP <- sequenza.fit(test)
sequenza.results(sequenza.extract = test, cp.table = CP, sample.id = args[2], out.dir=paste0(args[2],"-sequenza"))
print("finished sequenza")
