library("scarHRD")
args = commandArgs(T)
results = scar_score(args[1], "grch38", seqz=TRUE)
write.table(results, file=paste(args[2], "scarHRD.result.csv", sep="/"), sep=",", row.names=F, quote=T)

