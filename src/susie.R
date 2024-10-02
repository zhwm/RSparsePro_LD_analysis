library(susieR)
library(data.table)

args = commandArgs(trailingOnly=T)
zfile = args[1]
ldfile = args[2]

gwas = fread(zfile)
ld = as.matrix(fread(ldfile))
ld[is.na(ld)] = 0
model = susie_rss(gwas$Z, ld)

gwas$PIP = round(model$pip, 4)
gwas$converged = model$converged
gwas$cs = 0
if (length(model$sets$cs)>0) {
    for (i in 1:length(model$sets$cs)) {
        gwas$cs[model$sets$cs[[i]]] = i
    }
}
write.table(gwas, paste0(zfile,".susie.txt"), col.names=T, row.names=F, sep="\t", quote=F)