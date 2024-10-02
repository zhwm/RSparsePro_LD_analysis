library(susieR)
library(data.table)
library(PRROC)
args = commandArgs(trailingOnly=T)
zfile = args[1]
ztrue = args[2]
ldfile = args[3]

gwas = as.matrix(fread(zfile, header=T))
gtrue = as.matrix(fread(ztrue, header=T))
ld = as.matrix(fread(ldfile))

safe_susie = function(zscore, ld) {
  result = tryCatch(
    {
      model = susie_rss(zscore, ld)
      model
    },
    error = function(e) {
      return(NA)
    }
  )
  if (!is.list(result)) {
    exist = 0
    convrg = 0
    pip = rep('NA', length(zscore))
    csvec = rep('NA', length(zscore))
  } else {
    exist = 1
    convrg = ifelse(model$converged, 1, 0)
    pip = model$pip
    csvec = rep(0, length(zscore))
    if (length(model$sets$cs)>0) {
      for (i in 1:length(model$sets$cs)) {
        csvec[model$sets$cs[[i]]] = i
      }
    }
  }
  return(list(exist = exist, convrg = convrg, pip = pip, csvec = csvec))
}

out_matrix_z_pip = matrix(nrow = nrow(ld), ncol = ncol(gwas))
out_matrix_z_cs = matrix(nrow = nrow(ld), ncol = ncol(gwas))
colnames(out_matrix_z_pip) = colnames(out_matrix_z_cs) = colnames(gwas)

out_matrix_ztrue_pip = matrix(nrow = nrow(ld), ncol = ncol(gtrue))
out_matrix_ztrue_cs = matrix(nrow = nrow(ld), ncol = ncol(gtrue))
colnames(out_matrix_ztrue_pip) = colnames(out_matrix_ztrue_cs) = colnames(gtrue)

out_matrix_summary = matrix(ncol = 2, nrow = ncol(gtrue)+ncol(gwas))
rownames(out_matrix_summary) = c(colnames(gtrue), colnames(gwas))
colnames(out_matrix_summary) = c('exist', 'convrg')

for (ite in 0:49){
  res = safe_susie(gtrue[, paste0('ite', ite)], ld)
  out_matrix_summary[paste0('ite', ite), 'exist'] = res$exist
  out_matrix_summary[paste0('ite', ite), 'convrg'] = res$convrg
  out_matrix_ztrue_pip[, paste0('ite', ite)] = res$pip
  out_matrix_ztrue_cs[, paste0('ite', ite)] = res$csvec
  for (prop in c('0.001', '0.01', '0.1')){
    for (flip in c('-1.0', '-0.5', '0.0', '0.5')){
      res = safe_susie(gwas[, paste0('f', flip, '_p', prop, '_ite', ite)], ld)
      out_matrix_summary[paste0('f', flip, '_p', prop, '_ite', ite), 'exist'] = res$exist
      out_matrix_summary[paste0('f', flip, '_p', prop, '_ite', ite), 'convrg'] = res$convrg
      out_matrix_z_pip[, paste0('f', flip, '_p', prop, '_ite', ite)] = res$pip
      out_matrix_z_cs[, paste0('f', flip, '_p', prop, '_ite', ite)] = res$csvec
    }
  }
}
write.table(out_matrix_summary,paste0(zfile, ".susie.summary.txt"),col.names=T,row.names=F,sep="\t",quote=F)
write.table(out_matrix_ztrue_pip,paste0(zfile, ".susie.ztrue.pip.txt"),col.names=T,row.names=F,sep="\t",quote=F)
write.table(out_matrix_ztrue_cs,paste0(zfile, ".susie.ztrue.cs.txt"),col.names=T,row.names=F,sep="\t",quote=F)
write.table(out_matrix_z_pip,paste0(zfile, ".susie.z.pip.txt"),col.names=T,row.names=F,sep="\t",quote=F)
write.table(out_matrix_z_cs,paste0(zfile, ".susie.z.cs.txt"),col.names=T,row.names=F,sep="\t",quote=F)