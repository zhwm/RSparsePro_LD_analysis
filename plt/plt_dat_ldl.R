library(data.table)
library(ggplot2)

## Example
### APOB
apob_ld = as.matrix(fread("../dat/LDL/bed/AMR_2_20767461_21767461.ld"))
apob = read.table("../dat/LDL/formatted/AMR_2_20767461_21767461.txt.rsparsepro.txt",header = T)

pltdat = data.frame(pos = apob$POS,
                    pval = -log10(exp(1)) * pchisq(apob$Z^2, df = 1, lower.tail = F, log.p = T),
                    cs = ifelse(apob$cs == 0, "Other variants",
                                paste0("Credible set ",apob$cs)),
                    r2 = as.numeric(apob_ld[which(apob$cs==1)[which.max(apob$PIP[apob$cs==1])],]^2))

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = pval, color = cs)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = "") +
  scale_color_manual(values = c("darkred","grey")) +
  scale_x_continuous(breaks = c(20.6,21.0,21.4,21.8), limits = c(20.6,21.8)) +
  xlab(paste0("Chr2 (Mb)")) -> plt
ggsave("../fig/dat_ldl_apob_cs.pdf",plt,height=3,width=10)

ggplot(pltdat, aes(x = pos/1000000, y = pval, color = r2)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = expression(r^2)) +
  scale_color_stepsn(n.breaks = 6, colours = c("black","blue","green","orange","red")) +
  scale_x_continuous(breaks = c(20.6,21.0,21.4,21.8), limits = c(20.6,21.8)) +
  xlab(paste0("Chr2 (Mb)")) -> plt
ggsave("../fig/dat_ldl_apob_cs1.pdf",plt,height=3,width=10)

### GCKR
ld = as.matrix(fread("../dat/LDL/bed/EUR_2_27230940_28230940.ld"))
stats = read.table("../dat/LDL/formatted/EUR_2_27230940_28230940.txt.rsparsepro.txt",header = T)
pltdat = data.frame(pos = stats$POS,
                  pval = -log10(exp(1)) * pchisq(stats$Z^2, df = 1, lower.tail = F, log.p = T),
                  cs = ifelse(stats$cs == 0, "Other variants",
                              paste0("Credible set ",stats$cs)),
                  r2 = as.numeric(ld[which(stats$cs==1)[which.max(stats$PIP[stats$cs==1])],]^2))
ggplot(pltdat, aes(x = pos/1000000, y = pval, color = cs)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = "") +
  scale_color_manual(values = c("darkred","grey")) +
  scale_x_continuous(limits = c(27.7,27.8)) +
  xlab(paste0("Chr2 (Mb)")) -> plt
ggsave("../fig/dat_ldl_gckr_cs.pdf",plt,height = 3,width = 10)

ggplot(pltdat, aes(x = pos/1000000, y = pval, color = r2)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = expression(r^2)) +
  scale_color_stepsn(n.breaks = 6, colours = c("black","blue","green","orange","red")) +
  scale_x_continuous(limits = c(27.7,27.8)) +
  xlab(paste0("Chr2 (Mb)")) -> plt
ggsave("../fig/dat_ldl_gckr_cs1.pdf",plt,height = 3,width = 10)

### LDLRAP1

ld = as.matrix(fread("../dat/LDL/bed/EUR_1_25268937_26268937.ld"))
stats = read.table("../dat/LDL/formatted/EUR_1_25268937_26268937.txt.rsparsepro.txt",header = T)

pltdat = data.frame(pos = stats$POS,
                  pval = -log10(exp(1)) * pchisq(stats$Z^2, df = 1, lower.tail = F, log.p = T),
                  cs = ifelse(stats$cs == 0, "Other variants",
                            paste0("Credible set ",stats$cs)),
                  r2 = as.numeric(ld[which(stats$cs==1)[which.max(stats$PIP[stats$cs==1])],]^2),
                  r3 = as.numeric(ld[which(stats$cs==2)[which.max(stats$PIP[stats$cs==2])],]^2))

ggplot(pltdat, aes(x = pos/1000000, y = pval, color = cs)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = "") +
  scale_color_manual(values = c("darkred","darkgreen","grey")) +
  scale_x_continuous(limits = c(25.6,26.0)) +
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_ldl_rhce_cs.pdf",plt,height = 3,width = 10)

ggplot(pltdat, aes(x = pos/1000000, y = pval, color = r2)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = expression(r^2)) +
  scale_color_stepsn(n.breaks = 6, colours = c("black","blue","green","orange","red")) +
  scale_x_continuous(limits = c(25.6,26.0)) +
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_ldl_rhce_cs1.pdf",plt,height = 3,width = 10)

ggplot(pltdat, aes(x = pos/1000000, y = pval, color = r3)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = expression(r^2)) +
  scale_color_stepsn(n.breaks = 6, colours = c("black","blue","green","orange","red")) +
  scale_x_continuous(limits = c(25.6,26.0)) +
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_ldl_rhce_cs2.pdf",plt,height = 3,width = 10)

### NR1H4
ld = as.matrix(fread("../dat/LDL/bed/EAS_12_100334253_101334253.ld"))
stats = read.table("../dat/LDL/formatted/EAS_12_100334253_101334253.txt.rsparsepro.txt",header = T)

pltdat = data.frame(pos = stats$POS,
                  pval = -log10(exp(1)) * pchisq(stats$Z^2, df = 1, lower.tail = F, log.p = T),
                  cs = ifelse(stats$cs == 0, "Other variants",
                              paste0("Credible set ",stats$cs)),
                  r2 = as.numeric(ld[which(stats$cs==1)[which.max(stats$PIP[stats$cs==1])],]^2))

ggplot(pltdat, aes(x = pos/1000000, y = pval, color = cs)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = "") +
  scale_color_manual(values = c("darkred","grey")) +
  scale_x_continuous(limits = c(100.7,101)) +
  xlab(paste0("Chr12 (Mb)")) -> plt
ggsave("../fig/dat_ldl_nr1h4_cs.pdf",plt,height = 3,width = 10)

ggplot(pltdat, aes(x = pos/1000000, y = pval, color = r2)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = expression(r^2)) +
  scale_color_stepsn(n.breaks = 6, colours = c("black","blue","green","orange","red")) +
  scale_x_continuous(limits = c(100.7,101)) +
  xlab(paste0("Chr12 (Mb)")) -> plt
ggsave("../fig/dat_ldl_nr1h4_cs1.pdf",plt,height = 3,width = 10)

## Summary
stats = as.data.frame(fread("../doc/dat_LDL_ncs_convrg.txt"))

data = data.frame(
    category = c("SuSiE inference converged (431 loci)","SuSiE inference did not converge (5 loci)","SuSiE prior estimation failed (55 loci)"),
    count = c(sum(stats$sus_convrg), sum(stats$sus_convrg==0)-sum(stats$sus_exist==0), sum(stats$sus_exist==0))
)

options(repr.plot.width=10, repr.plot.height=6)
ggplot(data, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = c("darkblue","steelblue1","grey")) +
  labs(fill = "")

pltdat = data.frame(num = c(stats$rsp_ncs,
                        stats$sus_ncs),
                method = rep(c("RSparsePro","SuSiE"),each = nrow(stats)),
                group = ifelse(stats$sus_exist==0,"SuSiE prior estimation failed",
                                ifelse(stats$sus_convrg==1,"SuSiE inference converged","SuSiE inference did not converge")),
                ancestry = stats$ancestry)
pltdat$label = paste0(pltdat$method, " results\nin loci where\n", pltdat$group)
pltdat$label = factor(pltdat$label, levels = unique(pltdat$label)[c(4,5,6,1,3,2)])

options(repr.plot.width=10, repr.plot.height=6)
ggplot(pltdat[pltdat$label!="SuSiE results\nin loci where\nSuSiE prior estimation failed",], aes(x = label, y = num, color = label)) +
    geom_violin(fill = "white") +
    geom_boxplot(width = 0.1, fill = "white") +
    theme_classic() +
    theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue","steelblue1","red","gold1","lightsalmon")) +
    labs(x = "", y = "Number of identified credible sets") +
    scale_y_continuous(breaks = seq(0,10,by = 2)) -> plt
ggsave("../fig/dat_ldl_csnum.pdf",plt,height = 6,width = 10)

anno = as.data.frame(fread("../doc/dat_LDL_anno_enrichment.txt"))

pltdat = data.frame(prop = c(anno$rsp_prop, anno$sus_prop), anno = anno$anno, group = anno$group, method = rep(c("RSparsePro","SuSiE"), each = nrow(anno)))
pltdat$label = paste0(pltdat$method, " results\nin loci where\n", pltdat$group)
pltdat$label = factor(pltdat$label, levels = unique(pltdat$label)[c(4,5,6,1,2,3)])

pltdat$anno[pltdat$anno=='HIGH'] = "High impact"
pltdat$anno[pltdat$anno=='MODERATE'] = "High/Moderate impact"
pltdat$anno[pltdat$anno=='LOW'] = "High/Moderate/Low impact"
pltdat$anno[pltdat$anno=='eqtl'] = "eQTL"
pltdat$anno[pltdat$anno=='sqtl'] = "sQTL"
pltdat$anno[pltdat$anno=='qtl'] = "e/sQTL"
pltdat$anno = factor(pltdat$anno, levels = c("High impact","High/Moderate impact","High/Moderate/Low impact","eQTL","sQTL","e/sQTL"))

pltdatsub = pltdat[pltdat$group!="SuSiE inference did not converge" & !is.na(pltdat$prop),]

options(repr.plot.width=3, repr.plot.height=6)
ggplot(pltdatsub[pltdatsub$anno %in% c("High impact","High/Moderate impact","High/Moderate/Low impact"),], aes(x = label, y = prop, color = label, shape = label)) +
  geom_point(size = 3) +
  facet_grid(anno~., scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
#  scale_color_manual(values = c("darkblue","red","lightsalmon")) +
      scale_color_manual(name = "", labels=c('SuSiE results\nin loci where\nSuSiE inference converged',
                                        'RSparsePro results\nin loci where\nSuSiE inference converged',
                                        'RSparsePro results\nin loci where\nSuSiE prior estimation failed'),
                                  values=c("darkblue","red","lightsalmon")) +
    scale_shape_manual(name = "", labels=c('SuSiE results\nin loci where\nSuSiE inference converged',
                                        'RSparsePro results\nin loci where\nSuSiE inference converged',
                                        'RSparsePro results\nin loci where\nSuSiE prior estimation failed'),
                                values=c(16,16,18)) +
  scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0,0.5))) + 
  labs(x = "", y = "Proportion of credible sets", color = "") -> plt
ggsave("../fig/dat_ldl_enrich_impact.pdf",plt,height = 6,width = 3)

options(repr.plot.width=3, repr.plot.height=6)
ggplot(pltdatsub[pltdatsub$anno %in% c("eQTL","sQTL","e/sQTL"),], aes(x = label, y = prop, color = label, shape = label)) +
  geom_point(size = 3) +
  facet_grid(anno~., scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
#  scale_color_manual(values = c("darkblue","red","lightsalmon")) +
      scale_color_manual(name = "", labels=c('SuSiE results\nin loci where\nSuSiE inference converged',
                                        'RSparsePro results\nin loci where\nSuSiE inference converged',
                                        'RSparsePro results\nin loci where\nSuSiE prior estimation failed'),
                                  values=c("darkblue","red","lightsalmon")) +
    scale_shape_manual(name = "", labels=c('SuSiE results\nin loci where\nSuSiE inference converged',
                                        'RSparsePro results\nin loci where\nSuSiE inference converged',
                                        'RSparsePro results\nin loci where\nSuSiE prior estimation failed'),
                                values=c(16,16,18)) +
  scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0,0.5))) + 
  labs(x = "", y = "Proportion of credible sets", color = "") -> plt
ggsave("../fig/dat_ldl_enrich_qtl.pdf",plt,height = 6,width = 3)
