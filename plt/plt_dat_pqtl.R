library(data.table)
library(ggplot2)

## Example
ld = as.matrix(fread("../dat/pQTL/bed/PCSK9_X5231_79.ld"))
pcsk9_fenland_rsp = read.table("../dat/pQTL/formatted/PCSK9_X5231_79_Fenland.txt.rsparsepro.txt", header = T)
pcsk9_fenland_sus = read.table("../dat/pQTL/formatted/PCSK9_X5231_79_Fenland.txt.susie.txt", header = T)
pcsk9_decode_rsp = read.table("../dat/pQTL/formatted/PCSK9_X5231_79_deCODE.txt.rsparsepro.txt", header = T)

pltdat = data.frame(pos = pcsk9_fenland_rsp$POS,
                    fen_pval = -log10(exp(1)) * pchisq(pcsk9_fenland_rsp$Z^2, df = 1, lower.tail = F, log.p = T),
                    dec_pval = -log10(exp(1)) * pchisq(pcsk9_decode_rsp$Z^2, df = 1, lower.tail = F, log.p = T),
                    cs_fen_sus = ifelse(pcsk9_fenland_sus$cs == 0, "Other variants",
                                paste0("Credible set ",pcsk9_fenland_sus$cs)),
                    cs_fen_rsp = ifelse(pcsk9_fenland_rsp$cs == 0, "Other variants",
                                paste0("Credible set ",pcsk9_fenland_rsp$cs)),
                    cs_dec_rsp = ifelse(pcsk9_decode_rsp$cs == 0, "Other variants",
                                paste0("Credible set ",pcsk9_decode_rsp$cs)),
                    cs_fen_sus_r1 = as.numeric(ld[which(pcsk9_fenland_sus$cs == 1)[which.max(pcsk9_fenland_sus$PIP[pcsk9_fenland_sus$cs == 1])],]^2),
                    cs_fen_sus_r2 = as.numeric(ld[which(pcsk9_fenland_sus$cs == 2)[which.max(pcsk9_fenland_sus$PIP[pcsk9_fenland_sus$cs == 2])],]^2),
                    cs_fen_sus_r3 = as.numeric(ld[which(pcsk9_fenland_sus$cs == 3)[which.max(pcsk9_fenland_sus$PIP[pcsk9_fenland_sus$cs == 3])],]^2),
                    cs_fen_rsp_r1 = as.numeric(ld[which(pcsk9_fenland_rsp$cs == 1)[which.max(pcsk9_fenland_rsp$PIP[pcsk9_fenland_rsp$cs == 1])],]^2),
                    cs_fen_rsp_r2 = as.numeric(ld[which(pcsk9_fenland_rsp$cs == 2)[which.max(pcsk9_fenland_rsp$PIP[pcsk9_fenland_rsp$cs == 2])],]^2),
                    cs_fen_rsp_r3 = as.numeric(ld[which(pcsk9_fenland_rsp$cs == 3)[which.max(pcsk9_fenland_rsp$PIP[pcsk9_fenland_rsp$cs == 3])],]^2),
                    cs_dec_rsp_r1 = as.numeric(ld[which(pcsk9_decode_rsp$cs == 1)[which.max(pcsk9_decode_rsp$PIP[pcsk9_decode_rsp$cs == 1])],]^2),
                    cs_dec_rsp_r2 = as.numeric(ld[which(pcsk9_decode_rsp$cs == 2)[which.max(pcsk9_decode_rsp$PIP[pcsk9_decode_rsp$cs == 2])],]^2),
                    cs_dec_rsp_r3 = as.numeric(ld[which(pcsk9_decode_rsp$cs == 3)[which.max(pcsk9_decode_rsp$PIP[pcsk9_decode_rsp$cs == 3])],]^2),
                    cs_dec_rsp_r4 = as.numeric(ld[which(pcsk9_decode_rsp$cs == 4)[which.max(pcsk9_decode_rsp$PIP[pcsk9_decode_rsp$cs == 4])],]^2))

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = fen_pval, color = cs_fen_sus)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) +
  xlim(54.9,56.1) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = "") +
  scale_color_manual(values = c("darkred","darkgreen","royalblue","grey")) +
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_fen_sus_cs.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = fen_pval, color = cs_fen_rsp)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) +
  xlim(54.9,56.1) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = "") +
  scale_color_manual(values = c("darkred","darkgreen","royalblue","grey")) +
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_fen_rsp_cs.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = dec_pval, color = cs_dec_rsp)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) +
  xlim(54.9,56.1) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = "") +
  scale_color_manual(values = c("darkred","darkgreen","royalblue","deeppink","grey")) +
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_dec_rsp_cs.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = fen_pval, color = cs_fen_sus_r1)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_fen_sus_cs1.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = fen_pval, color = cs_fen_sus_r2)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_fen_sus_cs2.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = fen_pval, color = cs_fen_sus_r3)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_fen_sus_cs3.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = fen_pval, color = cs_fen_rsp_r1)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_fen_rsp_cs1.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = fen_pval, color = cs_fen_rsp_r2)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_fen_rsp_cs2.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = fen_pval, color = cs_fen_rsp_r3)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_fen_rsp_cs3.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = dec_pval, color = cs_dec_rsp_r1)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_dec_rsp_cs2.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = dec_pval, color = cs_dec_rsp_r2)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_dec_rsp_cs1.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = dec_pval, color = cs_dec_rsp_r3)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_dec_rsp_cs3.pdf",plt,width = 10,height = 3)

options(repr.plot.width=10, repr.plot.height=3)
ggplot(pltdat, aes(x = pos/1000000, y = dec_pval, color = cs_dec_rsp_r4)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkgrey", lty = 2) +
  geom_point(size = 2) + 
  xlim(54.9,56.1) +
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
  xlab(paste0("Chr1 (Mb)")) -> plt
ggsave("../fig/dat_pqtl_dec_rsp_cs4.pdf",plt,width = 10,height = 3)

## Summary

stats = read.table("../doc/dat_pQTL_ncs_convrg_pairs.txt", sep='\t', header=T)
data = data.frame(
  category = c("SuSiE inference converged\n1,141 proteins from Fenland\n907 proteins from deCODE","SuSiE inference did not converge\n128 proteins from Fenland\n187 proteins from deCODE","SuSiE prior estimation failed\n84 proteins from Fenland\n259 proteins from deCODE"),
  count = c(sum(stats$fen_sus_convrg), sum(stats$fen_sus_exist)-sum(stats$fen_sus_convrg), sum(stats$fen_sus_exist==0),
            sum(stats$dec_sus_convrg), sum(stats$dec_sus_exist)-sum(stats$dec_sus_convrg), sum(stats$dec_sus_exist==0)),
  cohort = rep(c("Fenland","deCODE"),each = 3))
data$category = gsub("proteins","cis-pQTLs",data$category)

options(repr.plot.width=10, repr.plot.height=6)
ggplot(data[data$cohort=="Fenland",], aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("darkblue","steelblue1","grey")) +
  labs(fill = "1,353 cis-pQTLs from each cohort") -> plt
ggsave("../fig/dat_pqtl_pie_fen.pdf",plt,height = 6, width = 10)

options(repr.plot.width=10, repr.plot.height=6)
ggplot(data[data$cohort=="deCODE",], aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("darkblue","steelblue1","grey")) +
  labs(fill = "1,353 cis-pQTLs from each cohort") -> plt
ggsave("../fig/dat_pqtl_pie_dec.pdf",plt,height = 6, width = 10)

pltdat = data.frame(fen_ncs = c(stats$fen_rsp_ncs,
                            stats$fen_sus_ncs),
                dec_ncs = c(stats$dec_rsp_ncs,
                            stats$dec_sus_ncs),
                method = rep(c("RSparsePro","SuSiE"),each = nrow(stats)),
                group_dec = ifelse(stats$dec_sus_exist==0,"SuSiE prior estimation failed",
                                    ifelse(stats$dec_sus_convrg==1,"SuSiE inference converged","SuSiE inference did not converge")),
                group_fen = ifelse(stats$fen_sus_exist==0,"SuSiE prior estimation failed",
                                    ifelse(stats$fen_sus_convrg==1,"SuSiE inference converged","SuSiE inference did not converge")))
pltdat$label_fen = paste0(pltdat$method, " results\nin cis-pQTLs where\n", pltdat$group_fen)
pltdat$label_fen = factor(pltdat$label_fen, levels = unique(pltdat$label_fen)[c(4,5,6,1,3,2)])
pltdat$label_dec = paste0(pltdat$method, " results\nin cis-pQTLs where\n", pltdat$group_dec)
pltdat$label_dec = factor(pltdat$label_dec, levels = unique(pltdat$label_dec)[c(5,4,6,2,1,3)])

options(repr.plot.width=10, repr.plot.height=6)
ggplot(pltdat[pltdat$label_fen!="SuSiE results\nin cis-pQTLs where\nSuSiE prior estimation failed",], aes(x = label_fen, y = fen_ncs, color = label_fen)) +
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
ggsave("../fig/dat_pqtl_csnum_fen.pdf",plt,height=6,width=10)

options(repr.plot.width=10, repr.plot.height=6)
ggplot(pltdat[pltdat$label_dec!="SuSiE results\nin cis-pQTLs where\nSuSiE prior estimation failed",], aes(x = label_dec, y = dec_ncs, color = label_dec)) +
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
ggsave("../fig/dat_pqtl_csnum_dec.pdf",plt,height=6,width=10)

anno = read.table("../doc/dat_pQTL_anno_enrichment.txt", sep='\t', header=T)
pltdat = data.frame(prop = c(anno$fen_rsp_prop, anno$fen_sus_prop, anno$dec_rsp_prop, anno$dec_sus_prop), 
        anno = anno$anno, 
        group = anno$group, 
        method = rep(c("RSparsePro","SuSiE"), each = nrow(anno)),
        study = rep(c("Fenland", "deCODE"), each = 2*nrow(anno)))
pltdat$label = paste0(pltdat$method, " results\nin cis-pQTLs where\n", pltdat$group)
pltdat$label = factor(pltdat$label, levels = unique(pltdat$label)[c(4,5,6,1,3,2)])
pltdat$anno[pltdat$anno=='HIGH'] = "High impact"
pltdat$anno[pltdat$anno=='MODERATE'] = "High/Moderate impact"
pltdat$anno[pltdat$anno=='LOW'] = "High/Moderate/Low impact"
pltdat$anno[pltdat$anno=='eqtl'] = "eQTL"
pltdat$anno[pltdat$anno=='sqtl'] = "sQTL"
pltdat$anno[pltdat$anno=='qtl'] = "e/sQTL"
pltdat$anno = factor(pltdat$anno, levels = c("High impact","High/Moderate impact","High/Moderate/Low impact","eQTL","sQTL","e/sQTL"))
pltdat$study = factor(pltdat$study, levels = c("Fenland","deCODE"))

options(repr.plot.width=8, repr.plot.height=6)
ggplot(pltdat[pltdat$label!="SuSiE results\nin cis-pQTLs where\nSuSiE prior estimation failed" & pltdat$anno %in% c("High impact","High/Moderate impact","High/Moderate/Low impact"),], aes(x = label, y = prop, color = label, shape = label)) +
  geom_point(size = 4) +
  facet_grid(anno~study, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),strip.text = element_text(size = 13)) +
 # scale_color_manual(values = c("darkblue","steelblue1","red","gold1","lightsalmon")) +
      scale_color_manual(name = "", labels=c('SuSiE results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'SuSiE results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE prior estimation failed'),
                                  values=c("darkblue","steelblue1","red","gold1","lightsalmon")) +
    scale_shape_manual(name = "", labels=c('SuSiE results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'SuSiE results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE prior estimation failed'),
                                values=c(16,15,16,15,18)) +
  scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0,0.5))) + 
  labs(x = "", y = "Proportion of credible sets", color = "") -> plt
ggsave("../fig/dat_pqtl_enrich_impact.pdf",plt,height=8,width=6)

options(repr.plot.width=8, repr.plot.height=6)
ggplot(pltdat[pltdat$label!="SuSiE results\nin cis-pQTLs where\nSuSiE prior estimation failed" & pltdat$anno %in% c("eQTL","sQTL","e/sQTL"),], aes(x = label, y = prop, color = label, shape = label)) +
  geom_point(size = 4) +
  facet_grid(anno~study, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),strip.text = element_text(size = 13)) +
 # scale_color_manual(values = c("darkblue","steelblue1","red","gold1","lightsalmon")) +
      scale_color_manual(name = "", labels=c('SuSiE results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'SuSiE results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE prior estimation failed'),
                                  values=c("darkblue","steelblue1","red","gold1","lightsalmon")) +
    scale_shape_manual(name = "", labels=c('SuSiE results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'SuSiE results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE prior estimation failed'),
                                values=c(16,15,16,15,18)) +
  scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0,0.5))) + 
  labs(x = "", y = "Proportion of credible sets", color = "") -> plt
ggsave("../fig/dat_pqtl_enrich_qtl.pdf",plt,height=8,width=6)

rep = read.table("../doc/dat_pQTL_replication_rate.txt", sep='\t', header=T)
pltdat = data.frame(
    rate = c(rep$fen_rsp_matched, rep$fen_sus_matched, rep$dec_rsp_matched, rep$dec_sus_matched),
    group = rep$group,
    method = rep(c("RSparsePro","SuSiE"), each = nrow(rep)),
    study = rep(c("Fenland", "deCODE"), each = 2*nrow(rep)))
pltdat$label = paste0(pltdat$method, " results\nin cis-pQTLs where\n", pltdat$group)
pltdat$label = factor(pltdat$label, levels = unique(pltdat$label)[c(4,5,6,1,3,2)])
ggplot(pltdat[pltdat$label!="SuSiE results\nin cis-pQTLs where\nSuSiE prior estimation failed" & pltdat$study=="Fenland",], aes(x = label, y = rate, color = label, shape = label)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),strip.text = element_text(size = 13),
        title = element_text(size = 8)) +
  scale_color_manual(name = "", labels=c('SuSiE results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'SuSiE results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE prior estimation failed'),
                                  values=c("darkblue","steelblue1","red","gold1","lightsalmon")) +
    scale_shape_manual(name = "", labels=c('SuSiE results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'SuSiE results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE prior estimation failed'),
                                values=c(16,15,16,15,18)) +
  scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0,0.1))) + 
  labs(x = "", y = "", title = "Proportion of credible sets identified in Fenland replicated in deCODE", color = "") -> plt
ggsave("../fig/dat_pqtl_replicate_fen_dec.pdf",plt,height=2,width=5)

ggplot(pltdat[pltdat$label!="SuSiE results\nin cis-pQTLs where\nSuSiE prior estimation failed" & pltdat$study=="deCODE",], aes(x = label, y = rate, color = label, shape = label)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),strip.text = element_text(size = 13),
        title = element_text(size = 8)) +
  scale_color_manual(name = "", labels=c('SuSiE results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'SuSiE results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE prior estimation failed'),
                                  values=c("darkblue","steelblue1","red","gold1","lightsalmon")) +
    scale_shape_manual(name = "", labels=c('SuSiE results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'SuSiE results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference converged',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin cis-pQTLs where\nSuSiE prior estimation failed'),
                                values=c(16,15,16,15,18)) +
  scale_y_continuous(limits = c(0,NA),expand = expansion(mult = c(0,0.1))) + 
  labs(x = "", y = "", title = "Proportion of credible sets identified in deCODE replicated in Fenland", color = "") -> plt
ggsave("../fig/dat_pqtl_replicate_dec_fen.pdf",plt,height=2,width=5)

