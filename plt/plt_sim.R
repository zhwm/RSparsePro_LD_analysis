library(data.table)
library(ggplot2)
library(ggrepel)

cat_summary = read.table('../doc/sim_pwr_cov_categorized.txt', sep='\t', header=T)
cat_summary$label = paste0(cat_summary$Method, " results\nin simulated replicates where\n", cat_summary$Category)
cat_summary$flip = factor(paste0("r = ", cat_summary$Deviation_Factor), levels = paste0("r = ",c(0.5,0,-0.5,-1)))
cat_summary$prop = factor(paste0("p = ", cat_summary$Mismatch_Proportion), levels = paste0("p = ",c(0,0.001,0.01,0.1)))

sumstats5 = cat_summary[(cat_summary$Number_of_Causal_Variants==5) & (cat_summary$Deviation_Factor!=1) & (cat_summary$Cases!=0) & (cat_summary$label!='SuSiE results\nin simulated replicates where\nSuSiE prior estimation failed'),]
sumstats3 = cat_summary[(cat_summary$Number_of_Causal_Variants==3) & (cat_summary$Deviation_Factor!=1) & (cat_summary$Cases!=0) & (cat_summary$label!='SuSiE results\nin simulated replicates where\nSuSiE prior estimation failed'),]
sumstats5$label = factor(sumstats5$label, levels = c('SuSiE results\nin simulated replicates where\nSuSiE inference converged', 'SuSiE results\nin simulated replicates where\nSuSiE inference did not converge', 'RSparsePro results\nin simulated replicates where\nSuSiE inference converged', 'RSparsePro results\nin simulated replicates where\nSuSiE inference did not converge', 'RSparsePro results\nin simulated replicates where\nSuSiE prior estimation failed'))
sumstats3$label = factor(sumstats3$label, levels = c('SuSiE results\nin simulated replicates where\nSuSiE inference converged', 'SuSiE results\nin simulated replicates where\nSuSiE inference did not converge', 'RSparsePro results\nin simulated replicates where\nSuSiE inference converged', 'RSparsePro results\nin simulated replicates where\nSuSiE inference did not converge', 'RSparsePro results\nin simulated replicates where\nSuSiE prior estimation failed'))

options(repr.plot.width=16, repr.plot.height=10)
set.seed(4)
ggplot(sumstats3, aes(x = Power, y = Coverage, color=label, shape=label, label=Cases)) +
    geom_point(size = 4) +
    geom_text_repel(show.legend = F, size = 6, nudge_x = -0.05, nudge_y = 0.05) +
    theme_bw() +
    theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14)) +
    scale_color_manual(name = "", labels=c('SuSiE results\nin simulated replicates where\nSuSiE inference converged',
                                        'SuSiE results\nin simulated replicates where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE inference converged',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE prior estimation failed'),
                                  values=c("darkblue","steelblue1","red","gold1","lightsalmon")) +
    scale_shape_manual(name = "", labels=c('SuSiE results\nin simulated replicates where\nSuSiE inference converged',
                                        'SuSiE results\nin simulated replicates where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE inference converged',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE prior estimation failed'),
                                values=c(16,15,16,15,18)) +
    #scale_color_manual(values = c("darkblue","steelblue1","red","gold1","lightsalmon")) +
    #scale_shape_manual(values = c(16,15,18)) +
    labs(x = "Power", y = "Coverage", color = "", shape="") +
    facet_grid(prop ~ flip) -> plt 
ggsave(paste0("../fig/sim_pwr_cov_3.pdf"), plt, height = 10, width = 16)

options(repr.plot.width=16, repr.plot.height=10)
set.seed(4)
ggplot(sumstats5, aes(x = Power, y = Coverage, color=label, shape=label, label=Cases)) +
    geom_point(size = 4) +
    geom_text_repel(show.legend = F, size = 6, nudge_x = -0.05, nudge_y = 0.05) +
    theme_bw() +
    theme(legend.position = "bottom",
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14)) +
    scale_color_manual(name = "", labels=c('SuSiE results\nin simulated replicates where\nSuSiE inference converged',
                                        'SuSiE results\nin simulated replicates where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE inference converged',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE prior estimation failed'),
                                  values=c("darkblue","steelblue1","red","gold1","lightsalmon")) +
    scale_shape_manual(name = "", labels=c('SuSiE results\nin simulated replicates where\nSuSiE inference converged',
                                        'SuSiE results\nin simulated replicates where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE inference converged',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE inference did not converge',
                                        'RSparsePro results\nin simulated replicates where\nSuSiE prior estimation failed'),
                                values=c(16,15,16,15,18)) +
    #scale_color_manual(values = c("darkblue","steelblue1","red","gold1","lightsalmon")) +
    #scale_shape_manual(values = c(16,15,18)) +
    labs(x = "Power", y = "Coverage", color = "", shape="") +
    facet_grid(prop ~ flip) -> plt 
ggsave(paste0("../fig/sim_pwr_cov_5.pdf"), plt, height = 10, width = 16)

summary = read.table('../doc/sim_pwr_cov.txt', sep='\t', header=T)
summary$flip = factor(paste0("r = ", summary$Deviation_Factor), levels = paste0("r = ",c(0.5,0,-0.5,-1)))
summary$prop = factor(paste0("p = ", summary$Mismatch_Proportion), levels = paste0("p = ",c(0,0.001,0.01,0.1)))

summary_mismatch3 = summary[(summary$Deviation_Factor!=1 & summary$Number_of_Causal_Variants==3),]

options(repr.plot.width=4, repr.plot.height=7)
ggplot(summary_mismatch3, aes(x = flip, y = SuSiE_Convergence_Rate)) +
  geom_point(size = 3, color = "darkblue") +
  facet_grid(prop~.) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ylab("Proportion of simulated replicates in which SuSiE inference converged") +
  xlab("") -> plt
ggsave(paste0("../fig/sim_conv_rate_3.pdf"), plt, height = 7, width = 4)

summary_mismatch5 = summary[(summary$Deviation_Factor!=1 & summary$Number_of_Causal_Variants==5),]

ggplot(summary_mismatch5, aes(x = flip, y = SuSiE_Convergence_Rate)) +
  geom_point(size = 3, color = "darkblue") +
  facet_grid(prop~.) +
  theme_bw() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14)) +
  ylab("Proportion of simulated replicates in which SuSiE inference converged") +
  xlab("") -> plt
ggsave(paste0("../fig/sim_conv_rate_5.pdf"), plt, height = 7, width = 4)

df_match = read.table('../doc/sim_match.txt', sep='\t', header=T)
df_mismatch = read.table('../doc/sim_mismatch.txt', sep='\t', header=T)

plt_mismatch = data.frame(ncs = c(df_mismatch$sus_ncs, df_mismatch$rsp_ncs), 
                        csz = c(df_mismatch$sus_csz, df_mismatch$rsp_csz), 
                        k = df_mismatch$kc,
                        group = ifelse(df_mismatch$exist==0, "SuSiE prior estimation failed",
                                    ifelse(df_mismatch$convrg==1,"SuSiE inference converged","SuSiE inference did not converge")),
                        flip = sapply(df_mismatch$flip, function(x) as.numeric(gsub("f", "", x))),
                        prop = sapply(df_mismatch$prop, function(x) as.numeric(gsub("p", "", x))),
                        method = rep(c("SuSiE","RSparsePro"), each = nrow(df_mismatch)))

plt_mismatch$flip = factor(paste0("r = ", plt_mismatch$flip), levels = paste0("r = ",c(0.5,0,-0.5,-1)))
plt_mismatch$prop = factor(paste0("p = ", plt_mismatch$prop), levels = paste0("p = ",c(0,0.001,0.01,0.1)))

plt_mismatch$label = paste0(plt_mismatch$method, " results\nin simulated replicates where\n", plt_mismatch$group)
plt_mismatch = plt_mismatch[plt_mismatch$label!="SuSiE results\nin simulated replicates where\nSuSiE prior estimation failed",]
plt_mismatch$label = factor(plt_mismatch$label, levels = c('SuSiE results\nin simulated replicates where\nSuSiE inference converged', 'SuSiE results\nin simulated replicates where\nSuSiE inference did not converge', 'RSparsePro results\nin simulated replicates where\nSuSiE inference converged', 'RSparsePro results\nin simulated replicates where\nSuSiE inference did not converge', 'RSparsePro results\nin simulated replicates where\nSuSiE prior estimation failed'))

options(repr.plot.width=18, repr.plot.height=10)
ggplot(plt_mismatch[plt_mismatch$k==3,], aes(x = label, y = ncs, color = label)) +
  geom_hline(yintercept = 3, lty = 2, col = "darkgrey") +
  geom_violin(fill = "white") +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("darkblue","steelblue1","red","gold1","lightsalmon")) +
  labs(x = "", y = "Number of identified credible sets", color = "") +
  scale_y_continuous(breaks = seq(0,10,by = 2)) +
  facet_grid(prop ~ flip) -> plt 
ggsave(paste0("../fig/sim_ncs_3.pdf"), plt, height = 10, width = 16)

options(repr.plot.width=18, repr.plot.height=10)
ggplot(plt_mismatch[plt_mismatch$k==5,], aes(x = label, y = ncs, color = label)) +
  geom_hline(yintercept = 5, lty = 2, col = "darkgrey") +
  geom_violin(fill = "white") +
  geom_boxplot(width = 0.1, fill = "white") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_color_manual(values = c("darkblue","steelblue1","red","gold1","lightsalmon")) +
  labs(x = "", y = "Number of identified credible sets", color = "") +
  scale_y_continuous(breaks = seq(0,10,by = 2)) +
  facet_grid(prop ~ flip) -> plt 
ggsave(paste0("../fig/sim_ncs_5.pdf"), plt, height = 10, width = 16)

## Example
library(data.table)
library(ggplot2)

ld = as.matrix(fread("../sim/Locus1.ld"))
btrue = fread("../sim/Locus1-3-0.0001-500000-btrue.txt")
zfile = fread("../sim/Locus1-3-0.0001-500000-z.txt")
zdiff = fread("../sim/Locus1-3-0.0001-500000-z.zdiff.txt")
zcs = fread("../sim/Locus1-3-0.0001-500000-z.cs.txt")
ztrue = fread("../sim/Locus1-3-0.0001-500000-ztrue.txt")

btruevec = btrue$ite16!=0
ldidx = which(btruevec)[2]
ztruevec = ztrue$ite16
zhatvec =zfile$`f-0.5_p0.1_ite16`
zdiffvec = zhatvec-zdiff$`f-0.5_p0.1_ite16`
matchvec = ifelse(zhatvec==ztruevec, 'Matched', 'Mismatched')
csvec = zcs$`f-0.5_p0.1_ite16`

pltdat = data.frame(idx = 1:nrow(btrue),
                    pval = c(-log10(exp(1)) * pchisq(ztruevec^2, df = 1, lower.tail = F, log.p = T), 
                            -log10(exp(1)) * pchisq(zhatvec^2, df = 1, lower.tail = F, log.p = T),
                            -log10(exp(1)) * pchisq(zdiffvec^2, df = 1, lower.tail = F, log.p = T)),
                    mismatch = c(rep('Matched', nrow(btrue)), rep(matchvec, 2)),
                    category = rep(c("GWAS summary statistics", "GWAS summary statistics with LD mismatch", "Estimated GWAS summary statistics"),each = nrow(btrue)),
                    r2 = as.numeric(ld[,ldidx])^2)

pltdat$category = factor(pltdat$category, levels = c("GWAS summary statistics", "GWAS summary statistics with LD mismatch", "Estimated GWAS summary statistics"))

options(repr.plot.width=6, repr.plot.height=6)
ggplot(pltdat[pltdat$category!="Estimated GWAS summary statistics",], aes(x = idx, y = pval, color = r2, shape=mismatch)) +
  facet_wrap(category~.,nrow = 3, strip.position = "top") +
  geom_point(size = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") +
  ylab(expression(-log[10](p-value))) +
  labs(color = expression(r^2)) +
  scale_color_stepsn(n.breaks = 6, colours = c("black","blue","green","orange","red")) +
  scale_shape_manual(values=c(16, 4)) +
  xlab("") -> plt 
ggsave(paste0("../fig/sim_example_gwas.pdf"), plt, height = 6, width = 6)

scplt = data.frame(ztrue = ztruevec, zhat = zhatvec, zdiff = zdiffvec, mismatch = matchvec, r2 = as.numeric(ld[,ldidx])^2)
scplt_ordered = scplt[order(scplt$r2), ]
options(repr.plot.width=3, repr.plot.height=3)
ggplot(scplt_ordered, aes(x = ztrue, y = zhat, shape = mismatch, color = r2)) + 
    geom_point(size = 2) + 
    theme_classic() + 
    scale_color_stepsn(n.breaks = 6, colours = c("black", "blue", "green", "orange", "red")) +
    labs(color = expression(r^2), shape = "") +
    scale_shape_manual(values = c(16, 4)) + 
    xlab("z-scores") + 
    ylab("z-scores with LD mismatch") +
    theme(axis.title = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") -> plt 
ggsave(paste0("../fig/sim_example_mismatch.pdf"), plt, height = 3, width = 3)

scplt$cs = "others"
for (i in 1:max(csvec)) {
    scplt$cs[csvec==i] = paste0("credible set ", i)
}
scplt_cs_order = scplt[order(scplt$cs, decreasing = T),]
options(repr.plot.width=3, repr.plot.height=3)
ggplot(scplt_cs_order, aes(x = ztrue, y = zdiff, shape=mismatch, color=cs)) + 
    geom_point(size = 2) + 
    theme_classic() + 
    scale_color_manual(values = c("darkred","darkgreen","royalblue", "grey")) +
    #scale_color_stepsn(n.breaks = 6, colours = c("black","blue","green","orange","red")) +
    labs(color = "", shape="") +
    scale_shape_manual(values=c(16, 4)) + 
    xlab("z-scores") + 
    ylab("Estimated z-scores")  +
    theme(axis.title = element_text(size = 14),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.position = "none") -> plt 
ggsave(paste0("../fig/sim_example_cs.pdf"), plt, height = 3, width = 3)