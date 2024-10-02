import pandas as pd
import numpy as np

df_summary = pd.read_csv("../doc/dat_pQTL_ncs_convrg_pairs.txt", sep='\t')
df_summary['fen_group'] = np.select([df_summary['fen_sus_exist'] == 0, df_summary['fen_sus_convrg'] == 1], ["SuSiE prior estimation failed", "SuSiE inference converged"], default="SuSiE inference did not converge")
df_summary['dec_group'] = np.select([df_summary['dec_sus_exist'] == 0, df_summary['dec_sus_convrg'] == 1], ["SuSiE prior estimation failed", "SuSiE inference converged"], default="SuSiE inference did not converge")
df_summary.index = df_summary['id']

df_rsp_fen_anno = pd.read_csv('../doc/dat_pQTL_fen_rsp_cs_anno.txt', sep='\t')
df_rsp_fen_anno['group'] = [df_summary.loc[i.replace("_Fenland.txt.rsparsepro.txt", ""), 'fen_group'] for i in df_rsp_fen_anno['loci']]
df_sus_fen_anno = pd.read_csv('../doc/dat_pQTL_fen_sus_cs_anno.txt', sep='\t')
df_sus_fen_anno['group'] = [df_summary.loc[i.replace("_Fenland.txt.susie.txt", ""), 'fen_group'] for i in df_sus_fen_anno['loci']]
df_rsp_dec_anno = pd.read_csv('../doc/dat_pQTL_dec_rsp_cs_anno.txt', sep='\t')
df_rsp_dec_anno['group'] = [df_summary.loc[i.replace("_deCODE.txt.rsparsepro.txt", ""), 'dec_group'] for i in df_rsp_dec_anno['loci']]
df_sus_dec_anno = pd.read_csv('../doc/dat_pQTL_dec_sus_cs_anno.txt', sep='\t')
df_sus_dec_anno['group'] = [df_summary.loc[i.replace("_deCODE.txt.susie.txt", ""), 'dec_group'] for i in df_sus_dec_anno['loci']]

res_lst = []
for c in ["SuSiE inference converged", "SuSiE prior estimation failed", "SuSiE inference did not converge"]:
        for a in ["HIGH", "MODERATE", "LOW", "eqtl", "sqtl", "qtl"]:
            fen_rsp_prop = np.mean(df_rsp_fen_anno[df_rsp_fen_anno['group']==c].groupby(['cs', 'loci'])[a].sum()!=0)
            fen_sus_prop = np.mean(df_sus_fen_anno[df_sus_fen_anno['group']==c].groupby(['cs', 'loci'])[a].sum()!=0)
            dec_rsp_prop = np.mean(df_rsp_dec_anno[df_rsp_dec_anno['group']==c].groupby(['cs', 'loci'])[a].sum()!=0)
            dec_sus_prop = np.mean(df_sus_dec_anno[df_sus_dec_anno['group']==c].groupby(['cs', 'loci'])[a].sum()!=0)
            res_lst.append({"group": c, "anno": a, "fen_rsp_prop": fen_rsp_prop, "fen_sus_prop": fen_sus_prop, "dec_rsp_prop": dec_rsp_prop, "dec_sus_prop": dec_sus_prop})

pd.DataFrame(res_lst).to_csv("../doc/dat_pQTL_anno_enrichment.txt", sep='\t', header=True, index=False)

rep_lst = []
for c in ["SuSiE inference converged", "SuSiE prior estimation failed", "SuSiE inference did not converge"]:
    fen_rsp_matched = df_summary.loc[df_summary['fen_group']==c]['rsp_matched'].sum()/df_summary.loc[df_summary['fen_group']==c]['fen_rsp_ncs'].sum()
    fen_sus_matched = df_summary.loc[df_summary['fen_group']==c]['sus_matched'].sum()/df_summary.loc[df_summary['fen_group']==c]['fen_sus_ncs'].sum()
    dec_rsp_matched = df_summary.loc[df_summary['dec_group']==c]['rsp_matched'].sum()/df_summary.loc[df_summary['dec_group']==c]['dec_rsp_ncs'].sum()
    dec_sus_matched = df_summary.loc[df_summary['dec_group']==c]['sus_matched'].sum()/df_summary.loc[df_summary['dec_group']==c]['dec_sus_ncs'].sum()
    rep_lst.append({"group": c,
                    "fen_rsp_matched": fen_rsp_matched,
                    "fen_sus_matched": fen_sus_matched,
                    "dec_rsp_matched": dec_rsp_matched,
                    "dec_sus_matched": dec_sus_matched})

pd.DataFrame(rep_lst).to_csv("../doc/dat_pQTL_replication_rate.txt", sep='\t', header=True, index=False)