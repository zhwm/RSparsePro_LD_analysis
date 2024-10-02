import pandas as pd
import numpy as np

df_summary = pd.read_csv("../doc/dat_LDL_ncs_convrg.txt", sep='\t')
df_summary['group'] = np.select([df_summary['sus_exist'] == 0, df_summary['sus_convrg'] == 1], ["SuSiE prior estimation failed", "SuSiE inference converged"], default="SuSiE inference did not converge")
df_summary.index = df_summary['ancestry'] + '_' + df_summary['chr'].astype(str) + '_' + df_summary['start'].astype(str) + '_' + df_summary['end'].astype(str) 
df_rsp_cs_anno = pd.read_csv('../doc/dat_LDL_rsp_cs_anno.txt', sep='\t')
df_rsp_cs_anno['group'] = [df_summary.loc[i.replace(".txt.rsparsepro.txt", ""), 'group'] for i in df_rsp_cs_anno['loci']]
df_sus_cs_anno = pd.read_csv('../doc/dat_LDL_sus_cs_anno.txt', sep='\t')
df_sus_cs_anno['group'] = [df_summary.loc[i.replace(".txt.susie.txt", ""), 'group'] for i in df_sus_cs_anno['loci']]
res_lst = []
for c in ["SuSiE inference converged", "SuSiE prior estimation failed", "SuSiE inference did not converge"]:
        for a in ["HIGH", "MODERATE", "LOW", "eqtl", "sqtl", "qtl"]:
            sus_prop = np.mean(df_sus_cs_anno[df_sus_cs_anno['group']==c].groupby(['cs', 'loci'])[a].sum()!=0)
            rsp_prop = np.mean(df_rsp_cs_anno[df_rsp_cs_anno['group']==c].groupby(['cs', 'loci'])[a].sum()!=0)
            res_lst.append({"group": c, "anno": a, "sus_prop": sus_prop, "rsp_prop": rsp_prop})

pd.DataFrame(res_lst).to_csv("../doc/dat_LDL_anno_enrichment.txt", sep='\t', header=True, index=False)