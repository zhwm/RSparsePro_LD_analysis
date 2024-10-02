import pandas as pd
import numpy as np
import os

def get_summary(row):
    chr, start, end, ans = row
    rsp_file = '../dat/LDL/formatted/{}_{}_{}_{}.txt.rsparsepro.txt'.format(ans, chr, start, end)
    sus_file = '../dat/LDL/formatted/{}_{}_{}_{}.txt.susie.txt'.format(ans, chr, start, end)
    rsp = pd.read_csv(rsp_file, sep='\t')
    if (rsp['cs']==0).all():
        rsp_cs = None
        rsp_ncs = 0
    else:
        rsp_ncs = rsp['cs'].max()
        rsp_cs = rsp.loc[rsp['cs']!=0].copy()
        rsp_cs['loci'] = rsp_file.split('/')[4]
    if os.path.exists(sus_file):
        sus_exist = 1
        sus = pd.read_csv(sus_file, sep='\t')
        sus_ncs = sus['cs'].max()
        sus_cs = sus.loc[sus['cs']!=0].copy()
        sus_cs['loci'] = sus_file.split('/')[4]
        if sus['converged'][0]:
            sus_convrg = 1
        else:
            sus_convrg = 0
    else:
        sus_exist = 0
        sus_convrg = 0
        sus_ncs = 0
        sus_cs = None
    return rsp_ncs, sus_ncs, sus_exist, sus_convrg, rsp_cs, sus_cs

afr_lst = pd.read_csv('../dat/LDL/loci/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.leads.txt.gtf', sep='\t', header=None)
afr_lst['ancestry'] = 'AFR'
amr_lst = pd.read_csv('../dat/LDL/loci/LDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz.leads.txt.gtf', sep='\t', header=None)
amr_lst['ancestry'] = 'AMR'
eas_lst = pd.read_csv('../dat/LDL/loci/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz.leads.txt.gtf', sep='\t', header=None)
eas_lst['ancestry'] = 'EAS'
eur_lst = pd.read_csv('../dat/LDL/loci/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.leads.txt.gtf', sep='\t', header=None)
eur_lst['ancestry'] = 'EUR'
sas_lst = pd.read_csv('../dat/LDL/loci/LDL_INV_SAS_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.leads.txt.gtf', sep='\t', header=None)
sas_lst['ancestry'] = 'SAS'
all_lst = pd.concat([afr_lst, amr_lst, eas_lst, eur_lst, sas_lst])
all_lst.columns = ['chr', 'start', 'end', 'ancestry']

rsp_ncs_lst = []
sus_ncs_lst = []
sus_exist_lst = []
sus_convrg_lst = []
rsp_cs_lst = []
sus_cs_lst = []
for idx, row in all_lst.iterrows():
    rsp_ncs, sus_ncs, sus_exist, sus_convrg, rsp_cs, sus_cs = get_summary(row)
    rsp_ncs_lst.append(rsp_ncs)
    sus_ncs_lst.append(sus_ncs)
    sus_exist_lst.append(sus_exist)
    sus_convrg_lst.append(sus_convrg)
    rsp_cs_lst.append(rsp_cs)
    sus_cs_lst.append(sus_cs)

all_lst['rsp_ncs'] = rsp_ncs_lst
all_lst['sus_ncs'] = sus_ncs_lst
all_lst['sus_exist'] = sus_exist_lst
all_lst['sus_convrg'] = sus_convrg_lst

all_lst.to_csv('../doc/dat_LDL_ncs_convrg.txt', sep='\t', index=False, header=True)

df_rsp_cs = pd.concat(rsp_cs_lst)
df_rsp_cs.to_csv('../doc/dat_LDL_rsp_cs.txt', sep='\t', index=False, header=True)
df_sus_cs = pd.concat(sus_cs_lst)
df_sus_cs.to_csv('../doc/dat_LDL_sus_cs.txt', sep='\t', index=False, header=True)
vep_rsid = set()
for df in [df_rsp_cs, df_sus_cs]:
    vep_rsid.update(df['RSID'].values)
with open('../doc/dat_LDL_vep.txt', 'w') as file:
    for rsid in vep_rsid:
        file.write(f"{rsid}\n")