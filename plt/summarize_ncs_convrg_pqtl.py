import pandas as pd
import numpy as np
import os

def get_pairs(df1, df2):
    if df1 is None or len(df1)==0 or df2 is None or len(df2)==0:
        matched_cs = 0
    else:
        matched_cs = pd.merge(df1[['RSID', 'cs']], df2[['RSID', 'cs']], on='RSID', how='inner')['cs_x'].nunique()
    return matched_cs

def load_file_rsp(dir):
    rsp = pd.read_csv(dir, sep='\t')
    if (rsp['cs']==0).all():
        rsp_cs = None
        
        ncs = 0
    else:
        rsp_cs = rsp.loc[rsp['cs']!=0].copy()
        rsp_cs['loci'] = dir.split('/')[4]
        ncs = rsp_cs['cs'].max()
    return rsp_cs, ncs

def load_file_sus(dir):
    if os.path.exists(dir):
        sus_exist = 1
        df_cs, ncs = load_file_rsp(dir)
        if df_cs['converged'].all():
            sus_convrg = 1
        else:
            sus_convrg = 0
    else:
        df_cs, ncs = None, 0
        sus_exist = 0
        sus_convrg = 0
    return sus_exist, sus_convrg, df_cs, ncs

def get_data(dec, fen):
    fen_rsp_file = '../dat/pQTL/{}.rsparsepro.txt'.format(fen)
    dec_rsp_file = '../dat/pQTL/{}.rsparsepro.txt'.format(dec)
    fen_sus_file = '../dat/pQTL/{}.susie.txt'.format(fen)
    dec_sus_file = '../dat/pQTL/{}.susie.txt'.format(dec)
    fen_rsp_cs, fen_rsp_ncs = load_file_rsp(fen_rsp_file)
    dec_rsp_cs, dec_rsp_ncs = load_file_rsp(dec_rsp_file)
    fen_sus_exist, fen_sus_convrg, fen_sus_cs, fen_sus_ncs = load_file_sus(fen_sus_file)
    dec_sus_exist, dec_sus_convrg, dec_sus_cs, dec_sus_ncs = load_file_sus(dec_sus_file)
    rsp_matched = get_pairs(fen_rsp_cs, dec_rsp_cs)
    sus_matched = get_pairs(fen_sus_cs, dec_sus_cs)
    summary = {'id': '_'.join(dec.split('/')[1].split('_')[0:3]),
                'fen_rsp_ncs': fen_rsp_ncs,
                'dec_rsp_ncs': dec_rsp_ncs,
                'fen_sus_exist': fen_sus_exist,
                'fen_sus_convrg': fen_sus_convrg,
                'fen_sus_ncs': fen_sus_ncs,
                'dec_sus_exist': dec_sus_exist,
                'dec_sus_convrg': dec_sus_convrg,
                'dec_sus_ncs': dec_sus_ncs,
                'rsp_matched': rsp_matched,
                'sus_matched': sus_matched}
    return summary, fen_rsp_cs, dec_rsp_cs, fen_sus_cs, dec_sus_cs

lst = pd.read_csv('../dat/pQTL/paired.lst.all', sep='\t', header=None)

summary_lst = []
dec_rsp_cs_lst = []
fen_rsp_cs_lst = []
dec_sus_cs_lst = []
fen_sus_cs_lst = []
for idx in range(len(lst)):
    print(idx)
    summary, fen_rsp_cs, dec_rsp_cs, fen_sus_cs, dec_sus_cs = get_data(lst.iloc[idx,6], lst.iloc[idx,7])
    summary_lst.append(summary)
    dec_rsp_cs_lst.append(dec_rsp_cs)
    fen_rsp_cs_lst.append(fen_rsp_cs)
    dec_sus_cs_lst.append(dec_sus_cs)
    fen_sus_cs_lst.append(fen_sus_cs)

df_summary = pd.DataFrame(summary_lst)
df_summary.to_csv('../doc/dat_pQTL_ncs_convrg_pairs.txt', sep='\t', index=False, header=True)

df_dec_rsp_cs = pd.concat(dec_rsp_cs_lst)
df_dec_rsp_cs.to_csv('../doc/dat_pQTL_dec_rsp_cs.txt', sep='\t', index=False, header=True)
df_fen_rsp_cs = pd.concat(fen_rsp_cs_lst)
df_fen_rsp_cs.to_csv('../doc/dat_pQTL_fen_rsp_cs.txt', sep='\t', index=False, header=True)
df_dec_sus_cs = pd.concat(dec_sus_cs_lst)
df_dec_sus_cs.to_csv('../doc/dat_pQTL_dec_sus_cs.txt', sep='\t', index=False, header=True)
df_fen_sus_cs = pd.concat(fen_sus_cs_lst)
df_fen_sus_cs.to_csv('../doc/dat_pQTL_fen_sus_cs.txt', sep='\t', index=False, header=True)
vep_rsid = set()
for df in [df_dec_rsp_cs, df_fen_rsp_cs, df_dec_sus_cs, df_fen_sus_cs]:
    vep_rsid.update(df['RSID'].values)
with open('../doc/dat_pQTL_vep.txt', 'w') as file:
    for rsid in vep_rsid:
        file.write(f"{rsid}\n")