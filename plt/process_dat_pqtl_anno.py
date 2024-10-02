import pandas as pd
import numpy as np
import pickle
import os

def process_anno(dir):
    anno_dict = {}
    rs_dict = {}
    sym_dict = {}
    with open(dir) as anno:
        header = next(anno)
        for line in anno:
            ll = dict(zip(header.strip().split(), line.strip().split()))
            if ll['SYMBOL'] != '-':
                ch, pos = ll['Location'].split(':')
                var_id = 'chr{}_{}'.format(ch, pos.split('-')[0])
                if ll['#Uploaded_variation'] not in rs_dict:
                    rs_dict[ll['#Uploaded_variation']] = var_id
                if ll['SYMBOL'] not in sym_dict:
                    sym_dict[ll['SYMBOL']] = ll['Gene']
                if var_id not in anno_dict:
                    anno_dict[var_id] = {}
                if ll['SYMBOL'] not in anno_dict[var_id]:
                    anno_dict[var_id][ll['SYMBOL']] = set()
                anno_dict[var_id][ll['SYMBOL']].add(ll['IMPACT'])
    return anno_dict, rs_dict, sym_dict

anno_dict, rs_dict, sym_dict = process_anno('../doc/dat_pQTL_vep_annotated.txt')

with open('../dat/GTEx/sqtl_dict.pkl', 'rb') as file:
    sqtl_dict = pickle.load(file)

with open('../dat/GTEx/eqtl_dict.pkl', 'rb') as file:
    eqtl_dict = pickle.load(file)

def get_enrich_anno_qtl(anno_dict, rs_dict, sym_dict, eqtl_dict, sqtl_dict, path):
    df = pd.read_csv(path, sep='\t')
    df['gene'] = df['loci'].apply(lambda x: x.split('_')[0])
    df['anno'] = df.apply(lambda row: anno_dict.get(rs_dict.get(row['RSID']), {}).get(row['gene'], None),axis=1)
    df['HIGH'] = df['anno'].apply(lambda x: 'HIGH' in x if x is not None else False)
    df['MODERATE'] = df['anno'].apply(lambda x: 'HIGH' in x or 'MODERATE' in x if x is not None else False)
    df['LOW'] = df['anno'].apply(lambda x: 'HIGH' in x or 'MODERATE' in x or 'LOW' in x if x is not None else False)
    df['eqtl'] = df.apply(lambda row: 1 if eqtl_dict.get(rs_dict.get(row['RSID']), {}).get(sym_dict.get(row['gene'])) is not None else 0, axis=1)
    df['sqtl'] = df.apply(lambda row: 1 if sqtl_dict.get(rs_dict.get(row['RSID']), {}).get(sym_dict.get(row['gene'])) is not None else 0, axis=1)
    df['qtl'] = df.apply(lambda row: 1 if row['eqtl'] or row['sqtl'] else 0, axis=1)
    return df

df_fen_sus_anno = get_enrich_anno_qtl(anno_dict, rs_dict, sym_dict, eqtl_dict, sqtl_dict, '../doc/dat_pQTL_fen_sus_cs.txt')
df_fen_rsp_anno = get_enrich_anno_qtl(anno_dict, rs_dict, sym_dict, eqtl_dict, sqtl_dict, '../doc/dat_pQTL_fen_rsp_cs.txt')
df_dec_sus_anno = get_enrich_anno_qtl(anno_dict, rs_dict, sym_dict, eqtl_dict, sqtl_dict, '../doc/dat_pQTL_dec_sus_cs.txt')
df_dec_rsp_anno = get_enrich_anno_qtl(anno_dict, rs_dict, sym_dict, eqtl_dict, sqtl_dict, '../doc/dat_pQTL_dec_rsp_cs.txt')

df_fen_sus_anno.to_csv('../doc/dat_pQTL_fen_sus_cs_anno.txt', sep='\t', header=True, index=False)
df_fen_rsp_anno.to_csv('../doc/dat_pQTL_fen_rsp_cs_anno.txt', sep='\t', header=True, index=False)
df_dec_sus_anno.to_csv('../doc/dat_pQTL_dec_sus_cs_anno.txt', sep='\t', header=True, index=False)
df_dec_rsp_anno.to_csv('../doc/dat_pQTL_dec_rsp_cs_anno.txt', sep='\t', header=True, index=False)