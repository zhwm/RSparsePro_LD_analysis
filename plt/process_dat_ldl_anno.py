import pandas as pd
import numpy as np
import pickle

anno_dict = {}
rs_dict = {}

with open('../doc/dat_LDL_vep_annotated.txt') as anno:
    header = next(anno)
    for line in anno:
        ll = line.strip().split()
        ch, pos = ll[1].split(':')
        if ll[0] not in rs_dict:
            rs_dict[ll[0]] = 'chr{}_{}'.format(ch, pos.split('-')[0])
        if ll[0] in anno_dict:
            if ll[4] not in anno_dict[ll[0]]:
                anno_dict[ll[0]].append(ll[4])
        else:
            anno_dict[ll[0]] = [ll[4]]

with open('../dat/GTEx/sqtl_dict.pkl', 'rb') as file:
    sqtl_dict = pickle.load(file)

with open('../dat/GTEx/eqtl_dict.pkl', 'rb') as file:
    eqtl_dict = pickle.load(file)

def get_enrich_anno_qtl(anno_dict, rs_dict, eqtl_dict, sqtl_dict, path):
    df = pd.read_csv(path, sep='\t')
    df['anno'] = df['RSID'].apply(lambda x: anno_dict.get(x))
    df['HIGH'] = df['anno'].apply(lambda x: 'HIGH' in x if x is not None else False)
    df['MODERATE'] = df['anno'].apply(lambda x: 'HIGH' in x or 'MODERATE' in x if x is not None else False)
    df['LOW'] = df['anno'].apply(lambda x: 'HIGH' in x or 'MODERATE' in x or 'LOW' in x if x is not None else False)
    df['eqtl'] = df['RSID'].apply(lambda x: 1 if eqtl_dict.get(rs_dict.get(x)) is not None else 0)
    df['sqtl'] = df['RSID'].apply(lambda x: 1 if sqtl_dict.get(rs_dict.get(x)) is not None else 0)
    df['qtl'] = df.apply(lambda row: 1 if row['eqtl'] or row['sqtl'] else 0, axis=1)
    return df

df_rsp_cs_anno = get_enrich_anno_qtl(anno_dict, rs_dict, eqtl_dict, sqtl_dict, '../doc/dat_LDL_rsp_cs.txt')
df_sus_cs_anno = get_enrich_anno_qtl(anno_dict, rs_dict, eqtl_dict, sqtl_dict, '../doc/dat_LDL_sus_cs.txt')

df_rsp_cs_anno.to_csv('../doc/dat_LDL_rsp_cs_anno.txt', sep='\t', header=True, index=False)
df_sus_cs_anno.to_csv('../doc/dat_LDL_sus_cs_anno.txt', sep='\t', header=True, index=False)