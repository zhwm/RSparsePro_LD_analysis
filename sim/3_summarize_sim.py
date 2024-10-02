import pandas as pd
from sklearn.metrics import precision_recall_curve, auc
import numpy as np

def get_auprc(pip, gt):
    if np.isnan(pip).any():
        auprc = 0.0
    else:
        precision, recall, _ = precision_recall_curve(gt, pip)
        auprc = auc(recall, precision)
    return auprc

def get_ncs(csvec):
    ncs = csvec.max()
    return ncs

def get_ntcs(csvec, gt):
    if (csvec==0).all():
        ntcs = 0
    else:
        ntcs = csvec[(gt!=0) & (csvec!=0)].nunique()
    return ntcs

def get_csz(csvec):
    lst = list(csvec[csvec!=0].value_counts())
    if len(lst) == 0:
        csz = np.nan
    else:
        csz = np.median(lst)
    return csz

df_match_lst = []
df_mismatch_lst = []
for i in range(1, 11):
    for kc in [3, 5]:
        loci = 'Locus{}'.format(i)
        btrue = pd.read_csv('{}-{}-0.0001-500000-btrue.txt'.format(loci, kc), sep='\t')
        z_rsp_cs = pd.read_csv('{}-{}-0.0001-500000-z.cs.txt'.format(loci, kc), sep='\t')
        z_sus_cs = pd.read_csv('{}-{}-0.0001-500000-z.txt.susie.z.cs.txt'.format(loci, kc), sep='\t')
        z_rsp_pip = pd.read_csv('{}-{}-0.0001-500000-z.pip.txt'.format(loci, kc), sep='\t')
        z_sus_pip = pd.read_csv('{}-{}-0.0001-500000-z.txt.susie.z.pip.txt'.format(loci, kc), sep='\t')
        ztrue_rsp_cs = pd.read_csv('{}-{}-0.0001-500000-ztrue.cs.txt'.format(loci, kc), sep='\t')
        ztrue_sus_cs = pd.read_csv('{}-{}-0.0001-500000-z.txt.susie.ztrue.cs.txt'.format(loci, kc), sep='\t')
        ztrue_rsp_pip = pd.read_csv('{}-{}-0.0001-500000-ztrue.pip.txt'.format(loci, kc), sep='\t')
        ztrue_sus_pip = pd.read_csv('{}-{}-0.0001-500000-z.txt.susie.ztrue.pip.txt'.format(loci, kc), sep='\t')
        z_sus_summary = pd.read_csv('{}-{}-0.0001-500000-z.txt.susie.summary.txt'.format(loci, kc), sep='\t')
        df_match_summary = z_sus_summary.loc[0:49].copy()
        df_match_summary['sus_auprc'] = [get_auprc(ztrue_sus_pip['ite{}'.format(i)], btrue['ite{}'.format(i)]!=0) for i in range(50)]
        df_match_summary['rsp_auprc'] = [get_auprc(ztrue_rsp_pip['ite{}'.format(i)], btrue['ite{}'.format(i)]!=0) for i in range(50)]
        df_match_summary['sus_ncs'] = [get_ncs(ztrue_sus_cs['ite{}'.format(i)]) for i in range(50)]
        df_match_summary['rsp_ncs'] = [get_ncs(ztrue_rsp_cs['ite{}'.format(i)]) for i in range(50)]
        df_match_summary['sus_ntcs'] = [get_ntcs(ztrue_sus_cs['ite{}'.format(i)], btrue['ite{}'.format(i)]) for i in range(50)]
        df_match_summary['rsp_ntcs'] = [get_ntcs(ztrue_rsp_cs['ite{}'.format(i)], btrue['ite{}'.format(i)]) for i in range(50)]
        df_match_summary['sus_csz'] = [get_csz(ztrue_sus_cs['ite{}'.format(i)]) for i in range(50)]
        df_match_summary['rsp_csz'] = [get_csz(ztrue_rsp_cs['ite{}'.format(i)]) for i in range(50)]
        df_match_summary['loci'] = loci
        df_match_summary['kc'] = kc
        df_match_summary['ite'] = ['ite{}'.format(i) for i in range(50)]
        df_match_lst.append(df_match_summary)
        df_mismatch_summary = z_sus_summary.loc[50:].copy()
        df_mismatch_summary['sus_auprc'] = [get_auprc(z_sus_pip[f'f{flip}_p{prop}_ite{ite}'], btrue[f'ite{ite}'] != 0) for ite in range(50) for prop in [0.001, 0.01, 0.1] for flip in [-1.0, -0.5, 0.0, 0.5]]
        df_mismatch_summary['rsp_auprc'] = [get_auprc(z_rsp_pip[f'f{flip}_p{prop}_ite{ite}'], btrue[f'ite{ite}'] != 0) for ite in range(50) for prop in [0.001, 0.01, 0.1] for flip in [-1.0, -0.5, 0.0, 0.5]]
        df_mismatch_summary['sus_ncs'] = [get_ncs(z_sus_cs[f'f{flip}_p{prop}_ite{ite}']) for ite in range(50) for prop in [0.001, 0.01, 0.1] for flip in [-1.0, -0.5, 0.0, 0.5]]
        df_mismatch_summary['rsp_ncs'] = [get_ncs(z_rsp_cs[f'f{flip}_p{prop}_ite{ite}']) for ite in range(50) for prop in [0.001, 0.01, 0.1] for flip in [-1.0, -0.5, 0.0, 0.5]]
        df_mismatch_summary['sus_ntcs'] = [get_ntcs(z_sus_cs[f'f{flip}_p{prop}_ite{ite}'], btrue[f'ite{ite}']) for ite in range(50) for prop in [0.001, 0.01, 0.1] for flip in [-1.0, -0.5, 0.0, 0.5]]
        df_mismatch_summary['rsp_ntcs'] = [get_ntcs(z_rsp_cs[f'f{flip}_p{prop}_ite{ite}'], btrue[f'ite{ite}']) for ite in range(50) for prop in [0.001, 0.01, 0.1] for flip in [-1.0, -0.5, 0.0, 0.5]]
        df_mismatch_summary['sus_csz'] = [get_csz(z_sus_cs[f'f{flip}_p{prop}_ite{ite}']) for ite in range(50) for prop in [0.001, 0.01, 0.1] for flip in [-1.0, -0.5, 0.0, 0.5]]
        df_mismatch_summary['rsp_csz'] = [get_csz(z_rsp_cs[f'f{flip}_p{prop}_ite{ite}']) for ite in range(50) for prop in [0.001, 0.01, 0.1] for flip in [-1.0, -0.5, 0.0, 0.5]]
        df_mismatch_summary['loci'] = loci
        df_mismatch_summary['kc'] = kc
        df_mismatch_summary['flip'] = [i.split('_')[0] for i in z_rsp_cs.columns]
        df_mismatch_summary['prop'] = [i.split('_')[1] for i in z_rsp_cs.columns]
        df_mismatch_summary['ite'] = [i.split('_')[2] for i in z_rsp_cs.columns]
        df_mismatch_lst.append(df_mismatch_summary)

df_match = pd.concat(df_match_lst, ignore_index=True)
df_mismatch = pd.concat(df_mismatch_lst, ignore_index=True)

summary_dict = []
for k in [3, 5]:
    for f in [1, 0.5, 0.0, -0.5, -1.0]:
        if f==1:
            p = 0.0
            df = df_match.loc[df_match['kc']==k]
            coverage_sus = df['sus_ntcs'].sum()/df['sus_ncs'].sum()
            coverage_rsp = df['rsp_ntcs'].sum()/df['rsp_ncs'].sum()
            power_sus = df['sus_ntcs'].sum()/k/df.shape[0]
            power_rsp = df['rsp_ntcs'].sum()/k/df.shape[0]
            sus_convrg = np.mean(df['convrg'])
            summary_dict.append({'Number_of_Causal_Variants': k, 
                                'Deviation_Factor': f, 
                                'Mismatch_Proportion': p, 
                                'RSparsePro_Power': power_rsp, 
                                'RSparsePro_Coverage': coverage_rsp,
                                'SuSiE_Convergence_Rate': sus_convrg,
                                'SuSiE_Power': power_sus, 
                                'SuSiE_Coverage': coverage_sus})
        else:
            for p in ['0.001', '0.01', '0.1']:
                df = df_mismatch.loc[(df_mismatch['kc']==k) & (df_mismatch['prop']=='p{}'.format(p)) & (df_mismatch['flip']=='f{}'.format(f))]
                coverage_sus = df['sus_ntcs'].sum()/df['sus_ncs'].sum()
                coverage_rsp = df['rsp_ntcs'].sum()/df['rsp_ncs'].sum()
                power_sus = df['sus_ntcs'].sum()/k/df.shape[0]
                power_rsp = df['rsp_ntcs'].sum()/k/df.shape[0]
                msize_sus = np.nanmedian(df['sus_csz'])
                msize_rsp = np.nanmedian(df['rsp_csz'])
                sus_convrg = np.mean(df['convrg'])
                summary_dict.append({'Number_of_Causal_Variants': k, 
                                    'Deviation_Factor': f, 
                                    'Mismatch_Proportion': p, 
                                    'RSparsePro_Power': power_rsp, 
                                    'RSparsePro_Coverage': coverage_rsp,
                                    'SuSiE_Convergence_Rate': sus_convrg,
                                    'SuSiE_Power': power_sus, 
                                    'SuSiE_Coverage': coverage_sus})

summary_categorized_dict = []
for k in [3, 5]:
    for f in [1, 0.5, 0.0, -0.5, -1.0]:
        if f==1:
            p = 0.0
            for c in ['SuSiE inference converged', 'SuSiE inference did not converge', 'SuSiE prior estimation failed']:
                if c=='SuSiE inference converged':
                    df = df_match.loc[(df_match['kc']==k) & (df_match['exist']==1) & (df_match['convrg']==1)]
                elif c=='SuSiE inference did not converge':
                    df = df_match.loc[(df_match['kc']==k) & (df_match['exist']==1) & (df_match['convrg']==0)]
                else:
                    df = df_match.loc[(df_match['kc']==k) & (df_match['exist']==0) & (df_match['convrg']==0)]
                for m in ['RSparsePro', 'SuSiE']:
                    nc = df.shape[0]
                    if nc == 0:
                        power, coverage = np.nan, np.nan
                    else:
                        if m=='SuSiE':
                            power = df['sus_ntcs'].sum()/k/nc
                            coverage = df['sus_ntcs'].sum()/df['sus_ncs'].sum()
                        else:
                            power = df['rsp_ntcs'].sum()/k/nc
                            coverage = df['rsp_ntcs'].sum()/df['rsp_ncs'].sum()
                    summary_categorized_dict.append({'Number_of_Causal_Variants': k, 
                                'Deviation_Factor': f, 
                                'Mismatch_Proportion': p, 
                                'Category': c,
                                'Method': m,
                                'Power': power, 
                                'Coverage': coverage,
                                'Cases': nc})
        else:
            for p in ['0.001', '0.01', '0.1']:
                for c in ['SuSiE inference converged', 'SuSiE inference did not converge', 'SuSiE prior estimation failed']:
                    if c=='SuSiE inference converged':
                        df = df_mismatch.loc[(df_mismatch['kc']==k) & (df_mismatch['exist']==1) & (df_mismatch['convrg']==1) & (df_mismatch['prop']=='p{}'.format(p)) & (df_mismatch['flip']=='f{}'.format(f))]
                    elif c=='SuSiE inference did not converge':
                        df = df_mismatch.loc[(df_mismatch['kc']==k) & (df_mismatch['exist']==1) & (df_mismatch['convrg']==0) & (df_mismatch['prop']=='p{}'.format(p)) & (df_mismatch['flip']=='f{}'.format(f))]
                    else:
                        df = df_mismatch.loc[(df_mismatch['kc']==k) & (df_mismatch['exist']==0) & (df_mismatch['convrg']==0) & (df_mismatch['prop']=='p{}'.format(p)) & (df_mismatch['flip']=='f{}'.format(f))]
                    for m in ['RSparsePro', 'SuSiE']:
                        nc = df.shape[0]
                        if nc == 0:
                            power, coverage = np.nan, np.nan
                        else:
                            if m=='SuSiE':
                                power = df['sus_ntcs'].sum()/k/nc
                                coverage = df['sus_ntcs'].sum()/df['sus_ncs'].sum()
                            else:
                                power = df['rsp_ntcs'].sum()/k/nc
                                coverage = df['rsp_ntcs'].sum()/df['rsp_ncs'].sum()
                        summary_categorized_dict.append({'Number_of_Causal_Variants': k, 
                                'Deviation_Factor': f, 
                                'Mismatch_Proportion': p, 
                                'Category': c,
                                'Method': m,
                                'Power': power, 
                                'Coverage': coverage,
                                'Cases': nc})

df_mismatch.to_csv('../doc/sim_mismatch.txt', sep='\t', header=True, index=False)
df_match.to_csv('../doc/sim_match.txt', sep='\t', header=True, index=False)
pd.DataFrame(summary_categorized_dict).to_csv('../doc/sim_pwr_cov_categorized.txt', sep='\t', header=True, index=False)
pd.DataFrame(summary_dict).to_csv('../doc/sim_pwr_cov.txt', sep='\t', header=True, index=False)