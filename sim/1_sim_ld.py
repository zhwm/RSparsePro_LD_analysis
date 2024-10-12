import pandas as pd
import dask.array as da
import dask.dataframe as dd
from pandas_plink import read_plink
from scipy.special import softmax
from sklearn.metrics import precision_recall_curve, auc
import numpy as np
import os
import argparse
import logging

def parse_args():
    parser = argparse.ArgumentParser(description='Simulation Commands:')
    parser.add_argument('--bfile', type=str, default=None, help='path to plink file', required=True)
    parser.add_argument('--props', type=float, default=None, nargs='+', help='flip proportion', required=True)
    parser.add_argument('--flips', type=float, default=None, nargs='+', help='flip sign', required=True)
    parser.add_argument('--kc', type=int, default=None, help='number of causal variants', required=True)
    parser.add_argument('--ld', type=str, default=None, help='path to save ld matrix', required=True)
    parser.add_argument('--nite', type=int, default=50, help='number of replication')
    parser.add_argument('--N', type=int, default=500000, help='sample size')
    parser.add_argument('--h2', type=float, default=0.0001, help='per-variant heritability')
    parser.add_argument('--K', type=int, default=10, help='largest number of causal signals')
    parser.add_argument('--maxite', type=int, default=100, help='max number of iterations')
    parser.add_argument('--eps', type=float, default=1e-5, help='convergence criterion')
    parser.add_argument('--ubound', type=int, default=100000, help='upper bound for convergence')
    parser.add_argument('--cthres', type=float, default=0.95, help='attainable coverage threshold for effect groups')
    parser.add_argument('--eincre', type=float, default=1.5, help='adjustment for error parameter')
    parser.add_argument('--minldthres', type=float, default=0.7, help='ld within effect groups')
    parser.add_argument('--maxldthres', type=float, default=0.2, help='ld across effect groups')
    parser.add_argument('--varemax', type=float, default=100.0, help='max error parameter')
    parser.add_argument('--varemin', type=float, default=1e-3, help='min error parameter')
    args = parser.parse_args()
    return args

def get_geno(bdir, lddir):
    bim, fam, bed = read_plink(bdir)
    maf = da.nanmean(bed, axis=1, keepdims=True)
    std = da.nanstd(bed, axis=1, keepdims=True)
    geno = (bed - maf)/std
    geno_n = da.where(da.isnan(geno), 0, geno).compute()
    if os.path.exists(lddir):
        ld = pd.read_csv(lddir, sep=' ', header=None).values
    else:
        ld = np.corrcoef(geno_n)
        np.savetxt(lddir, ld, delimiter=' ', fmt='%.6f')
    return geno_n, ld

def get_z(kc, h2, Ntar, ld, geno_n):
    p, n = geno_n.shape
    cidx = np.random.choice(np.arange(p), size=kc, replace=False)
    while kc != 1 and abs(np.tril(ld[np.ix_(cidx, cidx)], -1)).max() > 0.1:
        cidx = np.random.choice(np.arange(p), size=kc, replace=False)
    btrue = np.zeros(p)
    btrue[cidx] = np.sqrt(h2)
    ynull = np.random.normal(scale=np.sqrt(1 - h2 * kc), size=n) # var(ynull) = 1-kh^2
    zhat = np.dot(ld, btrue*np.sqrt(Ntar)) + np.dot(geno_n, ynull)/np.sqrt(n) # bhat = Rb + Xty/n* sqrt(n/Ntar) 
    return zhat.round(4), btrue

class RSparsePro(object):
    def __init__(self, P, K, R, vare):
        self.p = P
        self.k = K
        self.vare = vare
        if vare != 0:
            self.mat = np.dot(R, np.linalg.inv(np.eye(self.p) + 1/vare * R))
        self.beta_mu = np.zeros([self.p, self.k])
        self.gamma = np.zeros([self.p, self.k])
        self.tilde_b = np.zeros((self.p,))
    
    def infer_q_beta(self, R):
        for k in range(self.k):
            idxall = [x for x in range(self.k)]
            idxall.remove(k)
            beta_all_k = (self.gamma[:, idxall] * self.beta_mu[:, idxall]).sum(axis=1)
            res_beta = self.tilde_b - np.dot(R, beta_all_k)
            self.beta_mu[:, k] = res_beta
            u = 0.5 * self.beta_mu[:, k] ** 2
            self.gamma[:, k] = softmax(u)
    
    def infer_tilde_b(self, bhat):
        if self.vare == 0:
            self.tilde_b = bhat
        else:
            beta_all = (self.gamma * self.beta_mu).sum(axis=1)
            self.tilde_b = np.dot(self.mat, (1/self.vare * bhat + beta_all))
    
    def train(self, bhat, R, maxite, eps, ubound):
        for ite in range(maxite):
            old_gamma = self.gamma.copy()
            old_beta = self.beta_mu.copy()
            old_tilde = self.tilde_b.copy()
            self.infer_tilde_b(bhat)
            self.infer_q_beta(R)
            diff_gamma = np.linalg.norm(self.gamma-old_gamma)
            diff_beta = np.linalg.norm(self.beta_mu - old_beta)
            diff_b = np.linalg.norm(self.tilde_b - old_tilde)
            all_diff = diff_gamma + diff_beta + diff_b
            print('Iteration-->{} . Diff_b: {:.1f} . Diff_s: {:.1f} . Diff_mu: {:.1f} . ALL: {:.1f}'.format(ite, diff_b, diff_gamma, diff_beta, all_diff))
            if all_diff < eps:
                logging.info("The RSparsePro algorithm has converged.")
                converged = True
                break
            if ite == (maxite - 1) or abs(all_diff) > ubound:
                logging.info("The RSparsePro algorithm didn't converge.")
                converged = False
                break
        return converged
    
    def get_PIP(self):
        return np.max((self.gamma), axis=1).round(4)
    
    def get_effect(self, cthres):
        vidx = np.argsort(-self.gamma, axis=1)
        matidx = np.argsort(-self.gamma, axis=0)
        mat_eff = np.zeros((self.p, self.k))
        for p in range(self.p):
            mat_eff[p, vidx[p, 0]] = self.gamma[p, vidx[p, 0]]
        mat_eff[mat_eff < 1/(self.p+1)] = 0
        csum = mat_eff.sum(axis=0).round(2)
        logging.info("Attainable coverage for effect groups: {}".format(csum))
        eff = {}
        eff_gamma = {}
        eff_mu = {}
        for k in range(self.k):
            if csum[k] >= cthres:
                p = 0
                while np.sum(mat_eff[matidx[0:p, k], k]) < cthres * csum[k]:
                    p = p + 1
                cidx = matidx[0:p, k].tolist()
                eff[k] = cidx
                eff_gamma[k] = mat_eff[cidx, k].round(4)
                eff_mu[k] = self.beta_mu[cidx, k].round(4)
        return eff, eff_gamma, eff_mu
    
    def get_zdiff(self, bhat):
        zdiff = bhat - self.tilde_b
        return zdiff.round(4)

def get_eff_maxld(eff, ld):
    idx = [i[0] for i in eff.values()]
    if len(eff)>1:
        maxld = np.abs(np.tril(ld[np.ix_(idx,idx)],-1)).max()
    else:
        maxld = 0.0
    return maxld

def get_eff_minld(eff, ld):
    if len(eff)==0:
        minld = 1.0
    else:
        minld = min([abs(ld[np.ix_(v, v)]).min() for _,v in eff.items()])
    return minld

def get_ordered(eff_mu):
    if len(eff_mu)>1:
        ordered = (list(eff_mu.keys())[-1] == len(eff_mu)-1)
    else:
        ordered = True
    return ordered

def adaptive_train(zscore, ld, K, maxite, eps, ubound, cthres, minldthres, maxldthres, eincre, varemax, varemin):
    vare = 0
    mc = False
    while (not mc) or (not get_ordered(eff_mu)) or (minld < minldthres) or (maxld > maxldthres):
        model = RSparsePro(len(zscore), K, ld, vare)
        mc = model.train(zscore, ld, maxite, eps, ubound)
        eff, eff_gamma, eff_mu = model.get_effect(cthres)
        maxld = get_eff_maxld(eff, ld)
        minld = get_eff_minld(eff, ld)
        print("Max ld across effect groups: {}.".format(maxld))
        print("Min ld within effect groups: {}.".format(minld))
        print("vare = {}".format(round(vare,4)))
        if vare > varemax or (len(eff)<2 and get_ordered(eff_mu)):
            model = RSparsePro(len(zscore), 1, ld, 0)
            mc = model.train(zscore, ld, maxite, eps, ubound)
            eff, eff_gamma, eff_mu = model.get_effect(cthres)
            break
        elif vare ==0:
            vare = varemin
        else:
            vare *= eincre
    zdiff = model.get_zdiff(zscore)
    PIP = model.get_PIP()
    return eff, eff_gamma, eff_mu, PIP, zdiff

def get_auprc(pip, gt):
    precision, recall, _ = precision_recall_curve(gt, pip)
    auprc = auc(recall, precision)
    return auprc

#def get_pwr_cov(eff, gt, ld):
#    cstrue = set(np.where(gt!=0)[0])
#    if len(eff)==0:
#        return 0, 1, 1, 0, 0
#    truep = cstrue.intersection(set(var for subset in eff.values() for var in subset))
#    power = len(truep) / len(cstrue)
#    coverage = len(truep) / len(eff)
#    minld = np.min([abs(ld[np.ix_(var, list(cstrue))]).max() for var in eff.values()])
#    nofcs = len(eff)
#    msz = np.median([len(subset) for subset in eff.values()]) 
#    return power, coverage, minld, nofcs, msz

def get_cs_vec(eff, p):
    csvec = np.zeros(p)
    for e,val in eff.items():
        csvec[val] = e+1
    return csvec

if __name__ == "__main__":
    args = parse_args()
    np.random.seed(42)
    geno_n, ld = get_geno(args.bfile, args.ld)
    out_matrix_btrue = np.zeros((geno_n.shape[0], args.nite))
    out_matrix_ztrue = np.zeros((geno_n.shape[0], args.nite))
    out_matrix_z = np.zeros((geno_n.shape[0], len(args.props) * args.nite * len(args.flips)))
    out_matrix_ztrue_pip = np.zeros((geno_n.shape[0], args.nite))
    out_matrix_z_pip = np.zeros((geno_n.shape[0], len(args.props) * args.nite * len(args.flips)))
    out_matrix_ztrue_cs = np.zeros((geno_n.shape[0], args.nite))
    out_matrix_z_cs = np.zeros((geno_n.shape[0], len(args.props) * args.nite * len(args.flips)))
    out_matrix_ztrue_zdiff = np.zeros((geno_n.shape[0], args.nite))
    out_matrix_z_zdiff = np.zeros((geno_n.shape[0], len(args.props) * args.nite * len(args.flips)))
    for ite in range(args.nite):
        ztrue, btrue = get_z(args.kc, args.h2, args.N, ld, geno_n)
        out_matrix_btrue[:, ite] = btrue
        out_matrix_ztrue[:, ite] = ztrue
        eff, eff_gamma, eff_mu, PIP, zdiff = adaptive_train(ztrue, ld, args.K, args.maxite, args.eps, args.ubound, args.cthres, args.minldthres, args.maxldthres, args.eincre, args.varemax, args.varemin)
        out_matrix_ztrue_pip[:, ite] = PIP
        out_matrix_ztrue_zdiff[:, ite] = zdiff
        out_matrix_ztrue_cs[:, ite] = get_cs_vec(eff, geno_n.shape[0])
        auprc = get_auprc(PIP, btrue!=0)
        for pidx, prop in enumerate(args.props):
            km = int(np.ceil(prop * len(ztrue)))
            midx = np.random.choice(np.arange(len(ztrue)), size=km, replace=False)
            for fidx, flip in enumerate(args.flips):
                start_idx = ite * len(args.props) * len(args.flips) + pidx * len(args.flips) + fidx
                zerr = ztrue.copy()
                zerr[midx] *= flip
                out_matrix_z[:, start_idx] = zerr
                eff, eff_gamma, eff_mu, PIP, zdiff = adaptive_train(zerr, ld, args.K, args.maxite, args.eps, args.ubound, args.cthres, args.minldthres, args.maxldthres, args.eincre, args.varemax, args.varemin)
                out_matrix_z_pip[:, start_idx] = PIP
                out_matrix_z_zdiff[:, start_idx] = zdiff
                out_matrix_z_cs[:, start_idx] = get_cs_vec(eff, geno_n.shape[0])
                auprc = get_auprc(PIP, btrue!=0)
    out_matrix_z_df = pd.DataFrame(out_matrix_z, columns=[f"f{flip}_p{prop}_ite{ite}" for ite in range(args.nite) for prop in args.props for flip in args.flips])
    out_matrix_ztrue_df = pd.DataFrame(out_matrix_ztrue, columns=[f"ite{ite}" for ite in range(args.nite)])
    out_matrix_btrue_df = pd.DataFrame(out_matrix_btrue, columns=[f"ite{ite}" for ite in range(args.nite)])    
    out_matrix_btrue_df.to_csv(f"{args.bfile}-{args.kc}-{args.h2}-{args.N}-btrue.txt", sep='\t', index=False)
    out_matrix_ztrue_df.to_csv(f"{args.bfile}-{args.kc}-{args.h2}-{args.N}-ztrue.txt", sep='\t', index=False)
    out_matrix_z_df.to_csv(f"{args.bfile}-{args.kc}-{args.h2}-{args.N}-z.txt", sep='\t', index=False)
    out_matrix_ztrue_pip_df = pd.DataFrame(out_matrix_ztrue_pip, columns=[f"ite{ite}" for ite in range(args.nite)])
    out_matrix_z_pip_df = pd.DataFrame(out_matrix_z_pip, columns=[f"f{flip}_p{prop}_ite{ite}" for ite in range(args.nite) for prop in args.props for flip in args.flips])
    out_matrix_ztrue_cs_df = pd.DataFrame(out_matrix_ztrue_cs, columns=[f"ite{ite}" for ite in range(args.nite)])
    out_matrix_z_cs_df = pd.DataFrame(out_matrix_z_cs, columns=[f"f{flip}_p{prop}_ite{ite}" for ite in range(args.nite) for prop in args.props for flip in args.flips])
    out_matrix_ztrue_zdiff_df = pd.DataFrame(out_matrix_ztrue_zdiff, columns=[f"ite{ite}" for ite in range(args.nite)])
    out_matrix_z_zdiff_df = pd.DataFrame(out_matrix_z_zdiff, columns=[f"f{flip}_p{prop}_ite{ite}" for ite in range(args.nite) for prop in args.props for flip in args.flips])
    out_matrix_ztrue_pip_df.to_csv(f"{args.bfile}-{args.kc}-{args.h2}-{args.N}-ztrue.pip.txt", sep='\t', index=False)
    out_matrix_z_pip_df.to_csv(f"{args.bfile}-{args.kc}-{args.h2}-{args.N}-z.pip.txt", sep='\t', index=False)
    out_matrix_ztrue_cs_df.to_csv(f"{args.bfile}-{args.kc}-{args.h2}-{args.N}-ztrue.cs.txt", sep='\t', index=False)
    out_matrix_z_cs_df.to_csv(f"{args.bfile}-{args.kc}-{args.h2}-{args.N}-z.cs.txt", sep='\t', index=False)
    out_matrix_ztrue_zdiff_df.to_csv(f"{args.bfile}-{args.kc}-{args.h2}-{args.N}-ztrue.zdiff.txt", sep='\t', index=False)
    out_matrix_z_zdiff_df.to_csv(f"{args.bfile}-{args.kc}-{args.h2}-{args.N}-z.zdiff.txt", sep='\t', index=False)
