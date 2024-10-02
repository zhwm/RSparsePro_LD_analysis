import pandas as pd
import gzip
import argparse

def get_bim_dict(bimdir):
    # get bim
    bim_dict = {}
    bim_dup = []
    with open(bimdir) as bim:
        for line in bim:
            ll = line.strip().split()
            if ll[1] in bim_dict and ll[1] not in bim_dup:
                bim_dup.append(ll[1])
            else:
                bim_dict[ll[1]]=(ll[1], ll[0], ll[3], ll[4], ll[5])
    # remove dups
    for drs in bim_dup:
        del bim_dict[drs]
    return bim_dict

def get_ss_dict(ssdir, bim_dict, RSID, BETA, SE, A1, A2, MAF):
    ss_dict = {}
    with gzip.open(ssdir, 'rt') as ss:
        header = next(ss)
        for line in ss:
            ll = dict(zip(header.strip().split(), line.strip().split()))
            if ll[RSID] in bim_dict:
                if bim_dict.get(ll[RSID])[3]==ll[A1].upper() and bim_dict.get(ll[RSID])[4]==ll[A2].upper():
                    ss_dict[ll[RSID]]=[round(float(ll[BETA]) / float(ll[SE]), 4), round(float(ll[MAF]), 4)]
                elif bim_dict.get(ll[RSID])[4]==ll[A1].upper() and bim_dict.get(ll[RSID])[3]==ll[A2].upper():
                    ss_dict[ll[RSID]]=[round(-float(ll[BETA]) / float(ll[SE]), 4), round(1-float(ll[MAF]), 4)]
    return ss_dict

def parse_args():
    parser = argparse.ArgumentParser(description='format Commands:')
    parser.add_argument('--bim', type=str, default=None, help='path to bim file', required=True)
    parser.add_argument('--rss', type=str, default=None, nargs='+', help='path to raw summary statistics', required=True)
    parser.add_argument('--save', type=str, default=None, nargs='+', help='path to save results', required=True)
    parser.add_argument('--RSID', type=str, default=None, nargs='+', help='RSID', required=True)
    parser.add_argument('--BETA', type=str, default=None, nargs='+', help='BETA', required=True)
    parser.add_argument('--SE', type=str, default=None, nargs='+', help='SE', required=True)
    parser.add_argument('--A1', type=str, default=None, nargs='+', help='effective allele', required=True)
    parser.add_argument('--A2', type=str, default=None, nargs='+', help='reference allele', required=True)
    parser.add_argument('--MAF', type=str, default=None, nargs='+', help='MAF', required=True)
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    bim_dict = get_bim_dict(args.bim)
    ss_lst = [get_ss_dict(rss, bim_dict, rs, beta, se, a1, a2, maf) for rss, rs, beta, se, a1, a2, maf in zip(args.rss, args.RSID, args.BETA, args.SE, args.A1, args.A2, args.MAF)]
    common_snps = set.intersection(*[set(d.keys()) for d in ss_lst])
    df_base = pd.DataFrame.from_dict({k: v for k, v in bim_dict.items() if k in common_snps}, orient='index', columns=['RSID', 'CHR', 'POS', 'A1', 'A2'])
    for i, outf in enumerate(args.save):
        z_df = pd.DataFrame.from_dict({k:v for k,v in ss_lst[i].items() if k in common_snps}, orient='index')
        df_base[['Z', 'MAF']] = z_df
        df_base.to_csv(outf, sep='\t', header=True, index=False)