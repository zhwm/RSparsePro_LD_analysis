import glob
import pickle
import gzip

elist = glob.glob('../dat/GTEx/GTEx_Analysis_v8_eQTL/*.signif*')
slist = glob.glob('../dat/GTEx/GTEx_Analysis_v8_sQTL/*signif*')

eqtl_dict = {}
for e in elist:
    print(e)
    with gzip.open(e, 'rt') as qtl:
        header = next(qtl)
        for line in qtl:
            ll = dict(zip(header.strip().split(), line.strip().split()))
            var_id = '_'.join(ll['variant_id'].split('_')[0:2])
            gene_id = ll['gene_id'].split('.')[0]
            tissue_id = e.split('v8')[1]
            if var_id not in eqtl_dict:
                eqtl_dict[var_id] = {}
            if gene_id not in eqtl_dict[var_id]:
                eqtl_dict[var_id][gene_id] = set()
            eqtl_dict[var_id][gene_id].add(tissue_id)

with open('../dat/GTEx/eqtl_dict.pkl', 'wb') as file:
    pickle.dump(eqtl_dict, file)

sqtl_dict = {}
for e in slist:
    print(e)
    with gzip.open(e, 'rt') as qtl:
        header = next(qtl)
        for line in qtl:
            ll = dict(zip(header.strip().split(), line.strip().split()))
            var_id = '_'.join(ll['variant_id'].split('_')[0:2])
            gene_id = ll['phenotype_id'].split(':')[4].split('.')[0]
            tissue_id = e.split('v8')[1]
            if var_id not in sqtl_dict:
                sqtl_dict[var_id] = {}
            if gene_id not in sqtl_dict[var_id]:
                sqtl_dict[var_id][gene_id] = set()
            sqtl_dict[var_id][gene_id].add(tissue_id)

with open('../dat/GTEx/sqtl_dict.pkl', 'wb') as file:
    pickle.dump(sqtl_dict, file)