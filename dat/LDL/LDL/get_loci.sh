## get loci
### AFR
zcat raw/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz | cut -f 6 | sed 1d | sort -gr | head -1 > loci/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.nmax
python get_lead_LDL.py --raw raw/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz --save loci/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.leads.txt --mafcol POOLED_ALT_AF --pcol pvalue --ncol N --chr CHROM --pos POS_b37 --a1 REF --a2 ALT --nmin 85162
while read a b c; do plink --bfile ~/projects/def-glettre/wmzh22/1KG/$a --keep ~/projects/def-glettre/wmzh22/1KG/AFR.fam --chr $a --from-bp $b --to-bp $c --maf 0.05 --make-bed --out bed/AFR_$a\_$b\_$c; done < loci/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.leads.txt.gtf 

### EUR
zcat raw/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz | cut -f 6 | sed 1d | sort -gr | head -1 > loci/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.nmax
python get_lead_LDL.py --raw raw/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz --save loci/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.leads.txt --mafcol POOLED_ALT_AF --pcol pvalue --ncol N --chr CHROM --pos POS_b37 --a1 REF --a2 ALT --nmin 1108160
while read a b c; do plink --bfile ~/projects/def-glettre/wmzh22/1KG/$a --keep ~/projects/def-glettre/wmzh22/1KG/EUR.fam --chr $a --from-bp $b --to-bp $c --maf 0.05 --make-bed --out bed/EUR_$a\_$b\_$c; done < loci/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.leads.txt.gtf

### EAS
zcat raw/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz | cut -f 6 | sed 1d | sort -gr | head -1 > loci/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz.nmax
python get_lead_LDL.py --raw raw/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz --save loci/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz.leads.txt --mafcol POOLED_ALT_AF --pcol pvalue --ncol N --chr CHROM --pos POS_b37 --a1 REF --a2 ALT --nmin 74328
while read a b c; do plink --bfile ~/projects/def-glettre/wmzh22/1KG/$a --keep ~/projects/def-glettre/wmzh22/1KG/EAS.fam --chr $a --from-bp $b --to-bp $c --maf 0.05 --make-bed --out bed/EAS_$a\_$b\_$c; done < loci/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz.leads.txt.gtf 

### SAS
zcat raw/LDL_INV_SAS_HRC_1KGP3_others_ALL.meta.singlevar.results.gz | cut -f 6 | sed 1d | sort -gr | head -1 > loci/LDL_INV_SAS_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.nmax
python get_lead_LDL.py --raw raw/LDL_INV_SAS_HRC_1KGP3_others_ALL.meta.singlevar.results.gz --save loci/LDL_INV_SAS_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.leads.txt --mafcol POOLED_ALT_AF --pcol pvalue --ncol N --chr CHROM --pos POS_b37 --a1 REF --a2 ALT --nmin 36426
while read a b c; do plink --bfile ~/projects/def-glettre/wmzh22/1KG/$a --keep ~/projects/def-glettre/wmzh22/1KG/SAS.fam --chr $a --from-bp $b --to-bp $c --maf 0.05 --make-bed --out bed/SAS_$a\_$b\_$c; done < loci/LDL_INV_SAS_HRC_1KGP3_others_ALL.meta.singlevar.results.gz.leads.txt.gtf

### AMR
zcat raw/LDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz | cut -f 6 | sed 1d | sort -gr | head -1 > loci/LDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz.nmax
python get_lead_LDL.py --raw raw/LDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz --save loci/LDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz.leads.txt --mafcol POOLED_ALT_AF --pcol pvalue --ncol N --chr CHROM --pos POS_b37 --a1 REF --a2 ALT --nmin 41436
while read a b c; do plink --bfile ~/projects/def-glettre/wmzh22/1KG/$a --keep ~/projects/def-glettre/wmzh22/1KG/AMR.fam --chr $a --from-bp $b --to-bp $c --maf 0.05 --make-bed --out bed/AMR_$a\_$b\_$c; done < loci/LDL_INV_HIS_1KGP3_ALL.meta.singlevar.results.gz.leads.txt.gtf
