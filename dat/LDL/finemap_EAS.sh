#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 10G            # memory pool per process
#SBATCH -o finemap_EAS.out    # STDOUT
#SBATCH -t 06:00:00            # time (D-HH:MM)

module load StdEnv/2020 plink/1.9b_6.21-x86_64 r/4.0.0
source ~/py3/bin/activate

while read a b c; do
  python format_ss.py --rss raw/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz --bim bed/EAS_$a\_$b\_$c\.bim --save formatted/EAS_$a\_$b\_$c\.txt --RSID rsID --BETA EFFECT_SIZE --SE SE --A1 ALT --A2 REF --MAF POOLED_ALT_AF;
  plink --bfile bed/EAS_$a\_$b\_$c --extract formatted/EAS_$a\_$b\_$c\.txt --matrix --r --out bed/EAS_$a\_$b\_$c;
  Rscript ../../src/susie.R formatted/EAS_$a\_$b\_$c\.txt bed/EAS_$a\_$b\_$c\.ld;
  python ../../src/rsparsepro_ld.py --z formatted/EAS_$a\_$b\_$c\.txt --ld bed/EAS_$a\_$b\_$c\.ld --save formatted/EAS_$a\_$b\_$c\.txt;
done < loci/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz.leads.txt.gtf
