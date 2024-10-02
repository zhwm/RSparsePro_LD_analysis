#!/bin/bash
#
#SBATCH --ntasks 1            # number of tasks
#SBATCH --mem 100G            # memory pool per process
#SBATCH -o run_paired_lst_all.out    # STDOUT
#SBATCH -t 24:00:00            # time (D-HH:MM)

module load StdEnv/2020 plink/1.9b_6.21-x86_64 r/4.0.0
source ~/py3/bin/activate

while read a b c d e f g h; do
  plink --bfile ~/scratch/UKB/$a --keep UKB_EUR_5000.fam --chr $a --from-bp $b --to-bp $c --make-bed --out $d;
  python ../../src/format_ss.py --bim $d\.bim --rss $e $f --save $g $h --RSID rsids rsid --BETA Beta Effect --SE SE StdErr --A1 effectAllele Allele1 --A2 otherAllele Allele2 --MAF ImpMAF Freq1;
  plink --bfile $d --extract $g --matrix --r --out $d;
  python ../../src/rsparsepro_ld.py --z $g --ld $d\.ld --save $g;
  python ../../src/rsparsepro_ld.py --z $h --ld $d\.ld --save $h;
  Rscript ../../src/susie.R $g $d\.ld;
  Rscript ../../src/susie.R $h $d\.ld
done < paired.lst.all