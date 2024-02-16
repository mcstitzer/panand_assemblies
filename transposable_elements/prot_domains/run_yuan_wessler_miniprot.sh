#!/bin/bash -login
#SBATCH -D /project/buckler_lab_panand/michelle.stitzer/prot_domains/
#SBATCH -o /project/buckler_lab_panand/michelle.stitzer/slurm-log/teprot-stdout-%j.txt
#SBATCH -e /project/buckler_lab_panand/michelle.stitzer/slurm-log/teprot-stderr-%j.txt
#SBATCH -J teprot
set -e
set -u



conda activate mcs_r



while read -r altfa name ploidy
do
if [ ! -f ${name}.ywminiprot.gff ]
then

miniprot -ut64 -N 100000 --outn=100000 --outs 0.8 --outc 0.8 --gff ../genomes/${altfa}.fasta yuan_wessler_2011.nocommentlines.fa > ${name}.ywminiprot.gff
miniprot -ut64 -N 100000 --outn=100000 --outs 0.8 --outc 0.8 --trans ../genomes/${altfa}.fasta yuan_wessler_2011.nocommentlines.fa > ${name}.ywminiprot.trans
grep -v -P "0\t0\t0" ${name}.ywminiprot.trans > ${name}.ywminiprot.trans.filt
python trans_to_fa.py ${name}.ywminiprot.trans.filt ${name} > ${name}.ywminiprot.trans.fa
fi
done < <( sed -n "${SLURM_ARRAY_TASK_ID}p" ../panand_sp_ploidy.txt )
