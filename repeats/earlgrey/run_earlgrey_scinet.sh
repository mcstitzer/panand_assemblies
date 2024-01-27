#!/bin/bash -login
#SBATCH -D /project/buckler_lab_panand/michelle.stitzer/earlgrey/
#SBATCH -o /project/buckler_lab_panand/michelle.stitzer/slurm-log/earlgrey-stdout-%j.txt
#SBATCH -e /project/buckler_lab_panand/michelle.stitzer/slurm-log/earlgrey-stderr-%j.txt
#SBATCH -J earlgrey
set -e
set -u



conda activate earlgrey



while read -r altfa name ploidy
do
if [ ! -f earlGreyOutputs/${name}_EarlGrey/${name}_summaryFiles/${name}.highLevelCount.txt ]
then

earlGrey -g ../genomes/${altfa}.fasta -s ${name} -o ./earlGreyOutputs -t 128

fi
done < <( set -n "${SLURM_ARRAY_TASK_ID}p" ../panand_sp_ploidy.txt )

