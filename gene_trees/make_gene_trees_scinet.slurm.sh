#!/bin/bash -login
#SBATCH -D /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/gene_trees/
#SBATCH -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/slurm-log/buildtrees-stdout-%j.txt
#SBATCH -e /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/slurm-log/buildtrees-stderr-%j.txt
#SBATCH -J buildtrees
set -e
set -u

## above for atlas


## ceres
##### # SBATCH -D /90daydata/buckler_lab_panand/michelle.stitzer/gene_tree_unalignedfa
##### # SBATCH -o /90daydata/buckler_lab_panand/michelle.stitzer/slurm-log/buildtrees-stdout-%j.txt
##### # SBATCH -e /90daydata/buckler_lab_panand/michelle.stitzer/slurm-log/buildtrees-stderr-%j.txt



#module load raxml
module load mafft
module load samtools

conda activate anchorwave_new

## running first with 4 cpu per tree, can fail due to too little memory (~12 gb?)
## rerun again, after
# find -size 0c -delete
## which gets rid of *.withPvCDS.aln.fa (and .aln.fa if it exists)
## then rerun with --mem=48000
## can go up if needed

## running all at max limit of 2 days on short partition

#ls *0.fa > Pv_falist.txt
SAMPLE_LIST=($(<../syntenic_anchors/sharedSyntenicAnchors.txt))
i=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}
## slurm array limit on ceres
## so for second batch , add 10000
#i=${SAMPLE_LIST[$((SLURM_ARRAY_TASK_ID+10000))]}

# if [[ "$i" != *"aln"* ]]
# then
# gene="${i%%.*}"

## add paspalum gene
## not checking anymore because Sheng-Kai wants the Pasp genomic region in there too (which i added to r script)
# if grep -q "^>$gene" ${gene}.fa; then
# echo "already there"
# else
# samtools faidx ../Pv.CDS.fa ${gene}.1.v3.1 >> ${gene}.fa
# fi


if [ ! -f ${gene}.aln.fa ]
then
mafft --genafpair --maxiterate 1000 --adjustdirection ${gene}.fa > ${gene}.aln.fa
sed -i 's/()//g' ${gene}.aln.fa
sed -i 's/_R_//g' ${gene}.aln.fa

fi

# ## remove CDS before doing the gene tree
# if [ ! -f ${gene}.aln.fa ]
# then
# sed  '/^Pavag/,+1d' ${gene}.withPvCDS.aln.fa  > ${gene}.aln.fa
# fi


if [ ! -f RAxML_bestTree.${gene} ]
then
raxmlHPC-PTHREADS-AVX -T 4 -m GTRGAMMA -p 12345 -x 12345 -# 100 -f a -s ${gene}.aln.fa -n ${gene}

fi
fi
