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
module load seqkit

module load miniconda3
source activate /project/buckler_lab_panand/michelle.stitzer/conda/envs/anchorwave_new

#module load miniconda3
#conda activate anchorwave_new

## running first with 4 cpu per tree, can fail due to too little memory (~12 gb?)
## rerun again, after
# find -size 0c -delete
## which gets rid of *.withPvCDS.aln.fa (and .aln.fa if it exists)
## then rerun with --mem=48000
## can go up if needed
### quick one chr first??
## find -size 0c -name "Pavag01G*" -delete

## running all at max limit of 2 days on short partition

#ls *0.fa > Pv_falist.txt
SAMPLE_LIST=($(<../syntenic_anchors/sharedSyntenicAnchors.txt))
i=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}]}
gene=$i
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

mkdir -p aligned

if [ ! -f aligned/${gene}.aln.fa ]
then

mafft --thread 12 --genafpair --maxiterate 1000 --adjustdirection ../syntenic_anchors/gene_tree_unalignedfa/${gene}.fa > aligned/${gene}.aln.fa
sed -i 's/()//g' aligned/${gene}.aln.fa
sed -i 's/_R_//g' aligned/${gene}.aln.fa

fi

# ## remove CDS before doing the gene tree
# if [ ! -f ${gene}.aln.fa ]
# then
# sed  '/^Pavag/,+1d' ${gene}.withPvCDS.aln.fa  > ${gene}.aln.fa
# fi

## switch to 3 for bigmem??
mkdir -p trees
if [ ! -f trees/RAxML_bestTree.${gene} ]
then

if [ -f trees/RAxML_info.${gene} ]
then
seqkit rmdup < aligned/${gene}.aln.fa > aligned/${gene}.aln.rmdup.fa
rm -f trees/RAxML*${gene}
raxmlHPC-PTHREADS-AVX2 -T 12 -m GTRGAMMA -p 12345 -x 12345 -# 100 -f a -s aligned/${gene}.aln.rmdup.fa -n ${gene} -w /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/gene_trees/trees/
else

rm -f trees/RAxML*${gene}
raxmlHPC-PTHREADS-AVX2 -T 12 -m GTRGAMMA -p 12345 -x 12345 -# 100 -f a -s aligned/${gene}.aln.fa -n ${gene} -w /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/gene_trees/trees/
fi
fi
