source /home/$USER/miniconda3/bin/activate ## on cbsu
source activate anchorwave_new

## paspalum genomes
reffa=../genomes/Pvaginatum_672_v3.0.fa
refgff=../genes/Pvaginatum_672_v3.1.gene.gff3
## set up gene anchors
anchorwave gff2seq -r $reffa -i $refgff -o Pv.CDS.fa
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 $reffa Pv.CDS.fa > Pv.CDS.sam
ref=Pv

## set up file with tab delimited columns:
# name
# genome fasta base
# ploidy (haploid)
## this is panand_sp_ploidy.txt

while read -r altfa name ploidy
do
ploidy=$((2*ploidy))
if [ ! -f ${name}-${ref}.CDS.sam ]
then
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 ../genomes/${altfa}.fasta ${ref}.CDS.fa > ${name}-${ref}.CDS.sam &
fi
if [ ! -f ${name}-${ref}-${ploidy} ]
then
## run anchors
anchorwave proali -t 10 -i $refgff -as ${ref}.CDS.fa -r $reffa -a ${name}-${ref}.CDS.sam -ar ${ref}.CDS.sam -s ../genomes/${altfa}.fasta -n ${name}-${ref}-${ploidy} -R $ploidy -Q 1 -ns & ##-o ${name}-${ref}.maf -f ${name}-${ref}.f.maf &
fi
done < ../panand_sp_ploidy.txt

