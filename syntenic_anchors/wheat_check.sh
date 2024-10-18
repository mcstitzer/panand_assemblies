
## run on scinet atlas
cd /project/buckler_lab_panand/michelle.stitzer/wheat_subgenomes/


### download wheat genome
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/triticum_aestivum/dna/Triticum_aestivum.IWGSC.dna.toplevel.fa.gz
gunzip Triticum_aestivum.IWGSC.dna.toplevel.fa.gz

### srun -A buckler_lab_panand -p atlas --ntasks-per-node=20 --time=1-00:00 --pty bash


### get anchors
conda activate anchorwave_new ## atlas scinet

## paspalum genomes
reffa=../Phytozome/PhytozomeV13/Pvaginatum/v3.1/assembly/Pvaginatum_672_v3.0.fa
refgff=../Phytozome/PhytozomeV13/Pvaginatum/v3.1/annotation/Pvaginatum_672_v3.1.gene.gff3
## set up gene anchors
anchorwave gff2seq -r $reffa -i $refgff -o Pv.CDS.fa
minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 $reffa Pv.CDS.fa > Pv.CDS.sam
ref=Pv

## set up file with tab delimited columns:
# name
# genome fasta base
# ploidy (haploid)
## this is panand_sp_ploidy.txt

hap1=Triticum_aestivum.IWGSC.dna.toplevel.fa
ploidy=3

minimap2 -x splice -t 20 -k 12 -a -p 0.4 -N 20 $hap1 ${ref}.CDS.fa > taesti-${ref}.CDS.sam
anchorwave proali -t 20 -i $refgff -as ${ref}.CDS.fa -r $reffa -a taesti-${ref}.CDS.sam -ar ${ref}.CDS.sam -s $hap1 -n taesti-${ref}-${ploidy} -R $ploidy -Q 1 -ns ##-o ${name}-${ref}.maf -f ${name}-${ref}.f.maf &

