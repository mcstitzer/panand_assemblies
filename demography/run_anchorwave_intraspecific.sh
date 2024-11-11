
conda activate anchorwave_recent
## want to make sure i'm using the maf fix!!!!
## ugh bioconda isn't working and i can't figure out pixi
./anchorwave/anchorwave/anchorwave # v 1.2.6

### for hconto
## scaffolds homologous to Chr01 from Paspalum
# scaf_5   Chr01  64211423
# scaf_7   Chr01  61978727
# scaf_6   Chr01  60417967
# scaf_11  Chr01  58807442




scp cbsublfs1:/data2/PanAnd2/assemblies/andropogoneae/private/Corteva_pop/Hc-AUB53-1-Draft-PanAnd-1.0.fasta.gz .

refscaf=scaf_5

queryscafs=[scaf_7, scaf_6, scaf_11]

samtools faidx Hc-AUB53-1-Draft-PanAnd-1.0.fasta $refscaf > ${genome}-${refscaf}.fa

for queryscaf in queryscafs
samtools faidx Hc-AUB53-1-Draft-PanAnd-1.0.fasta $queryscaf > ${genome}-${queryscaf}.fa

# Define the reference scaffold and query scaffolds array
refscaf="scaf_5"
queryscafs=("scaf_7" "scaf_6" "scaf_11")
genome=hconto
genomefa=Hc-AUB53-1-Draft-PanAnd-1.0.fasta

# Run samtools faidx for the reference scaffold
samtools faidx $genomefa $refscaf > "${genome}-${refscaf}.fa"

# Loop through each query scaffold and run samtools faidx
for queryscaf in "${queryscafs[@]}"
do
    samtools faidx $genomefa $queryscaf > "${genome}-${queryscaf}.fa"
done



## paspalum genomes
reffa=${genome}-${refscaf}.fa

## need a gff, so run miniprot!
singularity run -B $PWD --pwd $PWD /programs/EMBOSS/emboss.sif transeq -sequence ~/transfer/Pvaginatum_672_v3.1.cds_primaryTranscriptOnly.fa  -outseq Pvaginatum_primaryTranscriptOnly.aa.fa

/programs/miniprot-0.13/miniprot -t 20 --gff ${genome}-${refscaf}.fa Pvaginatum_primaryTranscriptOnly.aa.fa > $genome-$refscaf.gff

scp cbsublfs1:/data4/users/zrm22/HelixerRuns/annotations/Hc-AUB53_1-DRAFT-PanAnd-1.0_helixer.gff .
grep "scaf_5	" Hc-AUB53_1-DRAFT-PanAnd-1.0_helixer.gff > Hc-AUB53_1-DRAFT-PanAnd-1.0_helixer.scaf_5.gff


refgff=Hc-AUB53_1-DRAFT-PanAnd-1.0_helixer.scaf_5.gff
reffa=${genome}-${refscaf}.fa
## set up gene anchors
./anchorwave/anchorwave/anchorwave gff2seq -r $reffa -i $refgff -o $genome-$refscaf.CDS.fa
/programs/minimap2-2.28/minimap2 -x splice -t 10 -k 12 -a -p 0.4 -N 20 $reffa $genome-$refscaf.CDS.fa > $genome-$refscaf.CDS.sam
ref=$genome-$refscaf

## set up file with tab delimited columns:
# name
# genome fasta base
# ploidy (haploid)
## this is panand_sp_ploidy.txt

for queryscaf in "${queryscafs[@]}"
do

ploidy=1
if [ ! -f ${genome}-${queryscaf}-$genome-$refscaf.CDS.sam ]
then
/programs/minimap2-2.28/minimap2 -x splice -t 20 -k 12 -a -p 0.4 -N 20 ${genome}-${queryscaf}.fa $genome-$refscaf.CDS.fa > ${genome}-${queryscaf}-$genome-$refscaf.CDS.sam 
fi
if [ ! -f ${name}-${ref}-${ploidy} ]
then
## run anchors
./anchorwave/anchorwave/anchorwave proali -t 20 -i $refgff -as $genome-$refscaf.CDS.fa -r $reffa -a ${genome}-${queryscaf}-$genome-$refscaf.CDS.sam -ar $genome-$refscaf.CDS.sam -s ${genome}-${queryscaf}.fa -n ${genome}-${queryscaf}-$genome-$refscaf-1 -R $ploidy -Q 1 -ns -o ${genome}-${queryscaf}-$genome-$refscaf.maf -f ${genome}-${queryscaf}-$genome-$refscaf.f.maf 
fi
done
