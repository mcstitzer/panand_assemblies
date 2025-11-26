## this just has clr, so take arun's contigs and scaffold with hic!

conda activate panand_assemblies

hicR1=AV9240011_Sorghastrum_nutans_HiC_pairsR1.fq.gz 
hicR2=AV9240011_Sorghastrum_nutans_HiC_pairsR2.fq.gz 
six=snutan

scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/*/$hicR1 .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/*/$hicR2 .


conda activate haphic

bwa index Sn-CAM1369-DRAFT-PanAnd-1.0.fasta
bwa mem -5SP -t 128 Sn-CAM1369-DRAFT-PanAnd-1.0.fasta $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## well, it's allelic so i guess we say 40 chromsomes!!!
../HapHiC/haphic pipeline Sn-CAM1369-DRAFT-PanAnd-1.0.fasta ${six}.filtered.bam 40 --outdir ${six}

### hmmm not the greatest - adding a lot of small (100kb) contigs to ends of each chromosome
### just be aggressive :)
juicer post -o ${six}_aggressivecorrection ${six}_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../Sn-CAM1369-DRAFT-PanAnd-1.0.fasta
cd ../../../
./fa_to_dotplot.sh ${six}/${six}/04.build/${six}_aggressivecorrection.FINAL.fa ${six}aggressive 3




## try yahs? NO DON"T DO THIS!!!?!?!?!?!!?!?!?
# conda activate panand_assemblies
# bwa mem -t 128 -5SP Sn-CAM1369-DRAFT-PanAnd-1.0.fasta ${hicR1} ${hicR2} | samblaster | samtools sort -@ 16 -o ${six}_hic.sorted.bam 
# samtools index ${six}_hic.sorted.bam
# 
# 
# yahs Sn-CAM1369-DRAFT-PanAnd-1.0.fasta ${six}_hic.sorted.bam -o ${six}yahs -e GATC
# 
# # juicer pre smicroHap2.bin smicroHap2_scaffolds_final.agp smicro.asm.hic.hap2.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part
# # mv alignments_sorted.txt.part alignments_sorted.txt
# # (java -jar -Xmx32G ../juicer_tools.2.18.00.jar  pre alignments_sorted.txt out.hic.part smicroHap2_scaffolds_final.fa.fai ) && (mv out.hic.part out.hic)
# 
# juicer pre -a -o snutanyahsout_JBAT ../snutan_hic.sorted.bam snutanyahs_scaffolds_final.agp ../Sn-CAM1369-DRAFT-PanAnd-1.0.fasta.fai >out_JBAT.log 2>&1
#  
# 



# ./fa_to_dotplot.sh smicro/smicro_scaffolds_final.fa smicroHap1yahs 1
# ./fa_to_dotplot.sh smicro/smicroHap2_scaffolds_final.fa smicroHap1yahs 2
# 
# juicer pre smicroHap2.bin smicroHap2_scaffolds_final.agp smicro.asm.hic.hap2.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part
# mv alignments_sorted.txt.part alignments_sorted.txt
# (java -jar -Xmx32G ../juicer_tools.2.18.00.jar  pre alignments_sorted.txt out.hic.part smicroHap2_scaffolds_final.fa.fai ) && (mv out.hic.part out.hic)
# 
# juicer pre -a -o smicroHap1out_JBAT smicroHap1_hic.sorted.bam smicro_scaffolds_final.agp smicro.asm.hic.hap1.p_ctg.fa.fai >out_JBAT.log 2>&1
 

######### after juicebox!!!!

juicer post -o ${six}_aggressivecorrection ${six}_aggressivecorrectionTry2.out_JBAT.review.assembly out_JBAT.liftover.agp ../../Sn-CAM1369-DRAFT-PanAnd-1.0.fasta
cd ../../../
./fa_to_dotplot.sh ${six}/${six}/04.build/${six}_aggressivecorrection.FINAL.fa ${six}aggressive 2



## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
GENOME=${six}/${six}/04.build/${six}_aggressivecorrection.FINAL.fa 
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

## run trash2
conda activate TRASH_2
cd ${six}/${six}/04.build
mkdir -p ${six}_aggressivecorrection_TRASH2
Rscript /workdir/mcs368/panand_assemblies/repeats/TRASH_2/src/TRASH.R -f ${six}_aggressivecorrection.FINAL.fa -o /workdir/mcs368/panand_assemblies/hic/${six}/${six}/04.build/${six}_aggressivecorrection_TRASH2 -p 24

## try on atlas??
conda activate trash
mkdir -p ${six}_aggressivecorrection_TRASH2
#Rscript /workdir/mcs368/panand_assemblies/repeats/TRASH_2/src/TRASH.R -f ${six}_aggressivecorrection.FINAL.fa -o /workdir/mcs368/panand_assemblies/hic/${six}/${six}/04.build/${six}_aggressivecorrection_TRASH2 -p 24
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=1-00:00 --wrap="conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/${six}_aggressivecorrection.FINAL.fa -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/snutan_aggressivecorrection_TRASH2 -p 48"


## run tidk
conda activate tidk
cd ${six}/${six}/04.build/
tidk find --clade Poales --output ${six}_aggressivecorrection.tidk --dir . ${six}_aggressivecorrection.FINAL.fa 

## run helixer
ssh cbsugpu08
cd /workdir/mcs368/
scp cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}/04.build/${six}_aggressivecorrection.FINAL.fa .
singularity run --nv --bind $PWD --pwd $PWD /programs/helixer-0.3.5/helixer.sif Helixer.py --fasta-path ${six}_aggressivecorrection.FINAL.fa --lineage land_plant --gff-output-path ${six}_aggressivecorrection.FINAL.helixer.gff3
scp ${six}_aggressivecorrection.FINAL.helixer.gff3 cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}/04.build/${six}_aggressivecorrection.FINAL.helixer.gff3
### OR ON ATLAS
cp /workdir/mcs368/panand_assemblies/hic/${six}/${six}/04.build/${six}_aggressivecorrection.FINAL.fa ~/transfer/
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}Hap1_aggressivecorrection.FINAL.fa .
srun -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --pty bash
module load apptainer
## first time had to download models!
# /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif 
# Singularity> fetch_helixer_models.py --lineage land_plant
six=snutan
time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path ${six}_aggressivecorrection.FINAL.fa  --lineage land_plant --gff-output-path ${six}_aggressivecorrection.FINAL.helixer.gff3
scp ${six}Hap1_aggressivecorrection.FINAL.helixer.gff3 mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/


### gc content
samtools faidx ${six}_aggressivecorrection.FINAL.fa
bedtools makewindows -g ${six}_aggressivecorrection.FINAL.fa.fai -w 1000000 > ${six}_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi ${six}_aggressivecorrection.FINAL.fa -bed ${six}_aggressivecorrection.FINAL.1mbwindows.bed > ${six}_aggressivecorrection.FINAL.1mbnuccontent.bed
## oh, this also gets me n's!!! plot those too