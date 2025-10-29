
conda activate panand_assemblies
## first time through, hard coded for udigit

mkdir -p udigit
cd udigit
## get data
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/PanAnd_082025/AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/PanAnd_082025/AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz .

bwa index Ud-Pasquet1171-DRAFT-PanAnd-1.0.fasta

## 5SP deals with weird hic distribution of reads
bwa mem -t 64 -5SP Ud-Pasquet1171-DRAFT-PanAnd-1.0.fasta AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o udigit_hic.sorted.bam 
samtools index udigit_hic.sorted.bam

## get site positions of restriction enzyme on scaffolds
python ../juicer/misc/generate_site_positions.py DpnII GATC Ud-Pasquet1171-DRAFT-PanAnd-1.0.fasta 



## try a easier genome come on michelle
## smicro
mkdir -p smicro
cd smicro
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz .

bwa index Sm-PI203595-DRAFT-PanAnd-1.0.fasta

## 5SP deals with weird hic distribution of reads
bwa mem -t 64 -5SP Sm-PI203595-DRAFT-PanAnd-1.0.fasta AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o smicro_hic.sorted.bam 
samtools index smicro_hic.sorted.bam

## get site positions of restriction enzyme on scaffolds
python ../juicer/misc/generate_site_positions.py DpnII GATC Sm-PI203595-DRAFT-PanAnd-1.0.fasta 






