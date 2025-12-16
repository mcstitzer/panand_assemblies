## this should be a fun try - basically 0 het, no Bionano
## 1751 mb flow
conda activate panand_assemblies

hicR1=AV9240010_Rottboellia_tuberculosa_HiC_pairsR1.fq.gz
hicR2=AV9240010_Rottboellia_tuberculosa_HiC_pairsR2.fq.gz
six=rtuber

scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/*/$hicR1 .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/*/$hicR2 .
samtools bam2fq -@ 64 ../../../panand_snps/input/rtuber_reads_withQ.bam > ${six}_hifi.fastq


## check kmers first!
mkdir -p tmp
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 ${six}_hifi.fastq ${six}kmc tmp/
/programs/kmc-3.2.4/kmc_tools transform ${six}kmc histogram ${six}.histo -cx10000
../genomescope2.0/genomescope.R -i ${six}.histo -o ${six}_genomescope -k 21

export PATH=/programs/fastk-1.1/bin:$PATH
export PATH=/programs/smudgeplot-0.4.0/bin:$PATH
# run FastK to create a k-mer database
FastK -v -t4 -k31 -M16 -T4 ${six}_hifi.fastq -N${six}_FastK_Table
## fastk just unzips it so wait to bgzip for the first time...
bgzip ${six}_hifi.fastq

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot.py hetmers -L 12 -t 4 -o ${six}_kmerpairs --verbose ${six}_FastK_Table
# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot.py all -o ${six}_smudgeplot ${six}_kmerpairs_text.smu


## dual-scaf, hic assembly, DISABLE PURGING, THIS HAS NO HET!!
## unclear to me whether i shoudl use p_ctg then or the two haplotypes...
 hifiasm -o ${six}.asm -t64 -l 0 --dual-scaf --h1 $hicR1 --h2 $hicR2 ${six}_hifi.fastq.gz

awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap1.p_ctg.gfa > ${six}.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap2.p_ctg.gfa > ${six}.asm.hic.hap2.p_ctg.fa


## check assembly, how many contigs per chromosome??
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap1.p_ctg.fa" out="${six}.asm.hic.hap1.p_ctg.stats.txt"
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap2.p_ctg.fa" out="${six}.asm.hic.hap2.p_ctg.stats.txt"
## hap1 925 Mb, 882 scaffolds, N50 82 Mb 5 scaffold N50
## hap2 899 Mb, 158 scaffolds, N50 83 Mb 5 scaffold N50

## check out a dotplot
cd ..
./fa_to_dotplot.sh ${six}/${six}.asm.hic.hap1.p_ctg.fa ${six}Hap1 1
./fa_to_dotplot.sh ${six}/${six}.asm.hic.hap2.p_ctg.fa ${six}Hap2 1
cd ${six}




#########################
## scaffold each haplotype separately

conda activate haphic

bwa index ${six}.asm.hic.hap1.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.hap1.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}Hap1.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}Hap1.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}Hap1.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.hap1.p_ctg.fa ${six}Hap1.filtered.bam 9 --outdir ${six}Hap1

## hap2
bwa index ${six}.asm.hic.hap2.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.hap2.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}Hap2.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}Hap2.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}Hap2.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.hap2.p_ctg.fa ${six}Hap2.filtered.bam 9 --outdir ${six}Hap2


### then use corrected contact maps
# Generate the final FASTA file for the scaffolds
## from the 04 directory...
python ../../../juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py -a out_JBAT.review.${six}Hap2.assembly -f scaffolds.fa -s

###
juicer post -o ${six}Hap1_aggressivecorrection ${six}Hap1_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.hap1.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa ${six}Hap1aggressive 1


## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp

GENOME=${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa 
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

## run trash2
conda activate TRASH_2
mkdir -p ${six}Hap1_aggressivecorrection_TRASH2
Rscript /workdir/mcs368/panand_assemblies/repeats/TRASH_2/src/TRASH.R -f ${six}Hap1_aggressivecorrection.FINAL.fa -o /workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection_TRASH2 -p 48

## trash2 on atlas
## on atlas - not having error?!?!?!?!
cd /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats
six=rtuber
mkdir -p ${six}Hap1_aggressivecorrection_TRASH2
conda activate trash ## i guess do it before???
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=rtuber; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/${six}Hap1_aggressivecorrection.FINAL.fa -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}Hap1_aggressivecorrection_TRASH2 -p 46"


### or if trash2 is stupid and fails, do trash
conda activate panand_assemblies
mkdir -p /workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection_TRASH
/workdir/mcs368/panand_assemblies/repeats/TRASH/TRASH_run.sh --def /workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa --par 48 --o /workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection_TRASH



## run tidk
conda activate tidk
tidk find --clade Poales --output ${six}Hap1_aggressivecorrection.tidk --dir . ${six}Hap1_aggressivecorrection.FINAL.fa 

## run helixer
ssh cbsugpu08
cd /workdir/mcs368/
scp cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa .
singularity run --nv --bind $PWD --pwd $PWD /programs/helixer-0.3.5/helixer.sif Helixer.py --fasta-path ${six}Hap1_aggressivecorrection.FINAL.fa --lineage land_plant --gff-output-path ${six}Hap1_aggressivecorrection.FINAL.helixer.gff3
scp ${six}Hap1_aggressivecorrection.FINAL.helixer.gff3 cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.helixer.gff3

### OR ON ATLAS
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}Hap1_aggressivecorrection.FINAL.fa .
srun -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --pty bash
module load apptainer
## first time had to download models!
# /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif 
# Singularity> fetch_helixer_models.py --lineage land_plant
time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path rtuberHap1_aggressivecorrection.FINAL.fa  --lineage land_plant --gff-output-path rtuberHap1_aggressivecorrection.FINAL.helixer.gff3
scp ${six}Hap1_aggressivecorrection.FINAL.helixer.gff3 mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/


### gc content
conda activate panand_assemblies
samtools faidx ${six}Hap1_aggressivecorrection.FINAL.fa
bedtools makewindows -g ${six}Hap1_aggressivecorrection.FINAL.fa.fai -w 1000000 > ${six}Hap1_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi ${six}Hap1_aggressivecorrection.FINAL.fa -bed ${six}Hap1_aggressivecorrection.FINAL.1mbwindows.bed > ${six}Hap1_aggressivecorrection.FINAL.1mbnuccontent.bed





############# NOW DO IT ALL WITH HAP2!!!!
###
juicer post -o ${six}Hap2_aggressivecorrection ${six}Hap2_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.hap2.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection.FINAL.fa ${six}Hap2aggressive 1



## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
GENOME=${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection.FINAL.fa 
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

## run trash2
conda activate TRASH_2
cd ${six}/${six}Hap2/04.build
mkdir -p ${six}Hap2_aggressivecorrection_TRASH2
cd /workdir/mcs368/panand_assemblies/repeats/TRASH_2/src/
Rscript /workdir/mcs368/panand_assemblies/repeats/TRASH_2/src/TRASH.R -f /workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection.FINAL.fa -o /workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection_TRASH2 -p 24

## run tidk
conda activate tidk
tidk find --clade Poales --output ${six}Hap2_aggressivecorrection.tidk --dir . ${six}Hap2_aggressivecorrection.FINAL.fa 

## run helixer
ssh cbsugpu08
cd /workdir/mcs368/
scp cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection.FINAL.fa .
singularity run --nv --bind $PWD --pwd $PWD /programs/helixer-0.3.5/helixer.sif Helixer.py --fasta-path ${six}Hap2_aggressivecorrection.FINAL.fa --lineage land_plant --gff-output-path ${six}Hap2_aggressivecorrection.FINAL.helixer.gff3
scp ${six}Hap2_aggressivecorrection.FINAL.helixer.gff3 cbsuxm01:/workdir/mcs368/panand_assemblies/hic/${six}/${six}Hap2/04.build/${six}Hap2_aggressivecorrection.FINAL.helixer.gff3

### OR ON ATLAS
cp ${six}Hap2_aggressivecorrection.FINAL.fa ~/transfer/

scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}Hap2_aggressivecorrection.FINAL.fa .
srun -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --pty bash
module load apptainer
## first time had to download models!
# /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif 
# Singularity> fetch_helixer_models.py --lineage land_plant
time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path rtuberHap2_aggressivecorrection.FINAL.fa  --lineage land_plant --gff-output-path rtuberHap2_aggressivecorrection.FINAL.helixer.gff3
scp ${six}Hap2_aggressivecorrection.FINAL.helixer.gff3 mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/




### gc content
samtools faidx ${six}Hap2_aggressivecorrection.FINAL.fa
bedtools makewindows -g ${six}Hap2_aggressivecorrection.FINAL.fa.fai -w 1000000 > ${six}Hap2_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi ${six}Hap2_aggressivecorrection.FINAL.fa -bed ${six}Hap2_aggressivecorrection.FINAL.1mbwindows.bed > ${six}Hap2_aggressivecorrection.FINAL.1mbnuccontent.bed
