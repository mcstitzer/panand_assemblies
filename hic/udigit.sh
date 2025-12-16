## very low het, Bionano, a little nanopore 4.2% ks btwn dups, 5.4Gb
conda activate panand_assemblies

hicR1=AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz
hicR2=AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz
six=udigit

scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/*/$hicR1 .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/*/$hicR2 .
samtools bam2fq -@ 64 ../../../panand_snps/input/rtuber_reads_withQ.bam > ${six}_hifi.fastq


## check kmers first!
mkdir -p tmp
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 ${six}_hifi.fastq.gz ${six}kmc tmp/
/programs/kmc-3.2.4/kmc_tools transform ${six}kmc histogram ${six}.histo -cx10000
../genomescope2.0/genomescope.R -i ${six}.histo -o ${six}_genomescope -k 21

export PATH=/programs/fastk-1.1/bin:$PATH
export PATH=/programs/smudgeplot-0.4.0/bin:$PATH
# run FastK to create a k-mer database
FastK -v -t4 -k31 -M16 -T4 ${six}_hifi.fastq -N${six}_FastK_Table
## fastk just unzips it so wait to bgzip for the first time...
#bgzip ${six}_hifi.fastq

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot.py hetmers -L 12 -t 4 -o ${six}_kmerpairs --verbose ${six}_FastK_Table
# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot.py all -o ${six}_smudgeplot ${six}_kmerpairs_text.smu

## very very low het, so l0
hifiasm -o ${six}.asm --dual-scaf -l 0 -t64 --ul /workdir/mcs368/panand_assemblies/scaffolding/udigit/udigit.50000.fq.gz --h1 AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz --h2 AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz udigit_hifi.fastq.gz

awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap1.p_ctg.gfa > ${six}.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap2.p_ctg.gfa > ${six}.asm.hic.hap2.p_ctg.fa


## check assembly, how many contigs per chromosome??
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap1.p_ctg.fa" out="${six}.asm.hic.hap1.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap2.p_ctg.fa" out="${six}.asm.hic.hap2.p_ctg.stats.txt" overwrite='true'


cd ../
./fa_to_dotplot.sh udigit/udigitDip.asm.hic.hap1.p_ctg.fa udigitDipHap1 3
./fa_to_dotplot.sh udigit/udigitDip.asm.hic.hap2.p_ctg.fa udigitDipHap2 3
cd udigit



## very low het on this guy so go ahead and sdo haps separately


#########################
## scaffold each haplotype separately

conda activate haphic

bwa index ${six}.asm.hic.hap1.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.hap1.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}Hap1.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}Hap1.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}Hap1.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.hap1.p_ctg.fa ${six}Hap1.filtered.bam 30 --outdir ${six}Hap1

conda activate panand_assemblies
juicer post -o ${six}Hap1_aggressivecorrection ${six}Hap1_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.hap1.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa ${six}Hap1aggressive 3

juicer post -o ${six}Hap1_aggressivecorrection ${six}Hap1_aggressivecorrection_Try2.out_JBAT.review.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.hap1.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa ${six}Hap1aggressive 3




## hap2
bwa index ${six}.asm.hic.hap2.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.hap2.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}Hap2.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}Hap2.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}Hap2.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.hap2.p_ctg.fa ${six}Hap2.filtered.bam 30 --outdir ${six}Hap2







## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
GENOME=${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

cd ${six}/${six}Hap1/04.build/
### gc content
samtools faidx ${six}Hap1_aggressivecorrection.FINAL.fa
bedtools makewindows -g ${six}Hap1_aggressivecorrection.FINAL.fa.fai -w 1000000 > ${six}Hap1_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi ${six}Hap1_aggressivecorrection.FINAL.fa -bed ${six}Hap1_aggressivecorrection.FINAL.1mbwindows.bed > ${six}Hap1_aggressivecorrection.FINAL.1mbnuccontent.bed

## run tidk
conda activate tidk
tidk find --clade Poales --output ${six}Hap1_aggressivecorrection.tidk --dir . ${six}Hap1_aggressivecorrection.FINAL.fa 


cd ../../../
Rscript generate_subphaser_input_cmdline.R ${six}Hap1aggressive-Pv-6 3 ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa.fai


## prepare for atlas!!!!
cp ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa ~/transfer/
cp ${six}Hap1aggressive_subphaserinput.txt ~/transfer/

### ON ATLAS
six=udigit
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}Hap1_aggressivecorrection.FINAL.fa .

## trash2 on atlas
## on atlas - not having error?!?!?!?!
cd /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats
mkdir ${six}Hap1_aggressivecorrection_TRASH2
conda activate trash ## i guess do it before???
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=udigit; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/${six}Hap1_aggressivecorrection.FINAL.fa -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}Hap1_aggressivecorrection_TRASH2 -p 46"

### then helixer
cd ..
module load apptainer
## first time had to download models!
# /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif 
# Singularity> fetch_helixer_models.py --lineage land_plant
### submit as script!~!!!! 
#### AAAHHH IT DOESN'T SET variable as variable, swicht in wrap stamentment
##six=achine
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=udigit; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path ${six}Hap1_aggressivecorrection.FINAL.fa  --lineage land_plant --gff-output-path ${six}Hap1_aggressivecorrection.FINAL.helixer.gff3'


conda activate SubPhaser
cd subphaser
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}Hap1aggressive_subphaserinput.txt .

## generate subphaser in put through my script from anchorwave output (need to improve usability)
six=sscopaHap1
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=udigitHap1; subphaser -i ../${six}_aggressivecorrection.FINAL.fa -c ${six}aggressive_subphaserinput.txt -pre ${six}_aggressivecorrection -k 15 -f 2 -q 50 -nsg 2 -non_specific -p 46"


