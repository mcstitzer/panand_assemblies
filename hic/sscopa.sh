## 0.45% het, Bionano, 6% ks btwn dups, 2Gb
conda activate panand_assemblies



hicR1=AV9240003_Schizachyrium_scoparium_HiC_pairsR1.fq.gz
hicR2=AV9240003_Schizachyrium_scoparium_HiC_pairsR2.fq.gz
six=sscopa

scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/*/$hicR1 .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/*/$hicR2 .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/*/AV9590001_Schizachyrium_scoparium_20250822_01_R84064-1_A01.filter.bam .
samtools bam2fq -@ 64 AV9590001_Schizachyrium_scoparium_20250822_01_R84064-1_A01.filter.bam > ${six}_hifi.fastq


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


## dual-scaf, hic assembly, DISABLE PURGING, THIS HAS NO HET!!
## unclear to me whether i shoudl use p_ctg then or the two haplotypes...
 hifiasm -o ${six}.asm -t64 --dual-scaf --h1 $hicR1 --h2 $hicR2 ${six}_hifi.fastq.gz

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






## combine haplotypes and haphic


conda activate haphic

cat ${six}.asm.hic.hap1.p_ctg.fa ${six}.asm.hic.hap2.p_ctg.fa >${six}.asm.hic.bothhaps.p_ctg.fa

bwa index ${six}.asm.hic.bothhaps.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.bothhaps.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}BothHaps.bam
## duh it can't index if not sorted
##samtools index ${six}BothHaps.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}BothHaps.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}BothHaps.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.filtered.bam 40 --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHaps

## skip filtering? potential rec here https://github.com/zengxiaofei/HapHiC/issues/86
## use gfa to be sure not to combine haps
#../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.bam 40 --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsUnfiltered

## can't see homeologs in the hic map with them combined, but that's okay!
###
juicer post -o ${six}BothHaps_aggressivecorrection ${six}BothHaps_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.bothhaps.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}BothHaps/04.build/${six}BothHaps_aggressivecorrection.FINAL.fa ${six}BothHapsaggressive 3

## aamazing!!!

## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
GENOME=${six}/${six}BothHaps/04.build/${six}BothHaps_aggressivecorrection.FINAL.fa
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

cd ${six}/${six}BothHaps/04.build/
### gc content
samtools faidx ${six}BothHaps_aggressivecorrection.FINAL.fa
bedtools makewindows -g ${six}BothHaps_aggressivecorrection.FINAL.fa.fai -w 1000000 > ${six}BothHaps_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi ${six}BothHaps_aggressivecorrection.FINAL.fa -bed ${six}BothHaps_aggressivecorrection.FINAL.1mbwindows.bed > ${six}BothHaps_aggressivecorrection.FINAL.1mbnuccontent.bed

## run tidk
conda activate tidk
tidk find --clade Poales --output ${six}BothHaps_aggressivecorrection.tidk --dir . ${six}BothHaps_aggressivecorrection.FINAL.fa 


cd ../../../
Rscript generate_subphaser_input_cmdline.R ${six}BothHapsaggressive-Pv-6 4 ${six}/${six}BothHaps/04.build/${six}BothHaps_aggressivecorrection.FINAL.fa.fai


## prepare for atlas!!!!
cp ${six}BothHaps_aggressivecorrection.FINAL.fa ~/transfer/
cp ../../../sscopaBothHapsaggressive_subphaserinput.txt ~/transfer/

### ON ATLAS
six=sscopa
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}BothHaps_aggressivecorrection.FINAL.fa .

## trash2 on atlas
## on atlas - not having error?!?!?!?!
cd /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats
mkdir ${six}BothHaps_aggressivecorrection_TRASH2
conda activate trash ## i guess do it before???
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=sscopa; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/${six}BothHaps_aggressivecorrection.FINAL.fa -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}BothHaps_aggressivecorrection_TRASH2 -p 46"

### then helixer
cd ..
module load apptainer
## first time had to download models!
# /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif 
# Singularity> fetch_helixer_models.py --lineage land_plant
### submit as script!~!!!! 
#### AAAHHH IT DOESN'T SET variable as variable, swicht in wrap stamentment
##six=achine
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=sscopa; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path ${six}BothHaps_aggressivecorrection.FINAL.fa  --lineage land_plant --gff-output-path ${six}BothHaps_aggressivecorrection.FINAL.helixer.gff3'


conda activate SubPhaser
cd subphaser
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}BothHapsaggressive_subphaserinput.txt .

## generate subphaser in put through my script from anchorwave output (need to improve usability)
six=sscopaBothHaps
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=sscopaBothHaps; subphaser -i ../${six}_aggressivecorrection.FINAL.fa -c ${six}_aggressivecorrection_subphaserinput.txt -pre ${six}_aggressivecorrection -k 15 -f 2 -q 50 -nsg 2 -non_specific -p 46"












