## this should be interesting - 5.1% het, no Bionano
## 1698 mb flow
conda activate panand_assemblies

hicR1=AV9240004_Thelepogon_elegans_HiC_pairsR1.fq.gz
hicR2=AV9240004_Thelepogon_elegans_HiC_pairsR2.fq.gz
six=telega

scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240004_Thelepogon_elegans_HiC_pairsR1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240004_Thelepogon_elegans_HiC_pairsR2.fq.gz .
samtools bam2fq -@ 64 ../../../panand_snps/input/telega_reads_withQ.bam > ${six}_hifi.fastq


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


## dual-scaf, hic assembly
# hifiasm -o ${six}.asm -t64 --dual-scaf --h1 $hicR1 --h2 $hicR2 ${six}_hifi.fastq.gz
## i guess try lower s for high het??
hifiasm -o ${six}.asm -t64 -l 0 --dual-scaf --h1 $hicR1 --h2 $hicR2 ${six}_hifi.fastq.gz
## it's still pretty big (3gb per haploytpe, when flow suggests 1698mb for n)
## but i'm losing it on this genome. maybe it endoreduplicated, i don't know. 

awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap1.p_ctg.gfa > ${six}.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap2.p_ctg.gfa > ${six}.asm.hic.hap2.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.p_ctg.gfa > ${six}.asm.hic.p_ctg.fa


## check assembly, how many contigs per chromosome??
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap1.p_ctg.fa" out="${six}.asm.hic.hap1.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap2.p_ctg.fa" out="${six}.asm.hic.hap2.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.p_ctg.fa" out="${six}.asm.hic.p_ctg.stats.txt" overwrite='true'


## check out a dotplot
cd ..
./fa_to_dotplot.sh ${six}/${six}.asm.hic.hap1.p_ctg.fa ${six}Hap1 3
./fa_to_dotplot.sh ${six}/${six}.asm.hic.hap2.p_ctg.fa ${six}Hap2 3
cd ${six}


## combine haplotypes and haphic



cat ${six}.asm.hic.hap1.p_ctg.fa ${six}.asm.hic.hap2.p_ctg.fa >${six}.asm.hic.bothhaps.p_ctg.fa

bwa index ${six}.asm.hic.bothhaps.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.bothhaps.p_ctg.fa AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}BothHaps.bam
## duh it can't index if not sorted
##samtools index ${six}BothHaps.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}BothHaps.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}BothHaps.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
#../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.filtered.bam 20 --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHaps

## skip filtering? potential rec here https://github.com/zengxiaofei/HapHiC/issues/86
## use gfa to be sure not to combine haps
../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.bam 20 --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsUnfiltered



## to correct contigs with hic info
## --correct_nrounds 2 
## to remove hic links between alleles, so set to number of haplotypes present (polyploid)
## --remove_allelic_links 2

## try quickview 
#../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.filtered.bam 20 --quick_view --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsQuickView
../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.bam 20 --quick_view  --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsQuickViewUnfiltered


#########################
## scaffold each haplotype separately

conda activate haphic

bwa index ${six}.asm.hic.hap1.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.hap1.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}Hap1.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}Hap1.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}Hap1.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.hap1.p_ctg.fa ${six}Hap1.filtered.bam 20 --outdir ${six}Hap1

juicer post -o ${six}Hap1_aggressivecorrection ${six}Hap1_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.hap1.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa ${six}Hap1aggressive 3
### dang this thing is really tetraploid - redo scaffolding with 20 chrom





## hap2
bwa index ${six}.asm.hic.hap2.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.hap2.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}Hap2.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}Hap2.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}Hap2.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.hap2.p_ctg.fa ${six}Hap2.filtered.bam 10 --outdir ${six}Hap2


### then use corrected contact maps
# Generate the final FASTA file for the scaffolds
## from the 04 directory...

###
juicer post -o ${six}Hap1_aggressivecorrection ${six}Hap1_aggressivecorrection.out_JBAT.review.assembly.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.hap1.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa ${six}Hap1aggressive 3







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
Rscript generate_subphaser_input_cmdline.R ${six}Hap1aggressive-Pv-6 2 ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa.fai


## prepare for atlas!!!!
cp ${six}/${six}Hap1/04.build/${six}Hap1_aggressivecorrection.FINAL.fa ~/transfer/
cp ${six}Hap1aggressive_subphaserinput.txt ~/transfer/

### ON ATLAS
six=telega
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}Hap1_aggressivecorrection.FINAL.fa .

## trash2 on atlas
## on atlas - not having error?!?!?!?!
cd /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats
mkdir ${six}Hap1_aggressivecorrection_TRASH2
conda activate trash ## i guess do it before???
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=telega; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/${six}Hap1_aggressivecorrection.FINAL.fa -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}Hap1_aggressivecorrection_TRASH2 -p 46"

### then helixer
cd ..
module load apptainer
## first time had to download models!
# /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif 
# Singularity> fetch_helixer_models.py --lineage land_plant
### submit as script!~!!!! 
#### AAAHHH IT DOESN'T SET variable as variable, swicht in wrap stamentment
##six=achine
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=telega; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path ${six}Hap1_aggressivecorrection.FINAL.fa  --lineage land_plant --gff-output-path ${six}Hap1_aggressivecorrection.FINAL.helixer.gff3'


conda activate SubPhaser
cd subphaser
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}Hap1aggressive_subphaserinput.txt .

## generate subphaser in put through my script from anchorwave output (need to improve usability)
six=telegaHap1
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=telegaHap1; subphaser -i ../${six}_aggressivecorrection.FINAL.fa -c ${six}aggressive_subphaserinput.txt -pre ${six}_aggressivecorrection -k 15 -f 2 -q 50 -nsg 2 -non_specific -p 46"

