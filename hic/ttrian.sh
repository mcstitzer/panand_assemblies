## this should be relatively easy - 0.175% het, no Bionano
## 854 mb flow
conda activate panand_assemblies

hicR1=AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz
hicR2=AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz
six=ttrian

scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082625_PanAnd/AV9590002_Themeda_triandra_HiFi_20250822_01_R84064-1_B01.filter.bam .
samtools bam2fq -@ 64 AV9590002_Themeda_triandra_HiFi_20250822_01_R84064-1_B01.filter.bam > ${six}_hifi.fastq
bgzip ${six}_hifi.fastq



## check kmers first!
mkdir tmp
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 ${six}_hifi.fastq.gz ${six}kmc tmp/
/programs/kmc-3.2.4/kmc_tools transform ${six}kmc histogram ${six}.histo -cx10000

../genomescope2.0/genomescope.R -i ${six}.histo -o ${six}_genomescope -k 21

export PATH=/programs/fastk-1.1/bin:$PATH
export PATH=/programs/smudgeplot-0.4.0/bin:$PATH
# run FastK to create a k-mer database
FastK -v -t4 -k31 -M16 -T4 ${six}_hifi.fastq.gz -N${six}_FastK_Table

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot.py hetmers -L 12 -t 4 -o ${six}_kmerpairs --verbose ${six}_FastK_Table
# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot.py all -o ${six}_smudgeplot ${six}_kmerpairs_text.smu


## dual-scaf, hic assembly
# hifiasm -o ${six}.asm -t64 --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
### this looks really weird - big assembly, also looks like a tetraploid? two paths in each haplotype visible in hic....
### hmm, maybe it's not getting homozygous coverage peak right?
# hifiasm -o ${six}.asm -t64 --hom-cov 22 --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
## back to normal because it was duplicating everyting???
#hifiasm -o ${six}.asm -t64  --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
## still big! try --hg-size to set haploid size?
#hifiasm -o ${six}.asm -t64  --dual-scaf --hg-size 850m --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
## NOPE still big
##ummm chatgpt told me to use 55, so i'll try it...
# hifiasm -o ${six}.asm -t64 --hom-cov 48 --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz

### this is beautiful!!!
 hifiasm -o ${six}.asm -t64 -s 0.2 --dual-scaf --h1 AV9240007_Themeda_triandra_HiC_pairsR1.fq.gz --h2 AV9240007_Themeda_triandra_HiC_pairsR2.fq.gz ${six}_hifi.fastq.gz
## 1.6 Gb a piece!


awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap1.p_ctg.gfa > ${six}.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.hap2.p_ctg.gfa > ${six}.asm.hic.hap2.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' ${six}.asm.hic.p_ctg.gfa > ${six}.asm.hic.p_ctg.fa


## check assembly, how many contigs per chromosome??
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap1.p_ctg.fa" out="${six}.asm.hic.hap1.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.hap2.p_ctg.fa" out="${six}.asm.hic.hap2.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="${six}.asm.hic.p_ctg.fa" out="${six}.asm.hic.p_ctg.stats.txt" overwrite='true'



### I DON"T GET IT... this is still huge, and what is happenoing??
# i'm going to build a tree using these assemblies and see if they're what we think they are?
## or like ahve they polyploidized??

## check out a dotplot
cd ..
./fa_to_dotplot.sh ${six}/${six}.asm.hic.hap1.p_ctg.fa ${six}Hap1 1
./fa_to_dotplot.sh ${six}/${six}.asm.hic.hap2.p_ctg.fa ${six}Hap2 1
cd ${six}


## combine haplotypes and haphic



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
../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.bam 40 --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsUnfiltered



## to correct contigs with hic info
## --correct_nrounds 2 
## to remove hic links between alleles, so set to number of haplotypes present (polyploid)
## --remove_allelic_links 2

## try quickview 
#../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.filtered.bam 20 --quick_view --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsQuickView
../HapHiC/haphic pipeline ${six}.asm.hic.bothhaps.p_ctg.fa ${six}BothHaps.bam 40 --quick_view  --gfa "${six}.asm.hic.hap1.p_ctg.gfa,${six}.asm.hic.hap2.p_ctg.gfa" --outdir ${six}BothHapsQuickViewUnfiltered


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
cd ${six}/${six}Hap1/04.build/
bash juicebox.sh


## hap2
bwa index ${six}.asm.hic.hap2.p_ctg.fa
bwa mem -5SP -t 128 ${six}.asm.hic.hap2.p_ctg.fa $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}Hap2.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}Hap2.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}Hap2.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline ${six}.asm.hic.hap2.p_ctg.fa ${six}Hap2.filtered.bam 20 --outdir ${six}Hap2
cd ${six}/${six}Hap2/04.build/
bash juicebox.sh


### then use corrected contact maps
# Generate the final FASTA file for the scaffolds
## from the 04 directory...
##python ../../../juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py -a out_JBAT.review.${six}Hap2.assembly -f scaffolds.fa -s


##### download out_JBAT.hic and out_JBAT.assembly, and correct in juicebox


juicer post -o ${six}BothHapsUnfiltered_aggressivecorrection ${six}BothHapsUnfiltered_aggressivecorrection.out_JBAT.review.assembly.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.bothhaps.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}BothHapsUnfiltered/04.build/${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa ${six}BothHapsUnfilteredaggressive 3


juicer post -o ${six}BothHapsUnfiltered_aggressivecorrection ${six}BothHapsUnfilteredTry2_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../${six}.asm.hic.bothhaps.p_ctg.fa
cd ../../../
./fa_to_dotplot.sh ${six}/${six}BothHapsUnfiltered/04.build/${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa ${six}BothHapsUnfilteredaggressive 3




## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
GENOME=${six}/${six}BothHapsUnfiltered/04.build/${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

cd ${six}/${six}BothHapsUnfiltered/04.build/
### gc content
samtools faidx ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa
bedtools makewindows -g ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa.fai -w 1000000 > ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa -bed ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.1mbwindows.bed > ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.1mbnuccontent.bed

## run tidk
conda activate tidk
tidk find --clade Poales --output ${six}BothHapsUnfiltered_aggressivecorrection.tidk --dir . ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa 

## added 5 as a ploidy option!!
cd ../../../
Rscript generate_subphaser_input_cmdline.R ${six}BothHapsUnfilteredaggressive-Pv-6 5 ${six}/${six}BothHapsUnfiltered/04.build/${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa.fai


## prepare for atlas!!!!
cd ${six}/${six}BothHapsUnfiltered/04.build/
cp ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa ~/transfer/
cp ../../../${six}BothHapsUnfilteredaggressive_subphaserinput.txt ~/transfer/

### ON ATLAS
six=ttrian

scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa .

## trash2 on atlas
## on atlas - not having error?!?!?!?!
cd /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats
mkdir ${six}BothHapsUnfiltered_aggressivecorrection_TRASH2
conda activate trash ## i guess do it before???
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=ttrianBothHapsUnfiltered; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}BothHapsUnfiltered_aggressivecorrection_TRASH2 -p 46"

### then helixer
cd ..
module load apptainer
## first time had to download models!
# /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif 
# Singularity> fetch_helixer_models.py --lineage land_plant
### submit as script!~!!!! 
#### AAAHHH IT DOESN'T SET variable as variable, swicht in wrap stamentment
##six=achine
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=ttrian; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.fa  --lineage land_plant --gff-output-path ${six}BothHapsUnfiltered_aggressivecorrection.FINAL.helixer.gff3'


conda activate SubPhaser
cd subphaser
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}BothHapsUnfilteredaggressive_subphaserinput.txt .

## generate subphaser in put through my script from anchorwave output (need to improve usability)
six=ttrianBothHapsUnfiltered
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=ttrianBothHapsUnfiltered; subphaser -i ../${six}_aggressivecorrection.FINAL.fa -c ${six}_aggressivecorrection_subphaserinput.txt -pre ${six}_aggressivecorrection -k 15 -f 2 -q 50 -nsg 2 -non_specific -p 46"
## tried 5 subgenomes for 19 and 21... 2 for all others


#####also map hic to existing referenceand try toscaffold??!?!


conda activate haphic

bwa index Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta
bwa mem -5SP -t 96 Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta $hicR1 $hicR2 | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o ${six}DIPLOID.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam ${six}DIPLOID.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o ${six}DIPLOID.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## well, it's allelic so i guess we say 40 chromsomes!!!
../HapHiC/haphic pipeline Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta ${six}DIPLOID.filtered.bam 10 --outdir ${six}DIPLOID

## try quickview
../HapHiC/haphic pipeline Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta ${six}DIPLOID.filtered.bam 10 --outdir ${six}DIPLOIDqv --quick_view



### hmmm not the greatest - adding a lot of small (100kb) contigs to ends of each chromosome
### just be aggressive :)
juicer post -o ${six}DIPLOIDqv_aggressivecorrection ${six}DIPLOIDqv_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta

cd ../../../
./fa_to_dotplot.sh ${six}/${six}DIPLOIDqv/04.build/${six}DIPLOIDqv_aggressivecorrection.FINAL.fa ${six}DIPLOIDqv_aggressive 3

juicer post -o ${six}DIPLOIDqv_aggressivecorrection ${six}DIPLOIDqv_aggressivecorrectionTry2.out_JBAT.review.assembly out_JBAT.liftover.agp ../../Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta
cd ../../../
./fa_to_dotplot.sh ${six}/${six}DIPLOIDqv/04.build/${six}DIPLOIDqv_aggressivecorrection.FINAL.fa ${six}DIPLOIDqv_aggressive 3

juicer post -o ${six}DIPLOIDqv_aggressivecorrection ${six}DIPLOIDqv_aggressivecorrectionTry3.out_JBAT.review.assembly out_JBAT.liftover.agp ../../Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta
cd ../../../
./fa_to_dotplot.sh ${six}/${six}DIPLOIDqv/04.build/${six}DIPLOIDqv_aggressivecorrection.FINAL.fa ${six}DIPLOIDqv_aggressive 3

## kahy this is getting ridiculous....
juicer post -o ${six}DIPLOIDqv_aggressivecorrection ${six}DIPLOIDqv_aggressivecorrectionTry4.out_JBAT.review.assembly out_JBAT.liftover.agp ../../Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta
cd ../../../
./fa_to_dotplot.sh ${six}/${six}DIPLOIDqv/04.build/${six}DIPLOIDqv_aggressivecorrection.FINAL.fa ${six}DIPLOIDqv_aggressive 3



## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
GENOME=${six}/${six}DIPLOIDqv/04.build/${six}DIPLOIDqv_aggressivecorrection.FINAL.fa
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

cd ${six}/${six}DIPLOIDqv/04.build/

### gc content
samtools faidx ${six}DIPLOIDqv_aggressivecorrection.FINAL.fa
bedtools makewindows -g ${six}DIPLOIDqv_aggressivecorrection.FINAL.fa.fai -w 1000000 > ${six}DIPLOIDqv_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi ${six}DIPLOIDqv_aggressivecorrection.FINAL.fa -bed ${six}DIPLOIDqv_aggressivecorrection.FINAL.1mbwindows.bed > ${six}DIPLOIDav_aggressivecorrection.FINAL.1mbnuccontent.bed

## run tidk
conda activate tidk
tidk find --clade Poales --output ${six}DIPLOIDqv_aggressivecorrection.tidk --dir . ${six}DIPLOIDqv_aggressivecorrection.FINAL.fa 


## atlas

### OR ON ATLAS
cp ${six}DIPLOIDqv_aggressivecorrection.FINAL.fa ~/transfer/
##reload six on atlas
scp mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/${six}DIPLOIDqv_aggressivecorrection.FINAL.fa .
module load apptainer
### submit as script!~!!!! 
#### AAAHHH IT DOESN'T SET variable as variable, swicht in wrap stamentment
##six=achine
sbatch -A buckler_lab_panand -p gpu-a100 --gres=gpu:a100:1 --ntasks-per-node=16 --time=1-00:00 --wrap='six=ttrianDIPLOIDqv; module load apptainer; time apptainer exec --nv /project/buckler_lab_panand/zachary.miller/helixerDocker/helixer-docker_helixer_v0.3.2_cuda_11.8.0-cudnn8.sif Helixer.py --fasta-path ${six}_aggressivecorrection.FINAL.fa  --lineage land_plant --gff-output-path ${six}_aggressivecorrection.FINAL.helixer.gff3'
## once done, scp back!
scp ${six}_aggressivecorrection.FINAL.helixer.gff3 mcs368@cbsulogin2.biohpc.cornell.edu:~/transfer/

## trash2 on atlas
## on atlas - not having error?!?!?!?!
cd /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats
six=ttrianDIPLOIDqv
mkdir ${six}_aggressivecorrection_TRASH2
conda activate trash ## i guess do it before???
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=ttrianDIPLOIDqv; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/${six}_aggressivecorrection.FINAL.fa -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}_aggressivecorrection_TRASH2 -p 46"







