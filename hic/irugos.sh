## this should be relatively easy - 0.4% het, no Bionano
conda activate panand_assemblies


## dual-scaf, hic assembly
hifiasm -o irugos.asm -t64 --n-hap 2 --dual-scaf --h1 AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz --h2 AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz irugos_hifi.fastq.gz

## also try -l3 because it's pretty dang inbred
#hifiasm -o irugos.asm -t64 -l 3 --h1 AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz --h2 AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz irugos_hifi.fastq.gz

## why does this assembly have more heterozygous bp than homozygous?????!?!?!??!?!?!!!
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 irugos_hifi.fastq.gz irugoskmc tmp/
/programs/kmc-3.2.4/kmc_tools transform irugoskmc histogram irugos.histo -cx10000

../genomescope2.0/genomescope.R -i irugos.histo -o irugos_genomescope -k 21

## there are definitely two kmer peaks there. this is really werid.
six=irugos
export PATH=/programs/fastk-1.1/bin:$PATH
export PATH=/programs/smudgeplot-0.4.0/bin:$PATH
# run FastK to create a k-mer database
FastK -v -t4 -k31 -M16 -T4 ${six}_hifi.fastq.gz -N${six}_FastK_Table
## fastk just unzips it so wait to bgzip for the first time...
#bgzip ${six}_hifi.fastq

# Find all k-mer pairs in the dataset using hetmer module
smudgeplot.py hetmers -L 12 -t 4 -o ${six}_kmerpairs --verbose ${six}_FastK_Table
# this now generated `data/Scer/kmerpairs_text.smu` file;
# it's a flat file with three columns; covB, covA and freq (the number of k-mer pairs with these respective coverages)

# use the .smu file to infer ploidy and create smudgeplot
smudgeplot.py all -o ${six}_smudgeplot ${six}_kmerpairs_text.smu


awk '/^S/{print ">"$2"\n"$3}' irugos.asm.hic.hap1.p_ctg.gfa > irugos.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' irugos.asm.hic.hap2.p_ctg.gfa > irugos.asm.hic.hap2.p_ctg.fa
## what about unphased contigs??
#awk '/^S/{print ">"$2"\n"$3}' irugos.asm.hic.p_ctg.gfa > irugos.asm.hic.p_ctg.fa

## check assembly, how many contigs per chromosome??
/programs/bbmap-38.96/stats.sh -Xmx200g in="irugos.asm.hic.hap1.p_ctg.fa" out="irugos.asm.hic.hap1.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="irugos.asm.hic.hap2.p_ctg.fa" out="irugos.asm.hic.hap2.p_ctg.stats.txt" overwrite='true'
## hap1 619 Mb, 590 scaffolds, N50 3.6 Mb 45 scaffold N50
#/programs/bbmap-38.96/stats.sh -Xmx200g in="irugos.asm.hic.p_ctg.fa" out="irugos.asm.hic.p_ctg.stats.txt" overwrite='true'

## check out a dotplot
cd ..
./fa_to_dotplot.sh irugos/irugos.asm.hic.hap1.p_ctg.fa irugosHap1 1
./fa_to_dotplot.sh irugos/irugos.asm.hic.hap2.p_ctg.fa irugosHap2 1
./fa_to_dotplot.sh irugos/irugos.asm.hic.p_ctg.fa irugos_p_ctg 1

cd irugos

#### jointly scaffold gfas together
cat irugos.asm.hic.hap1.p_ctg.fa irugos.asm.hic.hap2.p_ctg.fa > irugos.asm.hic.bothhaps.p_ctg.fa
bwa index irugos.asm.hic.bothhaps.p_ctg.fa
bwa mem -5SP -t 64 irugos.asm.hic.bothhaps.p_ctg.fa AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o irugos_bothhapsHiC.bam

../HapHiC/utils/filter_bam irugos_bothhapsHiC.bam 1 --nm 3 --threads 14 | samtools view - -b -@ 14 -o irugos_bothhapsHiC.filtered.bam

cd irugosBothHaps/04.build
bash juicebox.sh



conda activate haphic
# HapHic scaffolding - UNFILTERED DIDN"T WORK!!! 
../HapHiC/haphic pipeline irugos.asm.hic.bothhaps.p_ctg.fa irugos_bothhapsHiC.filtered.bam 20 --threads 32  --outdir irugosBothHaps --gfa "irugos.asm.hic.hap1.p_ctg.gfa,irugos.asm.hic.hap2.p_ctg.gfa"

../HapHiC/haphic pipeline irugos.asm.hic.bothhaps.p_ctg.fa irugos_bothhapsHiC.filtered.bam 18 --threads 32  --outdir irugosBothHaps18 --gfa "irugos.asm.hic.hap1.p_ctg.gfa,irugos.asm.hic.hap2.p_ctg.gfa"


# Generate the final FASTA file for the scaffolds
conda activate panand_assemblies
juicer post -o irugosBothHaps_aggressivecorrection irugosumBothHaps_aggressivecorrection.out_JBAT.review.assembly out_JBAT.liftover.agp ../../irugos.asm.hic.bothhaps.p_ctg.fa


### OKAY YAY
### NOW DO STUFF!!!!


## ribosomes
## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
GENOME=irugos/irugosBothHaps/04.build/irugosBothHaps_aggressivecorrection.FINAL.fa
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENOME.tmp
barrnap --quiet --kingdom euk $GENOME.tmp > ${GENOME%.fa}.rrna.gff3
rm $GENOME.tmp

## run trash2
conda activate TRASH_2
cd irugos/irugosBothHaps/04.build
mkdir -p irugosBothHaps_aggressivecorrection_TRASH2
Rscript /workdir/mcs368/panand_assemblies/repeats/TRASH_2/src/TRASH.R -f /workdir/mcs368/panand_assemblies/hic/irugos/irugosBothHaps/04.build/irugosBothHaps_aggressivecorrection.FINAL.fa -o /workdir/mcs368/panand_assemblies/hic/irugos/irugosBothHaps/04.build/irugosBothHaps_aggressivecorrection_TRASH2 -p 24


### or if trash2 is stupid and fails, do trash
conda activate panand_assemblies
six=irugos
mkdir -p /workdir/mcs368/panand_assemblies/hic/${six}/${six}BothHaps/04.build/${six}BothHaps_aggressivecorrection_TRASH
cd /workdir/mcs368/panand_assemblies/repeats/TRASH/
/workdir/mcs368/panand_assemblies/repeats/TRASH/TRASH_run.sh --def /workdir/mcs368/panand_assemblies/hic/${six}/${six}BothHaps/04.build/${six}BothHaps_aggressivecorrection.FINAL.fa --par 48 --o /workdir/mcs368/panand_assemblies/hic/${six}/${six}BothHaps/04.build/${six}BothHaps_aggressivecorrection_TRASH


## trash2 on atlas
## on atlas - not having error?!?!?!?!
cd /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats
six=irugos
mkdir -p ${six}BothHaps_aggressivecorrection_TRASH2
conda activate trash ## i guess do it before???
sbatch -A buckler_lab_panand -p atlas --ntasks-per-node=48 --time=10-00:00 --wrap="six=irugos; conda activate trash; Rscript /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/TRASH_2/src/TRASH.R -f /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/${six}BothHaps_aggressivecorrection.FINAL.fa -o /project/buckler_lab_panand/michelle.stitzer/panand_assemblies/hic/repeats/${six}BothHaps_aggressivecorrection_TRASH2 -p 46"


## run tidk
conda activate tidk
tidk find --clade Poales --output irugosBothHaps_aggressivecorrection.tidk --dir . irugosBothHaps_aggressivecorrection.FINAL.fa 

## run helixer
ssh cbsugpu08
cd /workdir/mcs368/
nvidia-smi ## check if anybody's using gpu
scp cbsuxm01:/workdir/mcs368/panand_assemblies/hic/irugos/irugosBothHaps/04.build/irugosBothHaps_aggressivecorrection.FINAL.fa .
singularity run --nv --bind $PWD --pwd $PWD /programs/helixer-0.3.5/helixer.sif Helixer.py --fasta-path irugosBothHaps_aggressivecorrection.FINAL.fa --lineage land_plant --gff-output-path irugosBothHaps_aggressivecorrection.FINAL.helixer.gff3
scp irugosBothHaps_aggressivecorrection.FINAL.helixer.gff3 cbsuxm01:/workdir/mcs368/panand_assemblies/hic/irugos/irugosBothHaps/04.build/irugosBothHaps_aggressivecorrection.FINAL.helixer.gff3


### gc content
samtools faidx irugosBothHaps_aggressivecorrection.FINAL.fa
bedtools makewindows -g irugosBothHaps_aggressivecorrection.FINAL.fa.fai -w 1000000 > irugosBothHaps_aggressivecorrection.FINAL.1mbwindows.bed
bedtools nuc -fi irugosBothHaps_aggressivecorrection.FINAL.fa -bed irugosBothHaps_aggressivecorrection.FINAL.1mbwindows.bed > irugosBothHaps_aggressivecorrection.FINAL.1mbnuccontent.bed

## dotplot
./fa_to_dotplot.sh irugos/irugosBothHaps/04.build/irugosBothHaps_aggressivecorrection.FINAL.fa irugosBothHaps_aggressivecorrection 2















#########################
# scaffold each haplotype separately

conda activate haphic

bwa index irugos.asm.hic.hap1.p_ctg.fa
bwa mem -5SP -t 128 irugos.asm.hic.hap1.p_ctg.fa AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o irugosHap1.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam irugosHap1.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o irugosHap1.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline irugos.asm.hic.hap1.p_ctg.fa irugosHap1.filtered.bam 9 --outdir irugosHap1


### sort by reference https://github.com/zengxiaofei/HapHiC/issues/92
cd ${six}/${six}Hap1/04.build/
# The preset can be `asm5` if the reference genome is well-assembled from the same species
/programs/minimap2-2.30/minimap2  -x asm20 /workdir/mcs368/panand_assemblies/genomes/Pvaginatum_672_v3.0.fa ../../irugos.asm.hic.hap1.p_ctg.fa --secondary=no -t 96 -o asm_to_ref.paf

../../../HapHiC/haphic refsort scaffolds.raw.agp asm_to_ref.paf > scaffolds.refsort.agp

ln -s ../../irugos.asm.hic.hap1.p_ctg.fa .
samtools faidx irugos.asm.hic.hap1.p_ctg.fa
# This command will generate several files with the prefix 'out_JBAT'. Once out_JBAT.assembly has been fully generated, you can terminate the command (there's no need to wait for the generation of out_JBAT.txt, which can be slow)
/local/workdir/mcs368/panand_assemblies/hic/HapHiC/scripts/../utils/juicer pre -a -q 1 -o out_JBAT /local/workdir/mcs368/panand_assemblies/hic/irugos/irugosHap1.filtered.bam scaffolds.refsort.agp irugos.asm.hic.hap1.p_ctg.fa.fai >out_JBAT.log 2>&1


# By default, scaffolds are output based on the alphabetical order of the chromosome IDs of the reference genome
$ haphic refsort 04.build/scaffolds.raw.agp asm_to_ref.paf > scaffolds.refsort.agp
# You can specify the order by listing chromosome IDs of the reference genome separated by commas (no spaces)
$ haphic refsort 04.build/scaffolds.raw.agp asm_to_ref.paf --ref_order "chr1,chr2,chr3,chr4,..." > scaffolds.refsort.agp
# If you want to generate a new FASTA file (default name: `scaffolds.refsort.fa`) as well
$ haphic refsort 04.build/scaffolds.raw.agp asm_to_ref.paf --fasta asm.fa > scaffolds.refsort.agp
# If you want to run `haphic refsort` on manually curated `out_JBAT.FINAL.agp`
$ haphic refsort out_JBAT.FINAL.agp asm_to_ref.paf > scaffolds.refsort.agp




## hap2
bwa index irugos.asm.hic.hap2.p_ctg.fa
bwa mem -5SP -t 128 irugos.asm.hic.hap2.p_ctg.fa AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o irugosHap2.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam irugosHap2.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o irugosHap2.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline irugos.asm.hic.hap2.p_ctg.fa irugosHap2.filtered.bam 10 --outdir irugosHap2

### p_ctg

bwa index irugos.asm.hic.p_ctg.fa
bwa mem -5SP -t 128 irugos.asm.hic.p_ctg.fa AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o irugos_p_ctg.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam irugos_p_ctg.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o irugos_p_ctg.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline irugos.asm.hic.p_ctg.fa irugos_p_ctg.filtered.bam 10 --outdir irugos_p_ctg

cd irugos_p_ctg/04.build/
/programs/bbmap-38.96/stats.sh -Xmx200g in="scaffolds.fa" out="scaffolds.stats.txt" overwrite='true'

bash juicebox.sh
cd ../../../
./fa_to_dotplot.sh irugos/irugos_p_ctg/04.build/scaffolds.fa irugos_p_ctgHapHiC 1



## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' irugos.asm.hic.hap1.p_ctg.fa > irugos.asm.hic.hap1.p_ctg.tmp.fasta
barrnap --quiet --kingdom euk irugos.asm.hic.hap1.p_ctg.tmp.fasta > irugos.asm.hic.hap1.p_ctg.rrna.gff3
rm irugos.asm.hic.hap1.p_ctg.tmp.fasta
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' irugos.asm.hic.hap2.p_ctg.fa > irugos.asm.hic.hap2.p_ctg.tmp.fasta
barrnap --quiet --kingdom euk irugos.asm.hic.hap2.p_ctg.tmp.fasta > irugos.asm.hic.hap2.p_ctg.rrna.gff3
rm irugos.asm.hic.hap2.p_ctg.tmp.fasta
