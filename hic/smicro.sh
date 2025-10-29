## this should be relatively easy - 0.4% het, no Bionano
conda activate panand_assemblies

## check kmers first!
mkdir tmp
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 smicro_hifi.fastq.gz smicrokmc tmp/
/programs/kmc-3.2.4/kmc_tools transform smicrokmc histogram smicro.histo -cx10000

../genomescope2.0/genomescope.R -i smicro.histo -o smicro_genomescope -k 21


## dual-scaf, hic assembly
 hifiasm -o smicro.asm -t64 --dual-scaf --h1 AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz --h2 AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz smicro_hifi.fastq.gz

awk '/^S/{print ">"$2"\n"$3}' smicro.asm.hic.hap1.p_ctg.gfa > smicro.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' smicro.asm.hic.hap2.p_ctg.gfa > smicro.asm.hic.hap2.p_ctg.fa


## check assembly, how many contigs per chromosome??
/programs/bbmap-38.96/stats.sh -Xmx200g in="smicro.asm.hic.hap1.p_ctg.fa" out="smicro.asm.hic.hap1.p_ctg.stats.txt"
/programs/bbmap-38.96/stats.sh -Xmx200g in="smicro.asm.hic.hap2.p_ctg.fa" out="smicro.asm.hic.hap2.p_ctg.stats.txt"
## hap1 925 Mb, 882 scaffolds, N50 82 Mb 5 scaffold N50
## hap2 899 Mb, 158 scaffolds, N50 83 Mb 5 scaffold N50

## check out a dotplot
cd ..
./fa_to_dotplot.sh smicro/smicro.asm.hic.hap1.p_ctg.fa smicroHap1 1
./fa_to_dotplot.sh smicro/smicro.asm.hic.hap2.p_ctg.fa smicroHap2 1
cd smicro


## combine haplotypes and haphic



cat smicro.asm.hic.hap1.p_ctg.fa smicro.asm.hic.hap2.p_ctg.fa >smicro.asm.hic.bothhaps.p_ctg.fa

bwa index smicro.asm.hic.bothhaps.p_ctg.fa
bwa mem -5SP -t 128 smicro.asm.hic.bothhaps.p_ctg.fa AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o smicroBothHaps.bam
## duh it can't index if not sorted
##samtools index smicroBothHaps.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam smicroBothHaps.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o smicroBothHaps.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
#../HapHiC/haphic pipeline smicro.asm.hic.bothhaps.p_ctg.fa smicroBothHaps.filtered.bam 20 --gfa "smicro.asm.hic.hap1.p_ctg.gfa,smicro.asm.hic.hap2.p_ctg.gfa" --outdir smicroBothHaps

## skip filtering? potential rec here https://github.com/zengxiaofei/HapHiC/issues/86
## use gfa to be sure not to combine haps
../HapHiC/haphic pipeline smicro.asm.hic.bothhaps.p_ctg.fa smicroBothHaps.bam 20 --gfa "smicro.asm.hic.hap1.p_ctg.gfa,smicro.asm.hic.hap2.p_ctg.gfa" --outdir smicroBothHapsUnfiltered



## to correct contigs with hic info
## --correct_nrounds 2 
## to remove hic links between alleles, so set to number of haplotypes present (polyploid)
## --remove_allelic_links 2

## try quickview 
#../HapHiC/haphic pipeline smicro.asm.hic.bothhaps.p_ctg.fa smicroBothHaps.filtered.bam 20 --quick_view --gfa "smicro.asm.hic.hap1.p_ctg.gfa,smicro.asm.hic.hap2.p_ctg.gfa" --outdir smicroBothHapsQuickView
../HapHiC/haphic pipeline smicro.asm.hic.bothhaps.p_ctg.fa smicroBothHaps.bam 20 --quick_view  --gfa "smicro.asm.hic.hap1.p_ctg.gfa,smicro.asm.hic.hap2.p_ctg.gfa" --outdir smicroBothHapsQuickViewUnfiltered


#########################
## scaffold each haplotype separately

conda activate haphic

bwa index smicro.asm.hic.hap1.p_ctg.fa
bwa mem -5SP -t 128 smicro.asm.hic.hap1.p_ctg.fa AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o smicroHap1.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam smicroHap1.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o smicroHap1.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline smicro.asm.hic.hap1.p_ctg.fa smicroHap1.filtered.bam 10 --outdir smicroHap1

## hap2
bwa index smicro.asm.hic.hap2.p_ctg.fa
bwa mem -5SP -t 128 smicro.asm.hic.hap2.p_ctg.fa AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o smicroHap2.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam smicroHap2.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o smicroHap2.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline smicro.asm.hic.hap2.p_ctg.fa smicroHap2.filtered.bam 10 --outdir smicroHap2

## where are the repeats on the genomes??
## ugh telomeres mess it up so make temp
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' smicro.asm.hic.hap1.p_ctg.fa > smicro.asm.hic.hap1.p_ctg.tmp.fasta
barrnap --quiet --kingdom euk smicro.asm.hic.hap1.p_ctg.tmp.fasta > smicro.asm.hic.hap1.p_ctg.rrna.gff3
rm smicro.asm.hic.hap1.p_ctg.tmp.fasta
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' smicro.asm.hic.hap2.p_ctg.fa > smicro.asm.hic.hap2.p_ctg.tmp.fasta
barrnap --quiet --kingdom euk smicro.asm.hic.hap2.p_ctg.tmp.fasta > smicro.asm.hic.hap2.p_ctg.rrna.gff3
rm smicro.asm.hic.hap2.p_ctg.tmp.fasta


### then use corrected contact maps
# Generate the final FASTA file for the scaffolds
## from the 04 directory...
python ../../../juicebox_scripts/juicebox_scripts/juicebox_assembly_converter.py -a out_JBAT.review.smicroHap2.assembly -f scaffolds.fa -s


