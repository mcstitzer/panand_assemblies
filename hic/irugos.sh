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



awk '/^S/{print ">"$2"\n"$3}' irugos.asm.hic.hap1.p_ctg.gfa > irugos.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' irugos.asm.hic.hap2.p_ctg.gfa > irugos.asm.hic.hap2.p_ctg.fa
## what about unphased contigs??
awk '/^S/{print ">"$2"\n"$3}' irugos.asm.hic.p_ctg.gfa > irugos.asm.hic.p_ctg.fa

## check assembly, how many contigs per chromosome??
/programs/bbmap-38.96/stats.sh -Xmx200g in="irugos.asm.hic.hap1.p_ctg.fa" out="irugos.asm.hic.hap1.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="irugos.asm.hic.hap2.p_ctg.fa" out="irugos.asm.hic.hap2.p_ctg.stats.txt" overwrite='true'
## hap1 619 Mb, 590 scaffolds, N50 3.6 Mb 45 scaffold N50
/programs/bbmap-38.96/stats.sh -Xmx200g in="irugos.asm.hic.p_ctg.fa" out="irugos.asm.hic.p_ctg.stats.txt" overwrite='true'

## check out a dotplot
cd ..
./fa_to_dotplot.sh irugos/irugos.asm.hic.hap1.p_ctg.fa irugosHap1 1
./fa_to_dotplot.sh irugos/irugos.asm.hic.hap2.p_ctg.fa irugosHap2 1
./fa_to_dotplot.sh irugos/irugos.asm.hic.p_ctg.fa irugos_p_ctg 1

cd irugos



#########################
## scaffold each haplotype separately

conda activate haphic

bwa index irugos.asm.hic.hap1.p_ctg.fa
bwa mem -5SP -t 128 irugos.asm.hic.hap1.p_ctg.fa AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 64 -S -h -b -F 3340 -o irugosHap1.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam irugosHap1.bam 1 --nm 3 --threads 64 | samtools view - -b -@ 64 -o irugosHap1.filtered.bam

## filtering reduces 27 Gb bam to 165 mb =--- seems dangerous! I think it's because this is so low heterozygosity.

## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline irugos.asm.hic.hap1.p_ctg.fa irugosHap1.filtered.bam 10 --outdir irugosHap1

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
