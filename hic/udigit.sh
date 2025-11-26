## very low het, Bionano, a little nanopore 4.2% ks btwn dups, 5.4Gb
conda activate panand_assemblies

## check kmers first!
mkdir tmp
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 udigit_hifi.fastq.gz udigitkmc tmp/
/programs/kmc-3.2.4/kmc_tools transform udigitkmc histogram udigit.histo -cx10000

../genomescope2.0/genomescope.R -i udigit.histo -o udigit_genomescope -k 21

## i may just want to use the p_ctg on this one!!!
hifiasm -o udigit.asm --dual-scaf -t64 --ul /workdir/mcs368/panand_assemblies/scaffolding/udigit/udigit.50000.fq.gz --h1 AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz --h2 AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz udigit_hifi.fastq

awk '/^S/{print ">"$2"\n"$3}' udigit.asm.hic.hap1.p_ctg.gfa > udigit.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' udigit.asm.hic.hap2.p_ctg.gfa > udigit.asm.hic.hap2.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' udigit.asm.hic.p_ctg.gfa > udigit.asm.hic.p_ctg.fa

## check assembly, how many contigs per chromosome??
/programs/bbmap-38.96/stats.sh -Xmx200g in="udigit.asm.hic.hap1.p_ctg.fa" out="irugos.asm.hic.hap1.p_ctg.stats.txt" overwrite='true'
/programs/bbmap-38.96/stats.sh -Xmx200g in="irugos.asm.hic.hap2.p_ctg.fa" out="irugos.asm.hic.hap2.p_ctg.stats.txt" overwrite='true'
## hap1 619 Mb, 590 scaffolds, N50 3.6 Mb 45 scaffold N50
/programs/bbmap-38.96/stats.sh -Xmx200g in="irugos.asm.hic.p_ctg.fa" out="irugos.asm.hic.p_ctg.stats.txt" overwrite='true'


cd ../
./fa_to_dotplot.sh udigit/udigitDip.asm.hic.hap1.p_ctg.fa udigitDipHap1 3
./fa_to_dotplot.sh udigit/udigitDip.asm.hic.hap2.p_ctg.fa udigitDipHap2 3
cd udigit
