
samtools bam2fq -@ 64 /workdir/mcs368/panand_snps/input/smicro_reads_withQ.bam > smicro_hifi.fastq
hifiasm -o smicro.asm -t64 --h1 AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz --h2 AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz smicro_hifi.fastq

awk '/^S/{print ">"$2"\n"$3}' smicro.asm.hic.hap1.p_ctg.gfa > smicro.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' smicro.asm.hic.hap2.p_ctg.gfa > smicro.asm.hic.hap2.p_ctg.fa
./fa_to_dotplot.sh smicro/smicro.asm.hic.hap1.p_ctg.fa smicroHap1 1
./fa_to_dotplot.sh smicro/smicro.asm.hic.hap2.p_ctg.fa smicroHap2 1
## dotplots in ~/transfer/smicroHap1_dotplot.pdf
## remap hic to assemblies individually
bwa index smicro.asm.hic.hap1.p_ctg.fa
## 5SP deals with weird hic distribution of reads
bwa mem -t 64 -5SP smicro.asm.hic.hap1.p_ctg.fa AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o smicroHap1_hic.sorted.bam 
samtools index smicroHap1_hic.sorted.bam
bwa index smicro.asm.hic.hap2.p_ctg.fa
bwa mem -t 64 -5SP smicro.asm.hic.hap2.p_ctg.fa AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o smicroHap2_hic.sorted.bam 
samtools index smicroHap2_hic.sorted.bam

## aah to screen, so this is only hap2
# samblaster:
# samblaster: Pair Type        Type_ID_Count   %Type/All_IDs Dup_ID_Count  %Dups/Type_ID_Count  %Dups/All_Dups  %Dups/All_IDs
# samblaster: ---------------------------------------------------------------------------------------------------------------
# samblaster: Both Unmapped          797657         0.450              0           0.000             0.000          0.000
# samblaster: Orphan/Singleton      2133432         1.204         240547          11.275             0.325          0.136
# samblaster: Both Mapped         174210367        98.345       73733130          42.324            99.675         41.624
# samblaster: Total               177141456       100.000       73973677          41.760           100.000         41.760
# samblaster:
# samblaster: Marked    73973677 of  177141456 (41.760%) total read ids as duplicates using 2016188k memory in 4M52S(291.745S) CPU seconds and 51M27S(3087S) wall time.

## then i need to have figured out my pairtools pipeline!!!
samtools view -h smicroHap1_hic.sorted.bam | PretextMap -o smicroHap1.pretext --sortby length --sortorder descend --mapq 0
samtools view -h smicroHap2_hic.sorted.bam | PretextMap -o smicroHap2.pretext --sortby length --sortorder descend --mapq 0
 
yahs smicro.asm.hic.hap1.p_ctg.fa smicroHap1_hic.sorted.bam -o smicro -e GATC
yahs smicro.asm.hic.hap2.p_ctg.fa smicroHap2_hic.sorted.bam -o smicroHap2 -e GATC
./fa_to_dotplot.sh smicro/smicro_scaffolds_final.fa smicroHap1yahs 1
./fa_to_dotplot.sh smicro/smicroHap2_scaffolds_final.fa smicroHap1yahs 2

juicer pre smicroHap2.bin smicroHap2_scaffolds_final.agp smicro.asm.hic.hap2.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part
mv alignments_sorted.txt.part alignments_sorted.txt
(java -jar -Xmx32G ../juicer_tools.2.18.00.jar  pre alignments_sorted.txt out.hic.part smicroHap2_scaffolds_final.fa.fai ) && (mv out.hic.part out.hic)

juicer pre -a -o smicroHap1out_JBAT smicroHap1_hic.sorted.bam smicro_scaffolds_final.agp smicro.asm.hic.hap1.p_ctg.fa.fai >out_JBAT.log 2>&1
 
 #### ALSO TRY WITH DUAL-SCAF to scaffold haplotypes off of each other
 hifiasm -o smicroDS.asm -t64 --dual-scaf --h1 AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz --h2 AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz smicro_hifi.fastq

 
 
  
 ## let's try udigit??
samtools bam2fq -@ 64 /workdir/mcs368/panand_snps/input/udigit_reads_withQ.bam > udigit_hifi.fastq
## i have ultralong nanopore reads, hic, and hifi, and will have it output 6 haplotypes
hifiasm -o udigit.asm -t64 --n-hap 6 --ul /workdir/mcs368/panand_assemblies/scaffolding/udigit/udigit.50000.fq.gz --h1 AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz --h2 AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz udigit_hifi.fastq

### not sure this is deep enough hifi - it's only 17x (but of genome size, not haploid size)

awk '/^S/{print ">"$2"\n"$3}' udigit.asm.hic.hap1.p_ctg.gfa > udigit.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' udigit.asm.hic.hap2.p_ctg.gfa > udigit.asm.hic.hap2.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' udigit.asm.hic.hap3.p_ctg.gfa > udigit.asm.hic.hap3.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' udigit.asm.hic.hap4.p_ctg.gfa > udigit.asm.hic.hap4.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' udigit.asm.hic.hap5.p_ctg.gfa > udigit.asm.hic.hap5.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' udigit.asm.hic.hap6.p_ctg.gfa > udigit.asm.hic.hap6.p_ctg.fa
./fa_to_dotplot.sh udigit/udigit.asm.hic.hap1.p_ctg.fa udigitHap1 1
./fa_to_dotplot.sh udigit/udigit.asm.hic.hap2.p_ctg.fa udigitHap2 1
./fa_to_dotplot.sh udigit/udigit.asm.hic.hap3.p_ctg.fa udigitHap3 1
./fa_to_dotplot.sh udigit/udigit.asm.hic.hap4.p_ctg.fa udigitHap4 1
./fa_to_dotplot.sh udigit/udigit.asm.hic.hap5.p_ctg.fa udigitHap5 1
./fa_to_dotplot.sh udigit/udigit.asm.hic.hap6.p_ctg.fa udigitHap6 1






## also try assembling as a "diploid" becasue subgenomes are different enough??
hifiasm -o udigitDip.asm -t64 --ul /workdir/mcs368/panand_assemblies/scaffolding/udigit/udigit.50000.fq.gz --h1 AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz --h2 AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz udigit_hifi.fastq

awk '/^S/{print ">"$2"\n"$3}' udigitDip.asm.hic.hap1.p_ctg.gfa > udigitDip.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' udigitDip.asm.hic.hap2.p_ctg.gfa > udigitDip.asm.hic.hap2.p_ctg.fa
./fa_to_dotplot.sh udigit/udigitDip.asm.hic.hap1.p_ctg.fa udigitDipHap1 3
./fa_to_dotplot.sh udigit/udigitDip.asm.hic.hap2.p_ctg.fa udigitDipHap2 3

bwa index udigitDip.asm.hic.hap1.p_ctg.fa
bwa index udigitDip.asm.hic.hap2.p_ctg.fa
bwa mem -t 64 -5SP udigitDip.asm.hic.hap1.p_ctg.fa AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o udigitDipHap1_hic.sorted.bam 
samtools index udigitDipHap1_hic.sorted.bam
bwa mem -t 64 -5SP udigitDip.asm.hic.hap2.p_ctg.fa AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o udigitDipHap2_hic.sorted.bam 
samtools index udigitDipHap2_hic.sorted.bam
samtools view -h udigitDipHap1_hic.sorted.bam | PretextMap -o udigitDipHap1.pretext --sortby length --sortorder descend --mapq 0
samtools view -h udigitDipHap2_hic.sorted.bam | PretextMap -o udigitDipHap2.pretext --sortby length --sortorder descend --mapq 0
 
 samtools faidx udigitDip.asm.hic.hap1.p_ctg.fa
 samtools faidx udigitDip.asm.hic.hap2.p_ctg.fa
yahs udigitDip.asm.hic.hap1.p_ctg.fa udigitDipHap1_hic.sorted.bam -o udigitDipHap1 -e GATC
yahs udigitDip.asm.hic.hap2.p_ctg.fa udigitDipHap2_hic.sorted.bam -o udigitDipHap2 -e GATC
../fa_to_dotplot.sh udigit/udigitDipHap1_scaffolds_final.fa udigitDipHap1yahs 3
../fa_to_dotplot.sh udigit/udigitDipHap2_scaffolds_final.fa udigitDipHap1yahs 3



juicer pre udigitDipHap1.bin udigitDipHap1_scaffolds_final.agp udigitDip.asm.hic.hap1.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > udigitDipHap1.alignments_sorted.txt
java -jar -Xmx32G ../juicer_tools.2.18.00.jar  pre udigitDipHap1.alignments_sorted.txt udigitDipHap1.hic <( cut -f1,2 udigitDipHap1_scaffolds_final.fa.fai )

 

## try scoparium??
### just using the hifi data, sorry clr :(
## scp them over, then assemble

scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240003_Schizachyrium_scoparium_HiC_pairsR1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240003_Schizachyrium_scoparium_HiC_pairsR2.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082625_PanAnd/AV9590001_Schizachyrium_scoparium_20250822_01_R84064-1_A01.filter.bam .

samtools bam2fq -@ 64 AV9590001_Schizachyrium_scoparium_20250822_01_R84064-1_A01.filter.bam > sscopa_hifi.fastq

hifiasm -o sscopa4hap.asm -t64 --n-hap 4 --h1 AV9240003_Schizachyrium_scoparium_HiC_pairsR1.fq.gz --h2 AV9240003_Schizachyrium_scoparium_HiC_pairsR2.fq.gz sscopa_hifi.fastq
hifiasm -o sscopa2hap.asm -t64 --n-hap 2 --h1 AV9240003_Schizachyrium_scoparium_HiC_pairsR1.fq.gz --h2 AV9240003_Schizachyrium_scoparium_HiC_pairsR2.fq.gz sscopa_hifi.fastq

awk '/^S/{print ">"$2"\n"$3}' sscopa2hap.asm.hic.hap1.p_ctg.gfa > sscopa2hap.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' sscopa2hap.asm.hic.hap2.p_ctg.gfa > sscopa2hap.asm.hic.hap2.p_ctg.fa
./fa_to_dotplot.sh sscopa/sscopa2hap.asm.hic.hap1.p_ctg.fa sscopa2hapHap1 2
./fa_to_dotplot.sh sscopa/sscopa2hap.asm.hic.hap2.p_ctg.fa sscopa2hapHap2 2

bwa index sscopa2hap.asm.hic.hap1.p_ctg.fa
## 5SP deals with weird hic distribution of reads
bwa mem -t 64 -5SP sscopa2hap.asm.hic.hap1.p_ctg.fa AV9240003_Schizachyrium_scoparium_HiC_pairsR1.fq.gz AV9240003_Schizachyrium_scoparium_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o sscopaHap1_hic.sorted.bam 
samtools index sscopaHap1_hic.sorted.bam
bwa index sscopa2hap.asm.hic.hap2.p_ctg.fa
bwa mem -t 64 -5SP sscopa2hap.asm.hic.hap2.p_ctg.fa AV9240003_Schizachyrium_scoparium_HiC_pairsR1.fq.gz AV9240003_Schizachyrium_scoparium_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o sscopaHap2_hic.sorted.bam 
samtools index sscopaHap2_hic.sorted.bam

 samtools faidx sscopa2hap.asm.hic.hap1.p_ctg.fa
 samtools faidx sscopa2hap.asm.hic.hap2.p_ctg.fa
yahs sscopa2hap.asm.hic.hap1.p_ctg.fa sscopaHap1_hic.sorted.bam -o sscopaHap1 -e GATC
yahs sscopa2hap.asm.hic.hap2.p_ctg.fa sscopaHap2_hic.sorted.bam -o sscopaHap2 -e GATC
./fa_to_dotplot.sh sscopa/sscopaHap1_scaffolds_final.fa sscopaHap1yahs 3
./fa_to_dotplot.sh sscopa/sscopaHap2_scaffolds_final.fa sscopaHap1yahs 3



juicer pre sscopaHap1.bin sscopaHap1_scaffolds_final.agp sscopa2hap.asm.hic.hap1.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > sscopaHap1.alignments_sorted.txt
java -jar -Xmx32G ../juicer_tools.2.18.00.jar  pre sscopaHap1.alignments_sorted.txt sscopaHap1.hic <( cut -f1,2 sscopaHap1_scaffolds_final.fa.fai )


# awk '/^S/{print ">"$2"\n"$3}' sscopa4hap.asm.hic.hap1.p_ctg.gfa > sscopa4hap.asm.hic.hap1.p_ctg.fa
# awk '/^S/{print ">"$2"\n"$3}' sscopa4hap.asm.hic.hap2.p_ctg.gfa > sscopa4hap.asm.hic.hap2.p_ctg.fa
# awk '/^S/{print ">"$2"\n"$3}' sscopa4hap.asm.hic.hap1.p_ctg.gfa > sscopa4hap.asm.hic.hap3.p_ctg.fa
# awk '/^S/{print ">"$2"\n"$3}' sscopa4hap.asm.hic.hap2.p_ctg.gfa > sscopa4hap.asm.hic.hap4.p_ctg.fa
# ./fa_to_dotplot.sh sscopa/sscopa4hap.asm.hic.hap1.p_ctg.fa sscopa4hapHap1 2
# ./fa_to_dotplot.sh sscopa/sscopa4hap.asm.hic.hap2.p_ctg.fa sscopa4hapHap2 2
# ./fa_to_dotplot.sh sscopa/sscopa4hap.asm.hic.hap3.p_ctg.fa sscopa4hapHap3 2
# ./fa_to_dotplot.sh sscopa/sscopa4hap.asm.hic.hap4.p_ctg.fa sscopa4hapHap4 2


## etrips

## big so in two parts?
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240008_Elionurus_tripsacoides_HiC_pairsR1_PART1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240008_Elionurus_tripsacoides_HiC_pairsR1_PART2.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240008_Elionurus_tripsacoides_HiC_pairsR2_PART1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240008_Elionurus_tripsacoides_HiC_pairsR2_PART2.fq.gz .

cat AV9240008_Elionurus_tripsacoides_HiC_pairsR1_PART1.fq.gz AV9240008_Elionurus_tripsacoides_HiC_pairsR1_PART2.fq.gz > etrips.hic.read1.fq.gz
cat AV9240008_Elionurus_tripsacoides_HiC_pairsR2_PART1.fq.gz AV9240008_Elionurus_tripsacoides_HiC_pairsR2_PART2.fq.gz > etrips.hic.read2.fq.gz

rm AV9240008_Elionurus_tripsacoides_HiC_pairsR1_PART1.fq.gz
rm AV9240008_Elionurus_tripsacoides_HiC_pairsR1_PART2.fq.gz
rm AV9240008_Elionurus_tripsacoides_HiC_pairsR2_PART1.fq.gz
rm AV9240008_Elionurus_tripsacoides_HiC_pairsR2_PART2.fq.gz

samtools bam2fq -@ 64 /workdir/mcs368/panand_snps/input/etrips_reads_withQ.bam > etrips_hifi.fastq

## oh great! Cinta has lots of nanopore for this
scp cbsublfs1:/data1/PanAnd1/RawSeqData/nanopore/andropogoneae/MinION/A1021-0002_02-recall.5.0.16/A1021-0002_02-fastq.tgz .
scp cbsublfs1:/data1/PanAnd1/RawSeqData/nanopore/andropogoneae/MinION/A1021-0002_03-recall.5.0.16/A1021-0002_03-fastq.tgz .
scp cbsublfs1:/data1/PanAnd1/RawSeqData/nanopore/andropogoneae/MinION/A1021-0002_04-recall.5.0.16/A1021-0002_04-fastq.tgz .
scp cbsublfs1:/data1/PanAnd1/RawSeqData/nanopore/andropogoneae/MinION/A1021_01-recall.5.0.16/A1021_01-fastq.tgz .
scp cbsublfs1:/data1/PanAnd1/RawSeqData/nanopore/andropogoneae/MinION/A1021_05-recall.5.0.16/A1021_05-fastq.tgz .
scp cbsublfs1:/data1/PanAnd1/RawSeqData/nanopore/andropogoneae/MinION/A1021_06-recall.5.0.16/A1021_06-fastq.tgz .
scp cbsublfs1:/data1/PanAnd1/RawSeqData/nanopore/andropogoneae/MinION/A1021_07-recall.5.0.16/A1021_07-fastq.tgz .
scp cbsublfs1:/data1/PanAnd1/RawSeqData/nanopore/andropogoneae/MinION/A1021_08-recall.5.0.16/A1021_08-fastq.tgz .



# Base remote path
REMOTE_BASE="cbsublfs1:/data1/PanAnd1/RawSeqData/nanopore/andropogoneae/MinION"
# List of dataset directories and .tgz filenames (without full path)
FILES=(
    "A1021-0002_02-recall.5.0.16/A1021-0002_02-fastq.tgz"
    "A1021-0002_03-recall.5.0.16/A1021-0002_03-fastq.tgz"
    "A1021-0002_04-recall.5.0.16/A1021-0002_04-fastq.tgz"
    "A1021_01-recall.5.0.16/A1021_01-fastq.tgz"
    "A1021_05-recall.5.0.16/A1021_05-fastq.tgz"
    "A1021_06-recall.5.0.16/A1021_06-fastq.tgz"
    "A1021_07-recall.5.0.16/A1021_07-fastq.tgz"
    "A1021_08-recall.5.0.16/A1021_08-fastq.tgz"
)
# Loop through each file
for FILE in "${FILES[@]}"; do
    echo "Processing $FILE..."
    # scp from remote
    scp "${REMOTE_BASE}/${FILE}" .
    # Get the local filename, dataset name (strip .tgz)
    TGZ=$(basename "$FILE")
    SAMPLE=${TGZ%.tgz}
    # Extract
    tar -xzvf "$TGZ"
    # Run seqkit (assuming fastq files are in SAMPLE/basecalls/fastq_pass/)
    /programs/seqkit-0.15.0/seqkit seq -m 50000 "${SAMPLE}/pass/"*.fastq > "${SAMPLE}.50000.fq"
    # Run bbmap stats
    /programs/bbmap-38.96/stats.sh -Xmx200g in="${SAMPLE}.50000.fq" out="${SAMPLE}.50000.stats.txt"
    # Delete original tgz
    rm "$TGZ"
    echo "Done with $SAMPLE"
    echo "----------------------"
done

cat *.50000.fq > etrips.50000.fq
bgzip etrips.50000.fq
/programs/bbmap-38.96/stats.sh -Xmx200g in="etrips.50000.fq.gz" out="etrips.50000.stats.txt"

hifiasm -o etrips2hap.asm -t64 --n-hap 2 --ul etrips.50000.fq.gz --h1 etrips.hic.read1.fq.gz --h2 etrips.hic.read2.fq.gz etrips_hifi.fastq

awk '/^S/{print ">"$2"\n"$3}' etrips2hap.asm.hic.hap1.p_ctg.gfa > etrips2hap.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' etrips2hap.asm.hic.hap2.p_ctg.gfa > etrips2hap.asm.hic.hap2.p_ctg.fa
./fa_to_dotplot.sh etrips/etrips2hap.asm.hic.hap1.p_ctg.fa etrips2hapHap1 2
./fa_to_dotplot.sh etrips/etrips2hap.asm.hic.hap2.p_ctg.fa etrips2hapHap2 2

bwa index etrips2hap.asm.hic.hap1.p_ctg.fa
## 5SP deals with weird hic distribution of reads
bwa mem -t 64 -5SP etrips2hap.asm.hic.hap1.p_ctg.fa etrips.hic.read1.fq.gz etrips.hic.read2.fq.gz | samblaster | samtools sort -@ 16 -o etrips2hapHap1_hic.sorted.bam 
samtools index etrips2hapHap1.sorted.bam
bwa index etrips2hap.asm.hic.hap2.p_ctg.fa
bwa mem -t 64 -5SP etrips2hap.asm.hic.hap2.p_ctg.fa etrips.hic.read1.fq.gz etrips.hic.read2.fq.gz | samblaster | samtools sort -@ 16 -o etrips2hapHap2_hic.sorted.bam 
samtools index etrips2hapHap2_hic.sorted.bam

 samtools faidx etrips2hap.asm.hic.hap1.p_ctg.fa
 samtools faidx etrips2hap.asm.hic.hap2.p_ctg.fa
yahs etrips2hap.asm.hic.hap1.p_ctg.fa etrips2hapHap1_hic.sorted.bam -o etrips2hapHap1 -e GATC
yahs etrips2hap.asm.hic.hap2.p_ctg.fa etrips2hapHap2_hic.sorted.bam -o etrips2hapHap2 -e GATC
./fa_to_dotplot.sh etrips/etrips2hapHap1_scaffolds_final.fa etrips2hapHap1yahs 3
./fa_to_dotplot.sh etrips/etrips2hapHap2_scaffolds_final.fa etrips2hapHap2yahs 3



juicer pre etrips2hapHap1.bin etrips2hapHap1_scaffolds_final.agp etrips2hap.asm.hic.hap1.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > etrips2hapHap1.alignments_sorted.txt
samtools faidx etrips2hapHap1_scaffolds_final.fa
java -jar -Xmx32G ../juicer_tools.2.18.00.jar  pre etrips2hapHap1.alignments_sorted.txt etrips2hapHap1.hic <( cut -f1,2 etrips2hapHap1_scaffolds_final.fa.fai )




##### irugos


## big so in two parts?
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz .

samtools bam2fq -@ 64 /workdir/mcs368/panand_snps/input/irugos_reads_withQ.bam  > irugos_hifi.fastq

hifiasm -o irugos.asm -t64 --n-hap 2 --dual-scaf --h1 AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz --h2 AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz irugos_hifi.fastq

awk '/^S/{print ">"$2"\n"$3}' irugos.asm.hic.hap1.p_ctg.gfa > irugos.asm.hic.hap1.p_ctg.fa
awk '/^S/{print ">"$2"\n"$3}' irugos.asm.hic.hap2.p_ctg.gfa > irugos.asm.hic.hap2.p_ctg.fa
./fa_to_dotplot.sh irugos/irugos.asm.hic.hap1.p_ctg.fa irugosHap1 2
./fa_to_dotplot.sh irugos/irugos.asm.hic.hap2.p_ctg.fa irugosHap2 2

bwa index irugos.asm.hic.hap1.p_ctg.fa
## 5SP deals with weird hic distribution of reads
bwa mem -t 64 -5SP irugos.asm.hic.hap1.p_ctg.fa AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o irugosHap1_hic.sorted.bam 
samtools index irugosHap1.sorted.bam
bwa index irugos.asm.hic.hap2.p_ctg.fa
bwa mem -t 64 -5SP irugos.asm.hic.hap2.p_ctg.fa AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o irugosHap2_hic.sorted.bam 
samtools index irugosHap2_hic.sorted.bam

 samtools faidx irugos.asm.hic.hap1.p_ctg.fa
 samtools faidx irugos.asm.hic.hap2.p_ctg.fa
yahs irugos.asm.hic.hap1.p_ctg.fa irugosHap1_hic.sorted.bam -o irugosHap1 -e GATC
yahs irugos.asm.hic.hap2.p_ctg.fa irugosHap2_hic.sorted.bam -o irugosHap2 -e GATC
./fa_to_dotplot.sh irugos/irugosHap1_scaffolds_final.fa irugosHap1yahs 2
./fa_to_dotplot.sh irugos/irugosHap2_scaffolds_final.fa irugosHap2yahs 2



juicer pre irugosHap1.bin irugosHap1_scaffolds_final.agp irugos.asm.hic.hap1.p_ctg.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > irugosHap1.alignments_sorted.txt
samtools faidx irugosHap1_scaffolds_final.fa
java -jar -Xmx32G ../juicer_tools.2.18.00.jar  pre irugosHap1.alignments_sorted.txt irugosHap1.hic <( cut -f1,2 irugosHap1_scaffolds_final.fa.fai )


cat irugos.asm.hic.hap1.p_ctg.fa irugos.asm.hic.hap2.p_ctg.fa > irugos.asm.hic.bothhaps.p_ctg.fa
bwa index irugos.asm.hic.bothhaps.p_ctg.fa
bwa mem -5SP -t 64 irugos.asm.hic.bothhaps.p_ctg.fa AV9240009_Ischaemum_rugosum_HiC_pairsR1.fq.gz AV9240009_Ischaemum_rugosum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o irugos_bothhapsHiC.bam

../HapHiC/utils/filter_bam irugos_bothhapsHiC.bam 1 --nm 3 --threads 14 | samtools view - -b -@ 14 -o irugos_bothhapsHiC.filtered.bam

conda activate haphic
# HapHic scaffolding
../HapHiC/haphic pipeline irugos.asm.hic.bothhaps.p_ctg.fa irugos_bothhapsHiC.filtered.bam 20 --RE "GATC" --correct_nrounds 2 --remove_allelic_links 2 --threads 32 --processes 16 --max_inflation 7.0 --flank 0 --bin_size 1000
## try with fewer chromosomes?
../HapHiC/haphic pipeline irugos.asm.hic.bothhaps.p_ctg.fa irugos_bothhapsHiC.filtered.bam 18 --RE "GATC" --correct_nrounds 2 --remove_allelic_links 2 --threads 32 --processes 16 --max_inflation 7.0 --flank 0 --bin_size 1000 --outdir chr18
## try providing the gfa of each haplotype?
../HapHiC/haphic pipeline irugos.asm.hic.bothhaps.p_ctg.fa irugos_bothhapsHiC.filtered.bam 18 --RE "GATC" --correct_nrounds 2 --remove_allelic_links 2 --threads 32 --processes 16 --max_inflation 7.0 --flank 0 --bin_size 1000 --gfa "etrips2hap.asm.hic.hap1.p_ctg.gfa,etrips2hap.asm.hic.hap2.p_ctg.gfa" --outdir irugosGfa18chr

## try using the unfiltered bam?
../HapHiC/haphic pipeline irugos.asm.hic.bothhaps.p_ctg.fa irugos_bothhapsHiC.bam 20 --threads 32 --gfa "irugos.asm.hic.hap1.p_ctg.gfa,irugos.asm.hic.hap2.p_ctg.gfa" --outdir irugosUnfiltered
## try inflation 10.0, had only made 11 chromsomes with the last one
../HapHiC/haphic pipeline irugos.asm.hic.bothhaps.p_ctg.fa irugos_bothhapsHiC.bam 20 --threads 32 --gfa "irugos.asm.hic.hap1.p_ctg.gfa,irugos.asm.hic.hap2.p_ctg.gfa" --outdir irugosUnfiltered --inflation 10.0


bash juicebox.sh
./fa_to_dotplot.sh /workdir/mcs368/panand_assemblies/hic/irugos/04.build/scaffolds.fa irugosHapHiC 2



# Generate the final FASTA file for the scaffolds
../HapHiC/utils/juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp asm.fa











