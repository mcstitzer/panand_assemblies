
conda activate panand_assemblies
## first time through, hard coded for udigit

mkdir -p udigit
cd udigit
## get data
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/PanAnd_082025/AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/PanAnd_082025/AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz .

bwa index Ud-Pasquet1171-DRAFT-PanAnd-1.0.fasta

## 5SP deals with weird hic distribution of reads
bwa mem -t 64 -5SP Ud-Pasquet1171-DRAFT-PanAnd-1.0.fasta AV9240014_Urelytrum_digitatum_HiC_pairsR1.fq.gz AV9240014_Urelytrum_digitatum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o udigit_hic.sorted.bam 
samtools index udigit_hic.sorted.bam

## get site positions of restriction enzyme on scaffolds
python ../juicer/misc/generate_site_positions.py DpnII GATC Ud-Pasquet1171-DRAFT-PanAnd-1.0.fasta 

source /programs/miniconda3/bin/activate hic_qc
/programs/hic_qc/hic_qc.py -b smicro_hic.namesorted.bam -o smicro


../pairix/util/bam2pairs/bam2pairs -c smicro.chromsizes.txt smicro_hic.namesorted.bam smicro.hicup
java -jar -Xmx32G ../juicer_tools.2.18.00.jar pre smicro.hicup.bsorted.pairs smicro.bam2pairs.hic smicro.chromsizes.txt

## try a easier genome come on michelle
## smicro
mkdir -p smicro
cd smicro
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz .
scp cbsublfs1:/data4/Incoming/Corteva/Aug2025_PanAndtransfers/082525_PanAnd/AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz .

bwa index Sm-PI203595-DRAFT-PanAnd-1.0.fasta

## 5SP deals with weird hic distribution of reads
bwa mem -t 64 -5SP Sm-PI203595-DRAFT-PanAnd-1.0.fasta AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz | samblaster | samtools sort -@ 16 -o smicro_hic.sorted.bam 
samtools index smicro_hic.sorted.bam

## get site positions of restriction enzyme on scaffolds
python ../juicer/misc/generate_site_positions.py DpnII GATC Sm-PI203595-DRAFT-PanAnd-1.0.fasta 

yahs Sm-PI203595-DRAFT-PanAnd-1.0.fasta smicro_hic.sorted.bam -o smicro -e GATC

# samblaster: Pair Type        Type_ID_Count   %Type/All_IDs Dup_ID_Count  %Dups/Type_ID_Count  %Dups/All_Dups  %Dups/All_IDs
# samblaster: ---------------------------------------------------------------------------------------------------------------
# samblaster: Both Unmapped          790626         0.446              0           0.000             0.000          0.000
# samblaster: Orphan/Singleton      2131632         1.203         237014          11.119             0.323          0.134
# samblaster: Both Mapped         174219198        98.350       73135776          41.979            99.677         41.287
# samblaster: Total               177141456       100.000       73372790          41.420           100.000         41.420
# samblaster:
# samblaster: Marked    73372790 of  177141456 (41.420%) total read ids as duplicates using 2021648k memory in 6M6S(365.820S) CPU seconds and 54M39S(3279S) wall time.


(juicer pre smicro.bin smicro_scaffolds_final.agp Sm-PI203595-DRAFT-PanAnd-1.0.fasta.fai | sort -k2,2d -k6,6d -T ./ --parallel=64 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

(java -jar -Xmx32G ../juicer_tools.2.18.00.jar pre alignments_sorted.txt out.hic.part <( cut -f1,2 Sm-PI203595-DRAFT-PanAnd-1.0.fasta.fai ) ) && (mv out.hic.part out.hic)


# The out_JBAT.hic can be loaded into Juicebox along with out_JBAT.assembly for manual editing.
# 
# NOTE 3: if your total assembly size is larger than 2Gb, there will be a scale factor applied to the HiC contact map file to make it loadable by Juicebox. This scale factor can be found in the log file (out_JBAT.log) and always a power of two, e. g. [I::main_pre] scale factor: 4. You need to set this parameter in Juicebox through Assembly > Set Scale. Otherwise the HiC contact map and assembly file will not match.
# 
samtools sort -n -@ 64 -o smicro_hic.namesorted.bam smicro_hic.sorted.bam


juicer pre smicro_hic.namesorted.bam smicro_scaffolds_final.agp Sm-PI203595-DRAFT-PanAnd-1.0.fasta.fai -o smicroHiC -q 0
juicer_tools pre -q 0 smicroHiC.txt Sm-PI203595-DRAFT-PanAnd-1.0.fasta.fai smicroHiC.hic

awk '$5=="W"{len[$1]+=$8-$7+1} END{for (c in len) print c, len[c]}' smicroHiC.assembly.agp

cut -f1,2 Sm-PI203595-DRAFT-PanAnd-1.0.fasta.fai > smicro.chromsizes.txt
# Append assembly
total_len=$(awk '$5=="W"{len[$1]+=$8-$7+1} END{print len["assembly"]}' smicroHiC.assembly.agp)
echo -e "assembly\t${total_len}" >> smicro.chromsizes.txt


juicer_tools pre -q 0 smicro_hic.namesorted.bam smicro.chromsizes.txt smicroHiC.hic

juicer pre -q 1 smicro_hic.namesorted.bam smicro.chromsizes.txt smicroHiC.hic


pairtools parse --assembly Sm-PI203595-DRAFT-PanAnd-1.0.fasta --chroms-path smicro.chromsizes.txt smicro_hic.namesorted.bam \
    | pairtools sort \
    | pairtools dedup \
    > smicro.pairs
    
    
## from way above
java -Xmx32g -jar ../juicer_tools.2.18.00.jar pre -q 0 alignments_sorted.txt smicroHic.hic smicro.chromsizes.txt





java -jar -Xmx32G ../juicer_tools.2.18.00.jar pre -q 0 alignments_sorted.txt out.hic.part smicro.chromsizes.txt


###### this took maybe 2-3 hours??
pairtools parse -c smicro.chromsizes.txt smicro_hic.namesorted.bam | \
pairtools sort | \
pairtools dedup > reads.dedup.pairs

## 
pairtools select '(pair_type == "UU") and not (is_duplicate)' reads.dedup.pairs | \
awk 'BEGIN{OFS="\t"} {str1=($6=="+" ? 0 : 1); str2=($7=="+" ? 0 : 1); print str1, $2, $3, 0, str2, $4, $5, 1}' \
> alignments_for_juicer.txt

java -jar -Xmx64G ../juicer_tools.2.18.00.jar pre alignments_for_juicer.txt out.hic smicro.chromsizes.txt




## switch to pretext, coordinate sort
samtools view -h smicro_hic.sorted.bam | PretextMap -o smicro.pretext --sortby length --sortorder descend --mapq 0










