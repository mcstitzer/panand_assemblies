

conda activate haphic # or: source /path/to/conda/bin/activate haphic

##smicro.asm.hic.hap1.p_ctg.fa smicroHap1_hic.sorted.bam
## filter a little more
samtools sort -n -@ 64 -o smicroHap1_hic.namesorted.bam smicroHap1_hic.sorted.bam

samtools view smicroHap1_hic.namesorted.bam  -@ 14 -S -h -b -F 3340 -o smicroHap1_hic.haphic.bam

../HapHiC/utils/filter_bam smicroHap1_hic.haphic.bam 1 --nm 3 --threads 14 | samtools view - -b -@ 14 -o smicroHap1_hic.haphic.filtered.bam

../HapHiC/haphic pipeline smicro.asm.hic.hap1.p_ctg.fa smicroHap1_hic.haphic.filtered.bam 10

bash juicebox.sh

# Generate the final FASTA file for the scaffolds
../HapHiC/utils/juicer post -o out_JBAT out_JBAT.review.assembly out_JBAT.liftover.agp asm.fa

# For HapHiC scaffolding result
../HapHiC/haphic plot scaffolds.raw.agp HiC.filtered.bam
# For the AGP file generated after manual curation in Juicebox
../HapHiC/haphic plot out_JBAT.FINAL.agp HiC.filtered.bam



### reorder based on paspalum chromososmes?
# The preset can be `asm5` if the reference genome is well-assembled from the same species
$ minimap2 -x asm20 ref.fa asm.fa --secondary=no -t 28 -o asm_to_ref.paf
# `haphic refsort` can also be compatible with other aligners, like wfmash
$ wfmash ref.fa asm.fa -m -n 1 -S 1 -t 28 | cut -f 1-6,8- > asm_to_ref.paf

# By default, scaffolds are output based on the alphabetical order of the chromosome IDs of the reference genome
$ haphic refsort 04.build/scaffolds.raw.agp asm_to_ref.paf > scaffolds.refsort.agp
# You can specify the order by listing chromosome IDs of the reference genome separated by commas (no spaces)
$ haphic refsort 04.build/scaffolds.raw.agp asm_to_ref.paf --ref_order "chr1,chr2,chr3,chr4,..." > scaffolds.refsort.agp
# If you want to generate a new FASTA file (default name: `scaffolds.refsort.fa`) as well
$ haphic refsort 04.build/scaffolds.raw.agp asm_to_ref.paf --fasta asm.fa > scaffolds.refsort.agp
# If you want to run `haphic refsort` on manually curated `out_JBAT.FINAL.agp`
$ haphic refsort out_JBAT.FINAL.agp asm_to_ref.paf > scaffolds.refsort.agp



##################################################
##ACTUALLY - haphic needs both alleles in one!!!
cat smicro.asm.hic.hap1.p_ctg.fa smicro.asm.hic.hap2.p_ctg.fa >smicro.asm.hic.bothhaps.p_ctg.fa

bwa index smicro.asm.hic.bothhaps.p_ctg.fa
bwa mem -5SP -t 64 smicro.asm.hic.bothhaps.p_ctg.fa AV9240013_Schizachyrium_microstachyum_HiC_pairsR1.fq.gz AV9240013_Schizachyrium_microstachyum_HiC_pairsR2.fq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o smicroBothHaps.bam
## duh it can't index if not sorted
##samtools index smicroBothHaps.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam smicroBothHaps.bam 1 --nm 3 --threads 14 | samtools view - -b -@ 14 -o smicroBothHaps.filtered.bam
## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline smicro.asm.hic.bothhaps.p_ctg.fa smicroBothHaps.filtered.bam 20 --gfa "smicro.asm.hic.hap1.p_ctg.gfa,smicro.asm.hic.hap2.p_ctg.gfa" --outdir smicroBothHaps







#######

cat etrips2hap.asm.hic.hap1.p_ctg.fa etrips2hap.asm.hic.hap2.p_ctg.fa > etrips2hap.asm.hic.bothhaps.p_ctg.fa

bwa index etrips2hap.asm.hic.bothhaps.p_ctg.fa
bwa mem -5SP -t 64 etrips2hap.asm.hic.bothhaps.p_ctg.fa etrips.hic.read1.fq.gz etrips.hic.read2.fq.gz | samblaster | samtools view - -@ 14 -S -h -b -F 3340 -o etripsBothHaps.bam
## duh it can't index if not sorted
##samtools index smicroBothHaps.bam

# (2) Filter the alignments with MAPQ 1 (mapping quality ≥ 1) and NM 3 (edit distance < 3)
../HapHiC/utils/filter_bam etripsBothHaps.bam 1 --nm 3 --threads 14 | samtools view - -b -@ 14 -o smicroBothHaps.filtered.bam
## and can use the gfa to deal with depths and mis-pairing
../HapHiC/haphic pipeline etrips2hap.asm.hic.bothhaps.p_ctg.fa etripsBothHaps.filtered.bam 20 --gfa "etrips2hap.asm.hic.hap1.p_ctg.gfa,etrips2hap.asm.hic.hap2.p_ctg.gfa" --outdir etripsBothHaps













