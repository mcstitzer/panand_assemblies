docker1 run --rm \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \ **Replace this string with exactly one of the following [WGS,WES,PACBIO,ONT_R104,HYBRID_PACBIO_ILLUMINA]**
  --ref=/input/Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta \
  --reads=/input/avirgi_reads-mapped.sorted.bam \
  --output_vcf=/output/avirgi_reads.vcf \
  --output_gvcf=/output/avirgi_reads.gvcf \
  --num_shards=200 \ **This will use all your cores to run make_examples. Feel free to change.**
  --logging_dir=/output/logs 
  
  
  
  
  
  
  \ **Optional. This saves the log output for each stage separately.
  --haploid_contigs="chrX,chrY" \ **Optional. Heterozygous variants in these contigs will be re-genotyped as the most likely of reference or homozygous alternates. For a sample with karyotype XY, it should be set to "chrX,chrY" for GRCh38 and "X,Y" for GRCh37. For a sample with karyotype XX, this should not be used.
  --par_regions_bed="/input/GRCh3X_par.bed" \ **Optional. If --haploid_contigs is set, then this can be used to provide PAR regions to be excluded from genotype adjustment. Download links to this files are available in this page.
  --dry_run=false **Default is false. If set to true, commands will be printed out but not executed.
  
  
  
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta \
  --reads=/input/avirgi_reads_withQ.bam \
  --output_vcf=/output/avirgi_reads.vcf.gz \
  --output_gvcf=/output/avirgi_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs
  
  
scp mcs368@cbsublfs1:/data4/Incoming/IowaStateU/rawdata_cinta_2024-06-07/avirgi/m64029_200818_154952.merged.hifi.bam .
scp mcs368@cbsublfs1:/data4/Incoming/IowaStateU/rawdata_cinta_2024-06-07/avirgi/m64029_200829_155609.merged.hifi.bam .



samtools fastq m64029_200818_154952.merged.hifi.bam | minimap2 -ax map-pb -t 160 Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta - | samtools sort -@ 40 -o m64029_200818_154952.merged.hifi.sorted.bam
samtools fastq m64029_200829_155609.merged.hifi.bam | minimap2 -ax map-pb -t 160 Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta - | samtools sort -@ 40 -o m64029_200829_155609.merged.hifi.sorted.bam

samtools merge -@ 200 avirgi_reads_withQ.bam m64029_200818_154952.merged.hifi.sorted.bam m64029_200829_155609.merged.hifi.sorted.bam

# Index the merged BAM file
samtools index avirgi_reads_withQ.bam



#####
shortname=avirgi
refgen=Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta
bam1=/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/AJ718_Andropogon_virginicus/Andropogon_virginicus.CELL01.hifi.bam
name1=$(basename "$bam1" .bam)
scp mcs368@cbsublfs1:${bam1} .
samtools fastq ${name1}.bam | minimap2 -ax map-pb -t 160 $refgen - | samtools sort -@ 40 -o ${name1}.sorted.bam

bam2=/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/AJ718_Andropogon_virginicus/Andropogon_virginicus.CELL02.hifi.bam
name2=$(basename "$bam2" .bam)
scp mcs368@cbsublfs1:${bam2} .
samtools fastq ${name2}.bam | minimap2 -ax map-pb -t 160 ${refgen} - | samtools sort -@ 40 -o ${name2}.sorted.bam

samtools merge -@ 200 ${shortname}_reads_withQ.bam ${name1}.sorted.bam ${name2}.sorted.bam

# Index the merged BAM file
samtools index ${shortname}_reads_withQ.bam



bcftools query -f '%POS\t%INFO/END\n' <(zcat avirgi_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
#897926186
bcftools view -v snps -g het <(zcat avirgi_reads.g.vcf.gz) | grep -v "^#" | wc -l
#2008
bcftools view -v snps -g hom <(zcat avirgi_reads.g.vcf.gz) | grep -v "^#" | wc -l
#8020



####
shortname=ppanic
refgen=Pi-Clark-DRAFT-PanAnd-1.0.fasta
bam1=/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Andro1-Pogonatherum_paniceum/ccs/m54334U_210629_201606.ccs.bam
name1=$(basename "$bam1" .bam)
scp mcs368@cbsublfs1:${bam1} .
samtools fastq ${name1}.bam | minimap2 -ax map-pb -t 160 $refgen - | samtools sort -@ 40 -o ${name1}.sorted.bam

mv ${name1}.sorted.bam ${shortname}_reads_withQ.bam
#samtools merge -@ 200 ${shortname}_reads_withQ.bam ${name1}.sorted.bam ${name2}.sorted.bam

# Index the merged BAM file
samtools index ${shortname}_reads_withQ.bam


  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ppanic_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ppanic_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ppanic_reads.g.vcf.gz) | grep -v "^#" | wc -l
#686561928
#4291952
#4659390






####
shortname=rrottb
refgen=Rr-Malcomber3106-DRAFT-PanAnd-1.0.fasta
bam1=/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Andro2-Rhytacne_rottboelloides/ccs/m54334U_210715_024754.hifi_reads.bam
name1=$(basename "$bam1" .bam)
scp mcs368@cbsublfs1:${bam1} .
samtools fastq ${name1}.bam | minimap2 -ax map-pb -t 160 $refgen - | samtools sort -@ 40 -o ${name1}.sorted.bam

mv ${name1}.sorted.bam ${shortname}_reads_withQ.bam
#samtools merge -@ 200 ${shortname}_reads_withQ.bam ${name1}.sorted.bam ${name2}.sorted.bam

# Index the merged BAM file
samtools index ${shortname}_reads_withQ.bam


  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
#1392445289
#379315
#536126



####
shortname="udigit"
refgen="Ud-Pasquet1171-DRAFT-PanAnd-1.0.fasta"
bam_files=(
    "/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Urelytrum/Urelytrum.CELL01.hifi.bam"
    "/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Urelytrum/Urelytrum.CELL02.hifi.bam"
    "/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Urelytrum/Urelytrum.CELL03.hifi.bam"
    "/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Urelytrum/Urelytrum.CELL04.hifi.bam"
    "/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Urelytrum/Urelytrum.CELL05.hifi.bam"
    "/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Urelytrum/Urelytrum.CELL06.hifi.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
#4408421481
#379315
#536126


####
shortname="smicro"
refgen="Sm-PI203595-DRAFT-PanAnd-1.0.fasta"
bam_files=(
    "/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-3/Grass-3__Schizachyrium_condensatum-microstachyum-PI_203595_115pM_01/CCS/m64041_210111_223935.hifi_reads.bam"
    "/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-3/Grass-3__Schizachyrium_condensatum-microstachyum-PI_203595_120pM_03/CCS/m64041_210119_014433.hifi_reads.bam"
    "/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-3/Grass-3__Schizachyrium_condensatum-microstachyum-PI_203595_130pM_04/CCS/m64041_210208_025008.hifi_reads.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l



##################
shortname="achine"
refgen="Ac-Pasquet1232-DRAFT-PanAnd-1.0.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-1/Andropogon_chinensis_Pasquet_1232_01/CCS/m64041_201230_230219.hifi_reads.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-1/Andropogon_chinensis_Pasquet_1232_02/CCS/m64041_210101_052407.hifi_reads.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-1/Andropogon_chinensis-Pasquet_1232_120pM_03/CCS/m64041_210105_070557.hifi_reads.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-1/Andropogon_chinensis-Pasquet_1232_140pM_04/CCS/m64041_210106_131439.hifi_reads.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l






##################
shortname="tdacs1"
refgen="Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-SouthernTrip/SouthernTrip.CELL01.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-SouthernTrip/SouthernTrip.CELL02.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-SouthernTrip/SouthernTrip.CELL03.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-SouthernTrip/SouthernTrip.CELL04.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-SouthernTrip/SouthernTrip.CELL05.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-SouthernTrip/SouthernTrip.CELL06.hifi.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l

## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta.fai  -w 100000) -b - -c > ~/transfer/tdacs1_snp_summary_100kb.bed
## fine do properly and get number genotyped sites
zcat ${shortname}_reads.g.vcf.gz | \
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' 2>/dev/null | \
awk '{OFS="\t"; if ($2 <= $3) print $1, $2-1, $3}' | \
bedtools intersect -a <(bedtools makewindows -g ../input/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a.fasta.fai -w 100000) -b - -c > ~/transfer/tdacs1_genotyped_sites_summary_100kb.bed
# Step 3: Combine the results and calculate normalized SNP density
paste ~/transfer/tdacs1_snp_summary_100kb.bed ~/transfer/tdacs1_genotyped_sites_summary_100kb.bed | \
awk '{if ($8 != 0) print $1, $2, $3, $4, $8, $4/$8; else print $1, $2, $3, $4, $8, "NA"}' > ~/transfer/tdacs1_normalized_snp_density_100kb.bed


##################
shortname="tdacn1"
refgen="Td-KS_B6_1-REFERENCE-PanAnd-2.0a.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL01.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL02.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL03.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL04.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL05.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL06.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL07.hifi.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/Td-KS_B6_1-REFERENCE-PanAnd-2.0a.fasta.fai  -w 100000) -b - -c > ~/transfer/tdacs1_snp_summary_100kb.bed


##################
shortname="tdacn1"
refgen="Td-KS_B6_1-REFERENCE-PanAnd-2.0a.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL01.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL02.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL03.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL04.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL05.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL06.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Tripsacum_dactyloides-NorthernTrip/NorthernTrip.CELL07.hifi.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed







##################
shortname="rtuber"
refgen="Rt-Layton_Zhong169-DRAFT-PanAnd-1.0.fasta"
bam_files=(
"/data4/PanAnd4/RawSeqData/Pacbio/andropogoneae/Andro5-Coelorachis_tuberculosa_NewLib/ccs/m54334Ue_211106_171621.hifi_reads.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed






######

shortname="etrips"
refgen="Et-Layton_Zhong168-DRAFT-PanAnd-1.0.fasta"
bam_files=(
"/data4/PanAnd4/RawSeqData/Pacbio/andropogoneae/Arizona/Elionorus_tripsacoides_HiFi_100pM_01/CCS/m64041_211030_034647.hifi_reads.bam"
"/data4/PanAnd4/RawSeqData/Pacbio/andropogoneae/Arizona/Elionorus_tripsacoides_HiFi_100pM_02/CCS/m64313e_211031_232409.hifi_reads.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed





######

shortname="irugos"
refgen="Ir-Pasquet1136-DRAFT-PanAnd-1.0.fasta"
bam_files=(
"/data4/PanAnd4/RawSeqData/Pacbio/andropogoneae/Andro6-Ischaemum_rugosum/ccs/m54334U_210907_163710.hifi_reads.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed



######

shortname="crefra"
refgen="Cr-AUB069-DRAFT-PanAnd-1.0.fasta"
bam_files=(
""
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed




######

shortname="zmhuet"
refgen="Zh-RIMHU001-REFERENCE-PanAnd-1.0.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-6/Grass-6__Zea_mays-huehuetenagensis_130pM_01/CCS/m64041_210416_171959.hifi_reads.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-6/Grass-6__Zea_mays-huehuetenagensis_140pM_02/CCS/m64041_210419_092651.hifi_reads.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-6/Grass-6__Zea_mays-huehuetenagensis_140pM_03/CCS/m64041_210420_153700.hifi_reads.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed




######

shortname="zluxur"
refgen="Zl-RIL003-REFERENCE-PanAnd-1.0.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-6/Grass-6__Zea_mays-huehuetenagensis_130pM_01/CCS/m64041_210416_171959.hifi_reads.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-6/Grass-6__Zea_mays-huehuetenagensis_140pM_02/CCS/m64041_210419_092651.hifi_reads.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Arizona/Grass-6/Grass-6__Zea_mays-huehuetenagensis_140pM_03/CCS/m64041_210420_153700.hifi_reads.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed






######

shortname="zTIL25"
refgen="Zx-TIL25-REFERENCE-PanAnd-1.0.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_mexicana-TIL25/TIL25.HiFi.CELL01.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_mexicana-TIL25/TIL25.HiFi.CELL02.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_mexicana-TIL25/TIL25.HiFi.CELL03.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_mexicana-TIL25/TIL25.HiFi.CELL04.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed


# samtools faidx Zx-TIL25-REFERENCE-PanAnd-1.0.fasta -r chrlist.txt -o Zx-TIL25-REFERENCE-PanAnd-1.0.chrOnly.fasta
# samtools faidx Zx-TIL25-REFERENCE-PanAnd-1.0.chrOnly.fasta 

shortname="zTIL25chr"
refgen="Zx-TIL25-REFERENCE-PanAnd-1.0.chrOnly.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_mexicana-TIL25/TIL25.HiFi.CELL01.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_mexicana-TIL25/TIL25.HiFi.CELL02.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_mexicana-TIL25/TIL25.HiFi.CELL03.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_mexicana-TIL25/TIL25.HiFi.CELL04.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
#    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed




####
 
shortname="zTIL18"
refgen="Zx-TIL18-REFERENCE-PanAnd-1.0.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL01.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL02.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL03.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL04.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL05.hifi.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed

####
 
 # samtools faidx Zx-TIL18-REFERENCE-PanAnd-1.0.fasta -r chrlist.txt -o Zx-TIL18-REFERENCE-PanAnd-1.0.chrOnly.fasta
# samtools faidx Zx-TIL18-REFERENCE-PanAnd-1.0.chrOnly.fasta 

shortname="zTIL18chr"
refgen="Zx-TIL18-REFERENCE-PanAnd-1.0.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL01.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL02.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL03.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL04.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_mexicana-TIL18/TIL18.CELL05.hifi.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
 #   scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed


####
 
shortname="telega"
refgen="Te-Pasquet1246-DRAFT-PanAnd-1.0.fasta"
bam_files=(
"/data4/PanAnd4/RawSeqData/Pacbio/andropogoneae/Andro3-Thelopogon_elegans/ccs/m54334U_210727_193135.hifi_reads.bam"
"/data4/PanAnd4/RawSeqData/Pacbio/andropogoneae/Andro3-Thelopogon_elegans_NewLib/ccs/m64310e_211025_162532.hifi_reads.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed




######
####
  # samtools faidx Zv-TIL01-REFERENCE-PanAnd-1.0.fasta -r chrlist.txt -o Zv-TIL01-REFERENCE-PanAnd-1.0.chrOnly.fasta
# samtools faidx Zv-TIL01-REFERENCE-PanAnd-1.0.chrOnly.fasta 

shortname="TIL01chr"
refgen="Zv-TIL01-REFERENCE-PanAnd-1.0.chrOnly.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_parviglumis-TIL01/TIL01.CELL01.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_parviglumis-TIL01/TIL01.CELL02.hifi.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zmay_parviglumis-TIL01/TIL01.CELL03.hifi.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed




######
####
  # samtools faidx Zd-Gigi-REFERENCE-PanAnd-1.0.fasta -r chrlist.txt -o Zd-Gigi-REFERENCE-PanAnd-1.0.chrOnly.fasta
# samtools faidx Zd-Gigi-REFERENCE-PanAnd-1.0.chrOnly.fasta 

shortname="zdgigichr"
refgen="Zd-Gigi-REFERENCE-PanAnd-1.0.chrOnly.fasta"
bam_files=(
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_diploperennis-GiGi/DIP.CELL01.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_diploperennis-GiGi/DIP.CELL02.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_diploperennis-GiGi/DIP.CELL03.bam"
"/data3/PanAnd3/RawSeqData/Pacbio/andropogoneae/Zea_diploperennis-GiGi/DIP.CELL04.bam"
)

sorted_bams=()

# Loop through each BAM file
for bam in "${bam_files[@]}"; do
    # Get the base name without extension
    name=$(basename "$bam" .bam)
    
    # Transfer BAM file via scp
    scp "mcs368@cbsublfs1:${bam}" .

    # Convert BAM to FASTQ, align with minimap2, and sort
    samtools fastq "${name}.bam" | minimap2 -ax map-pb -t 160 "$refgen" - | samtools sort -@ 40 -o "${name}.sorted.bam"
    
    # Collect the sorted BAM file names for merging later
    sorted_bams+=("${name}.sorted.bam")
done

# Merge all sorted BAM files
samtools merge -@ 200 "${shortname}_reads_withQ.bam" "${sorted_bams[@]}"

# Index the merged BAM file
samtools index "${shortname}_reads_withQ.bam"

# Optional: Clean up intermediate sorted BAM files
# rm "${sorted_bams[@]}"
  docker1 run \
  -v /workdir/mcs368/panand_snps/input:/input \
  -v /workdir/mcs368/panand_snps/output:/output \
  google/deepvariant:latest \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=PACBIO \
  --ref=/input/${refgen} \
  --reads=/input/${shortname}_reads_withQ.bam \
  --output_vcf=/output/${shortname}_reads.vcf.gz \
  --output_gvcf=/output/${shortname}_reads.g.vcf.gz \
  --num_shards=200 \
  --logging_dir=/output/logs

bcftools query -f '%POS\t%INFO/END\n' <(zcat ${shortname}_reads.g.vcf.gz) 2>/dev/null | grep -P '\d+' | awk '{if ($2 > $1) sum += ($2 - $1 + 1)} END {print sum}'
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l
bcftools view -v snps -g hom <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | wc -l


## plot where on chr!!
bcftools view -v snps -g het <(zcat ${shortname}_reads.g.vcf.gz) | grep -v "^#" | awk '{OFS="\t"; print $1, $2-1, $2}'  | bedtools intersect -a <(bedtools makewindows -g ../input/${refgen}.fai  -w 100000) -b - -c > ~/transfer/${shortname}_snp_summary_100kb.bed

