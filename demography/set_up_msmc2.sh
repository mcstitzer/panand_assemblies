


## xm01
conda activate pip


bcftools view -g ^miss -o filtered.vcf.gz -O z ../output/tdacn1_reads.vcf.gz
tabix -p vcf filtered.vcf.gz

bcftools view -m2 -M2 -v snps filtered.vcf.gz -Oz -o filtered_biallelic.vcf.gz
tabix filtered_biallelic.vcf.gz 


chr=chr18
python msmc-tools/generate_multihetsep.py --as_phased --chr ${chr} filtered_biallelic.vcf.gz > tdacn1.multihetsep.${chr}.txt 
./msmc2_Linux -t 2 -p 1*2+15*1+1*2 -o tdacn1.${chr}.msmc2 -I 0,1 <( awk '$3>0' tdacn1.multihetsep.${chr}.txt | head -100000 )


awk '$3>0' tdacn1.multihetsep.txt  | head -1000000 > tdacn1.multihetsep.1000k.txt

./msmc2_Linux -t 2 -p 1*2+15*1+1*2 -o tdacn1.chr1.1000k.msmc2 -I 0,1 tdacn1.multihetsep.1000k.txt 



## try for biggest etripsacoides chromosome

bcftools view -g ^miss -o filtered.vcf.gz -O z ../output/etrips_reads.vcf.gz
tabix -p vcf filtered.vcf.gz

bcftools view -m2 -M2 -v snps filtered.vcf.gz -Oz -o filtered_biallelic.vcf.gz
tabix filtered_biallelic.vcf.gz 


chr=scaf_38
python msmc-tools/generate_multihetsep.py --as_phased --chr ${chr} filtered_biallelic.vcf.gz > etrips.multihetsep.${chr}.txt 
./msmc2_Linux -t 2 -p 1*2+15*1+1*2 -o etrips.${chr}.msmc2 -I 0,1 <( awk '$3>0' etrips.multihetsep.${chr}.txt | head -500000 )


## try for biggest ppanic chromosome

bcftools view -g ^miss -o filtered.vcf.gz -O z ../output/ppanic_reads.vcf.gz
tabix -p vcf filtered.vcf.gz

bcftools view -m2 -M2 -v snps filtered.vcf.gz -Oz -o filtered_biallelic.vcf.gz
tabix filtered_biallelic.vcf.gz 


chr=ctg_1
python msmc-tools/generate_multihetsep.py --as_phased --chr ${chr} filtered_biallelic.vcf.gz > ppanic.multihetsep.${chr}.txt 
./msmc2_Linux -t 2 -p 1*2+15*1+1*2 -o ppanic.${chr}.msmc2 -I 0,1 <( awk '$3>0' ppanic.multihetsep.${chr}.txt | head -500000 )




## now, need a gvcf of h contortus!! start with just four homologs for the first chromsome
## will i be able to output a 4 sample msmc??? i think i can figure it out.
## maybe i need to merge gvcf so there are three samples???


