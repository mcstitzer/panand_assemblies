
#for species in $(cut -f1 panand_sp_ploidy.txt)
while read species six ploidy
do
  genome=$(cat ../genomes/$species.fasta.fai | awk '{sum+=$2} END{print sum}')
 # AHHHHHHHHHHHH i just found out htat there are assemblies with over 25% Ns in the sequence. This is WILD. Let's count actual sequenced base pairs now
  enns=$(/programs/kentUtils/bin/faCount ../genomes/$species.fasta | tail -1 | cut -f 7 )
#  gNA=$(grep -v alt genomes/$species.fasta.gz.fai | awk '{sum+=$2} END{print sum}')
#  tebp=$(bedtools genomecov -i <( bedtools merge -i repeats_final/${species}_EDTA.gff3 ) -g genomes/${species}.fasta.gz.fai | grep -v "genome" | awk '{ if ( $2=="1" ) { sum+=$3;}} END{print sum}')
#  tNA=$(bedtools genomecov -i <( bedtools merge -i repeats_final/${species}_EDTA.gff3 ) -g genomes/${species}.fasta.gz.fai | grep -v "genome" | grep -v "alt" | awk '{ if ( $2=="1" ) { sum+=$3;}} END{print sum}')
#  echo $species $genome $gNA $tebp $tNA
#  printf '%s\t%s\t%s\t%s\t%s\n' "$species" "$genome" "$gNA" "$tebp" "$tNA"
  printf '%s\t%s\t%s\t%s\n' "$species" "$six" "$genome" "$enns" #"$tebp"

done < "../panand_sp_ploidy.txt" > panand_assembly_sizes.txt

## while loop is doing weird things and printing genome instead of gNA variable?????
#### lolol i'm stupid - there is no alt scaf in Ac????
## there go three hours of my life
