conda activate miniprot ## on cbsu panand_htt/aa_panand/

while read species six ploidy
do
echo $( ls ../genes/${species}*.gff3 )
gffread -y ${six}.aa.fa -g ../genomes/${species}.fasta ../genes/${species}*.gff3

done < ../panand_sp_ploidy.txt

## scaffolds didn't work 
gffread -y sbicol.aa.fa -g ../genomes/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fasta ../genes/old/Sbicolor_454_v3.0.1_BTx623PA.1_final.gff3 


## also do paspalum!!
gffread -y pvagin.aa.fa -g ../genomes/Pvaginatum_672_v3.0.fa ../genes/old/Pvaginatum_672_v3.1.gene.gff3 

