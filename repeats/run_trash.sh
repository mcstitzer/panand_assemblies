

#./TRASH_run.sh --def ../../genomes/Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta --par 11 --o /workdir/mcs368/panand_htt/trash/avirgi/

## get in TRASH directory of the conda - weird path not changing if up one?

while read -r altfa name ploidy
do
if [ ! -f /workdir/mcs368/panand_htt/trash/${name}/TRASH_${altfa}.fasta.gff ]
then
mkdir -p /workdir/mcs368/panand_htt/trash/${name}
## run anchors
./TRASH_run.sh --def ../../genomes/${altfa}.fasta --par 5 --o /workdir/mcs368/panand_htt/trash/${name}/ &
#./TRASH_run.sh --def ../../genomes/Td-KS_B6_1-REFERENCE-PanAnd-1.0.fasta --par 11 --o /workdir/mcs368/panand_htt/trash/tdacts/
#anchorwave proali -t 10 -i $refgff -as ${ref}.CDS.fa -r $reffa -a ${name}-${ref}.CDS.sam -ar ${ref}.CDS.sam -s genomes/${altfa}.fasta -n ${name}-${ref}-${ploidy} -R $ploidy -Q 1 -ns & ##-o ${name}-${ref}.maf -f ${name}-${ref}.f.maf &
fi
done < ../../panand_sp_ploidy.txt


#TRASH_Av-Kellogg1287_8-REFERENCE-PanAnd-1.0.fasta.gff
#./TRASH_run.sh --def ../../genomes/Td-KS_B6_1-REFERENCE-PanAnd-1.0.fasta --par 11 --o /workdir/mcs368/panand_htt/trash/tdacts/
