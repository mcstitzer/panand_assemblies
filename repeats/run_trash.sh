

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


#wget https://github.com/vlothec/TRASH/raw/main/TRASH.v1.2.tar.gz
#tar -xzvf TRASH.v1.2.tar.gz

## lol i have to make the directory first!!! e.g. zTIL18 

./TRASH_run.sh --def ../Zx-TIL25-REFERENCE-PanAnd-1.0.fasta --par 48 --o /workdir/mcs368/panand_assemblies/repeats/zTIL25/


./TRASH_run.sh --def /workdir/mcs368/panand_assemblies/genomes/Zx-TIL18-REFERENCE-PanAnd-1.0.fasta --par 48 --o /workdir/mcs368/panand_assemblies/repeats/zTIL18/

./TRASH_run.sh --def /workdir/mcs368/panand_assemblies/genomes/Zn-PI615697-REFERENCE-PanAnd-1.0.fasta --par 48 --o /workdir/mcs368/panand_assemblies/repeats/znicar/


### rerun anything that is absent 
./TRASH_run.sh --def /workdir/mcs368/panand_assemblies/genomes/Vc-Pasquet1098-DRAFT-PanAnd-1.0.fasta --par 48 --o /workdir/mcs368/panand_assemblies/repeats/vcuspi/

./TRASH_run.sh --def /workdir/mcs368/panand_assemblies/genomes/Sn-CAM1369-DRAFT-PanAnd-1.0.fasta --par 48 --o /workdir/mcs368/panand_assemblies/repeats/snutan/
./TRASH_run.sh --def /workdir/mcs368/panand_assemblies/genomes/Cs-KelloggPI219580-DRAFT-PanAnd-1.0.fasta --par 48 --o /workdir/mcs368/panand_assemblies/repeats/cserru/
./TRASH_run.sh --def /workdir/mcs368/panand_assemblies/genomes/Tt-AUB21_1-DRAFT-PanAnd-1.0.fasta --par 48 --o /workdir/mcs368/panand_assemblies/repeats/ttrian/
./TRASH_run.sh --def /workdir/mcs368/panand_assemblies/genomes/Sm-PI203595-DRAFT-PanAnd-1.0.fasta --par 48 --o /workdir/mcs368/panand_assemblies/repeats/smicro/

