
## cbsu
export PATH=/programs/RepeatMasker_4-1-0:$PATH

mkdir -p repeatmask_tandems/
cd repeatmask_tandems
while read -r altfa name ploidy
do
if [ -f ../${name}_filteredTandemRepeats.TRASH.fa ] ## some genomes are so bad they don't have tandem repeats
then
mkdir -p ${name}_trrm
if [ ! -f ${name}_trrm/${altfa}.fasta.tbl ]
then
RepeatMasker -e ncbi -pa 30 -q -no_is -norna -nolow -div 40 -dir ${name}_trrm -gff -lib ../${name}_filteredTandemRepeats.TRASH.fa -cutoff 225 ../../genomes/${altfa}.fasta 
fi
fi
done < ../../panand_sp_ploidy.txt
