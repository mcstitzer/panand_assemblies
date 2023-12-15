

## use TR repeatmasked gff now

while read filename genome ploidy
do
if [ -f ${filename}_TandemRepeat.gff3 ]
then
## -A removes entire TE that overlaps with a tandem repeat - I think this is okay, because tandem boundaries aren't perfect
bedtools subtract -A -a ../../repeats_final/${filename}_EDTA.gff3 -b ${filename}_TandemRepeat.gff3 > ${filename}_EDTA.filteredTandemRepeat.gff3
bedtools sort -i <(cat ${filename}_EDTA.filteredTandemRepeat.gff3 ${filename}_TandemRepeat.gff3 ) > ${filename}_EDTATandemRepeat.gff3
fi
done < ../../panand_sp_ploidy.txt

## b73 is named differently??

genome=zmB735
filename=Zm-B73-REFERENCE-NAM-5.0
## icky the B73_chr1 is on these
bedtools subtract -A -a <( zcat ../../../maize_TE_load_fitness/genomes_and_annotations/tes/new_edta_tes/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.split.gff3.gz | sed 's/^B73_//g' ) -b ${filename}_TandemRepeat.gff3 > ${filename}_EDTA.filteredTandemRepeat.gff3
bedtools sort -i <(cat ${filename}_EDTA.filteredTandemRepeat.gff3 ${filename}_TandemRepeat.gff3 ) > ${filename}_EDTATandemRepeat.gff3


## cserru doesn't have tandem repeats (it's the nanopore assembly)
## just copy the unchanged EDTA output anyways
cp ../../repeats_final/Cs-KelloggPI219580-DRAFT-PanAnd-1.0_EDTA.gff3 Cs-KelloggPI219580-DRAFT-PanAnd-1.0_EDTATandemRepeat.gff3 

## okay, do for smicro, snutan, telega, ttrian, tzolop, vcuspi
cp ../../repeats_final/Sm-PI203595-DRAFT-PanAnd-1.0_EDTA.gff3 Sm-PI203595-DRAFT-PanAnd-1.0_EDTATandemRepeat.gff3 
cp ../../repeats_final/Sn-CAM1369-DRAFT-PanAnd-1.0_EDTA.gff3 Sn-CAM1369-DRAFT-PanAnd-1.0_EDTATandemRepeat.gff3 
cp ../../repeats_final/Te-Pasquet1246-DRAFT-PanAnd-1.0_EDTA.gff3 Te-Pasquet1246-DRAFT-PanAnd-1.0_EDTATandemRepeat.gff3 
cp ../../repeats_final/Tt-AUB21_1-DRAFT-PanAnd-1.0_EDTA.gff3 Tt-AUB21_1-DRAFT-PanAnd-1.0_EDTATandemRepeat.gff3 
#cp ../repeats_final/Tz-DC_05_58_3A-DRAFT-PanAnd-1.0_EDTA.gff3 Tz-DC_05_58_3A-DRAFT-PanAnd-1.0_EDTATandemRepeat.gff3 
cp ../../repeats_final/Vc-Pasquet1098-DRAFT-PanAnd-1.0_EDTA.gff3 Vc-Pasquet1098-DRAFT-PanAnd-1.0_EDTATandemRepeat.gff3 
cp ../../repeats_final/Sbicolor_454_v3.0.1_EDTA.gff3 Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel_EDTATandemRepeat.gff3 
