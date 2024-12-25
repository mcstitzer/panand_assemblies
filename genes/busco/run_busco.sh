for file in ./helixer_protein_aa/*.fasta
do
    busco -i $file -l poales_odb10 -o ${file%.fasta}_busco -m protein -c 200
done
