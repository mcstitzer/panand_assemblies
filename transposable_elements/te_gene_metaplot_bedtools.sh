

bedtools flank -i genes.bed -g genome.sizes -l 20000 -r 20000 > flanked_genes.bed
bedtools makewindows -b flanked_genes.bed -w 100 > windows.bed

bedtools coverage -a windows.bed -b tes.bed > coverage.txt
