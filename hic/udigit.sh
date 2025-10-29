## very low het, Bionano, a little nanopore 4.2% ks btwn dups, 5.4Gb
conda activate panand_assemblies

## check kmers first!
mkdir tmp
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 udigit_hifi.fastq.gz udigitkmc tmp/
/programs/kmc-3.2.4/kmc_tools transform udigitkmc histogram udigit.histo -cx10000

../genomescope2.0/genomescope.R -i udigit.histo -o udigit_genomescope -k 21


