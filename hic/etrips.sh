## this should be relatively hard - 0.7% het, Bionano, 5% ks btwn dups, 4.6Gb
conda activate panand_assemblies

## check kmers first!
mkdir tmp
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 etrips_hifi.fastq.gz etripskmc tmp/
/programs/kmc-3.2.4/kmc_tools transform etripskmc histogram etrips.histo -cx10000

../genomescope2.0/genomescope.R -i etrips.histo -o etrips_genomescope -k 21
