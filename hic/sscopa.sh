## 0.45% het, Bionano, 6% ks btwn dups, 2Gb
conda activate panand_assemblies

## check kmers first!
mkdir tmp
/programs/kmc-3.2.4/kmc -k21 -t10 -m64 -ci1 -cs10000 sscopa_hifi.fastq.gz sscopakmc tmp/
/programs/kmc-3.2.4/kmc_tools transform sscopakmc histogram sscopa.histo -cx10000

../genomescope2.0/genomescope.R -i sscopa.histo -o sscopa_genomescope -k 21


