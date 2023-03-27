# fasta-conversion
command-line tools for converting fasta files with population samples into other data formats

## Instructions
The `fasta2Strcuture.pl` script can convert a fasta file with multiple haplotypes per individual into an integer format accepted by STRUCTURE. The allelic dosage for a site is also calculated and this will generate several files suitable for PCA analyses. I use the resulting PCA matrices and pass those to adegenet for imputing the missing values and analysis.

An example command with the provided test data is 
```
perl fasta2Structure.pl -i test.fasta -o test -p test.popmap -c 0 -n 0 -mac 2 -minsites 20
```

### Some notes on inputs:
 1. The fasta file as a strict naming convention. The fasta header should be ">SP_HAPLOTYPE", where SP is the species or individual name and HAPLOTYPE is an integer representing the haplotype sampled from that individual. The "\_" delimits the name from the haplotype number, so an "\_" in an individual name should be avoided.
 2. The popmap file give the individual name, the populaion, and the ploidy level (optional) as a tab-delimited text file. Setting `-n 0` cause the ploidy level to be read from the file. Otherwise, `-n` can be set to an integer if all inividuals have the same ploidy. The population names are currently not used.
 3. A currently unused option is `-c 0`. Keeping 0 just means that an ID will be assigned for the STRUCTURE file based on its position in the alignment.
 4. The *minimum allele count* filter is `-mac`. Sometimes data are filtered for mac and sometimes not depending on the application. This is a simple way to incorporate a mac filter for analyses of sructure, assuming the variants going in are filtered for other quality metrics.
 5. A filter the *minimum number of sites* is `-minsites`. Same as mac, the initial data might tolerate more missing values for some other analyses. I do not think it is too bad to allow a few missing sites for PCAs, as long as the imputed values are based on a decent amount of data. This is based on the total number of sequences (individuals * haplotypes) rather than the number of individuals.

### Some notes on outputs:
 1. Running the above example should generate the following files
  - test.str (The file for STRUCTURE, with population labels excluded and arbitrary locus labels included)
  - test.pcacomp2.txt (PCA matrix with all homozygous ref as 0 and at least one alt as 1)
  - test.pcacomp3.txt (PCA matrix treated with homozygous ref as 0 and homozygous alt as 2. Anything in between is 1)
  - test.pcacomp5.txt (PCA matrix PCA matrix treated with homozygous ref as 0 and homozygous alt as 4. A 50:50 ratio of 0:1 for a site is 2. Anything less than 50:50 but greater than 100:0 is 1 and anyhing greater than 50:50 but less than 0:100 is 3)
  - test.pcaraw.txt (PCA matrix without any type of compression. If the highest ploidy level for an individual is a 6, you will see a 6 for a homozygous alt. This will not be comparable with mixed ploidies though, since a homozygous alt in a tetraploid will get a 4).
 2. In the test, notice that a site was dropped due to the -minsites filter. For polyploid data, the mac may need to be increased, since there may be some non-independence among the genotype calls. Maybe a homozygous alt in one tetraploid individual sneaks through erroneously. A mac of 5 would kick this site out. 

### Some opinions:
 1. I question how reliable genotype estimates and any single point estimate of allelic dosage at a site is. My opinion is that there is a lot of noise beyond tetraploids. In my experience, projecting higher ploidy levels onto lower dimensions is helpful for analyses of population structure.
 2. This does not mean getting the allele frequencies accurate is not important for other analyses such as demographic reconstruction. This area of research for polyploids is still emerging.
