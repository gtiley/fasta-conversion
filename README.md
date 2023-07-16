# fasta-conversion
command-line tools for converting fasta files with population samples into other data formats

## Instructions
The `fasta2Strcuture.pl` script can convert a fasta file with multiple haplotypes per individual into an integer format accepted by STRUCTURE. The allelic dosage for a site is also calculated and this will generate several files suitable for PCA analyses. I use the resulting PCA matrices and pass those to adegenet for imputing the missing values and analysis.

An example command with the provided test data is 
```
perl fasta2Structure.pl -i test.fasta -o test -p test.popmap -c 0 -n 0 -mac 2 -minsites 20
```

In cases where you need to thin the structure file to 1 SNP per locus to alleviate tight linkage, you can use
```
perl thinStructure.pl test
```
All that is needed is the prefix name of the prefix.str file, since the locus names and site positions should be in the header of the structure file (see notes below). Note that this does not work if you did not provide the coordinate data to `fasta2Sructure.pl` or if your loci are whole chromosomes.


### Some notes on inputs:
 1. The fasta file as a strict naming convention. The fasta header should be ">SP_HAPLOTYPE", where SP is the species or individual name and HAPLOTYPE is an integer representing the haplotype sampled from that individual. The "\_" delimits the name from the haplotype number, so an "\_" in an individual name should be avoided.
 2. The popmap file give the individual name, the populaion, and the ploidy level (optional) as a tab-delimited text file. Setting `-n 0` cause the ploidy level to be read from the file. Otherwise, `-n` can be set to an integer if all inividuals have the same ploidy. The population names are currently not used.
 3. There is an option for including site coordinates with `-c 0`. Keeping 0 just means that an ID will be assigned for the STRUCTURE file based on its position in the alignment. Otherwise you want to pass it a file with information about the locus name and site position. That can be a `\*.pos` file generate by the `vcf2fasta.pl` script or a RAxML-formatted partition file if starting from a concatenated fasta file. There should be no `\_` characters in the locus names as they will be used to delimit site positions in loci.
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

### The \*.pos file
This is a custom format that I use for various scripts. It is a simple two-column list of loci and sites in order of the fasta file from left to right that could look like this:
|------|------|
|group1|294708|
|group1|297223|
|group1|297229|
|------|------|
It is generated by the `vcf2fasta.pl` script that creates FASTA files from the VCF format. It can be invoked like this:
```
./vcf2fasta.pl -i murinusClade.keepBad1.mac3.FS6.vcf -o Mmur_all.fasta -t 2 -d 10000 -b 1\n";
```
The input vcf (-i) should not be gzipped. The output file (-o) is always FASTA format. The -t and -d flags control how to thin SNPs and the -b flag can be used to include all sites or only biallelic sites. This will generate a \*.pos file. In general, I do not recommend the script as there are many highly-cited programs that can do the same operations for you, but I am working on some other things where this was convenient. This will not generate the input you want for the PCA with polyploids since that requires a concatenated fasta files of the haplotype sequences with the "\_" character in individual names to denote haplotypes. But, that is the origin of what a \*.pos file is.




