#!/usr/bin/perl -w
use strict;
#----------------------------------------------------------------------------------------#
#George P. Tiley
#26 April 2022
#updated 19 July 2023 to be less restrictive on partition file format
#contact: g.tiley@kew.org
#Convert fasta file to Structure format
#----------------------------------------------------------------------------------------#


##########################################################################################
#--------------------------------Input Flag Handling-------------------------------------#
##########################################################################################
my @checkArgs = ("i","o","p","c","n","mac","minsites");
my %passedArgs = ();
if (scalar(@ARGV) == 0)
{
	die "/*--------INPUT PARAMETERS--------*/\n
	-i STRING <input fasta file>
	-o STRING <output str file>
	-p STRING <population map file>
	-c 0 | STRING <coordinate file>
	-n 2 | 0 <ploidy level>
	-mac INT <minor allele count>
	-minsites INT <minimum number of indidivuals with data>
	\n/*--------EXAMPLE COMMANDS--------*/\n
	Case 1: Genomics coordinates are known for all sites and all sequences are diploid with heterozygotes represented by ambiguity codes.\n
	./fasta2Structure.pl -i in.fasta -o out -p in.popmap -c in.pos -n 2 -mac 2 -minsites 10\n
	
	Case 2: Coordinates are unknown but individuals have variable ploidy. Ploidy levels are represented by the third column in the population map file. In this case, all lines in the fasta file should be haploid (phased) and the individual names should repeat in the fasta headers (does not have to be exact matches).\n
	./fasta2Structure.pl -i in.fasta -o out -p in.popmap -c 0 -n 0 -mac 2 -minsites 10\n
	
	Note: Fasta files are assumed to have the whole sequence per line - e.g. not 60 based per line with breaks\n
	
	Note: A RAxML-sytle partition file can be used in place of the coordinate file
	";
}	
elsif (scalar(@ARGV) > 0)
{
	for my $i (0..(scalar(@ARGV) - 1))
	{
		if ($ARGV[$i] eq "-i")
		{
			$passedArgs{i} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-o")
		{
			$passedArgs{o} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-p")
		{
			$passedArgs{p} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-c")
		{
			$passedArgs{c} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-n")
		{
			$passedArgs{n} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-mac")
		{
			$passedArgs{mac} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-minsites")
		{
			$passedArgs{minsites} = $ARGV[$i+1];
		}
	}
	foreach my $arg (@checkArgs)
	{
		if (! exists $passedArgs{$arg})
		{
			die "/*--------MISSING PARAMETER--------*/\nMissing command line argument: $arg\n\n";
		}
	}
}

my $inFile = $passedArgs{i};
my $outFilePrefix = $passedArgs{o};
my $popMap = $passedArgs{p};
my $coordfile = $passedArgs{c};
my $ploidyflag = $passedArgs{n};

open OUT1,'>',"$outFilePrefix.str";
open OUT2,'>',"$outFilePrefix.pcaraw.txt";
open OUT3,'>',"$outFilePrefix.pcacomp5.txt";
open OUT4,'>',"$outFilePrefix.pcacomp3.txt";
open OUT5,'>',"$outFilePrefix.pcacomp2.txt";

##########################################################################################
##########################################################################################

##########################################################################################
#---------------------------------Define Populations-------------------------------------#
##########################################################################################
my %populations = ();
my %individuals = ();
my $npops = 0;
my %popList = ();
my %ploidylevels = ();
open FH1, '<',"$popMap";
while (<FH1>)
{
	my $line = $_;
	chomp $line;
	if ($line =~ m/^\S+.+/)
	{
		my @temp = ();
		@temp = split (/\t/,$line);
		if ($ploidyflag == 2)
		{
			my $individual = $temp[0];
			my $population = $temp[1];
			$individuals{$individual} = $population;
			push @{$populations{$population}}, $individual;
			if (! exists $popList{$population})
			{
				$popList{$population} = 1;
				$npops++;
			}
		}
		elsif ($ploidyflag == 0)
		{
			my $individual = $temp[0];
			my $population = $temp[1];
			my $ploidy = $temp[2];
			$individuals{$individual} = $population;
			$ploidylevels{$individual} = $ploidy;
			push @{$populations{$population}}, $individual;
			if (! exists $popList{$population})
			{
				$popList{$population} = 1;
				$npops++;
			}
			print "read in $individual from $population with ploidy $ploidy in pop map file\n";
		}
		else
		{
			die ("Stopping: ploidy level (-n) $ploidyflag option not recognized. Must be either 2 or n.");
		}	
	}
}
close FH1;
##########################################################################################
##########################################################################################


##########################################################################################
#------------------------Read in coordinates or assign them -----------------------------#
##########################################################################################
my %coords = ();
my $site = 0;
if ($coordfile ne "0")
{
	if ($coordfile =~ m/\S+\.pos/)
	{
		open FH1,'<',"$coordfile";
		while(<FH1>)
		{
			if (/(\S+)\s+(\S+)/)
			{
				my $chromosome = $1;
				my $position = $2;
				$coords{chrom}{$site} = $chromosome;
				$coords{position}{$site} = $position;
				$site++;
			}
		}
		close FH1;
	}
	else
	{
		#DNA, L1.leftFlank = 1-365
		open FH1,'<',"$coordfile";
		while(<FH1>)
		{
			if (/DNA,\s*(\S+)\s*\=\s*(\d+)\-(\d+)/)
			{
				my $chromosome = $1;
				my $startPos = $2;
				my $endPos = $3;
				for my $position ($startPos..$endPos)
				{
					$coords{chrom}{$site} = $chromosome;
					$coords{position}{$site} = $position;
					$site++;
				}
			}
		}
		close FH1;
	}
}
elsif ($coordfile eq "0")
{
	my $correctLength = 0;
	open FH1,'<',"$inFile";
	while (<FH1>)
	{
		if (/^>\S+/)
		{
			#the individual is not important here
		}
		elsif (/(\S+)/)
		{
			my $seq = $1;
			if (length($seq) > $correctLength)
			{
				$correctLength = length($seq);
				for my $i (1..$correctLength)
				{
					$coords{chrom}{$site} = 1;
					$coords{position}{$site} = $i;
					$site++;
				}
			}
		}
	}	
}
##########################################################################################
##########################################################################################

##########################################################################################
#---------------------------------FASTA Processing---------------------------------------#
##########################################################################################
#separate special diploid case that can handle ambiguities
if ($ploidyflag == 2)
{
my %individualList = ();
my %seqData = ();
my $individual = "";
my @alleles = ("allele1","allele2");
open FH1,'<',"$inFile";
while(<FH1>)
{
	if (/^>(\S+)/)
	{
		$individual = $1;
#		print "$individual present in fasta\n";
		foreach my $tax (sort keys %individuals)
		{
#			print "\tchecking agaist $tax in pop map file\n";
			if (index($individual, $tax) >= 0)
			{
				if (! exists $individualList{$individual})
				{
					$individualList{$individual} = 1;
				}
			}	
		}
	}
	elsif (/(\S+)/)
	{
		my $seq = $1;
		if (exists $individualList{$individual})
		{
		my @temp = ();
		@temp = split(//,$seq);
		#print "$individual\n";
		for my $i (0..(scalar(@temp)-1))
		{
			for my $k (0..1)
			{
				#print "\t$i\t$k\t$temp[$i]\n";
				if ($k == 0)
				{
					if ($temp[$i] eq "A")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "A";
					}
					if ($temp[$i] eq "C")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "C";
					}
					if ($temp[$i] eq "G")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "G";
					}
					if ($temp[$i] eq "T")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "T";
					}
					if ($temp[$i] eq "M")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "A";
					}
					if ($temp[$i] eq "R")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "A";
					}
					if ($temp[$i] eq "W")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "A";
					}
					if ($temp[$i] eq "S")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "C";
					}
					if ($temp[$i] eq "Y")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "C";
					}
					if ($temp[$i] eq "K")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "G";
					}
					elsif ($temp[$i] ne "A" and $temp[$i] ne "C" and $temp[$i] ne "G" and $temp[$i] ne "T" and $temp[$i] ne "M" and $temp[$i] ne "R" and $temp[$i] ne "W" and $temp[$i] ne "S" and $temp[$i] ne "Y" and $temp[$i] ne "K")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "N";
					}
				}
				elsif ($k == 1)
				{
					if ($temp[$i] eq "A")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "A";
					}
					if ($temp[$i] eq "C")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "C";
					}
					if ($temp[$i] eq "G")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "G";
					}
					if ($temp[$i] eq "T")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "T";
					}
					if ($temp[$i] eq "M")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "C";
					}
					if ($temp[$i] eq "R")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "G";
					}
					if ($temp[$i] eq "W")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "T";
					}
					if ($temp[$i] eq "S")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "G";
					}
					if ($temp[$i] eq "Y")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "T";
					}
					if ($temp[$i] eq "K")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "T";
					}
					elsif ($temp[$i] ne "A" and $temp[$i] ne "C" and $temp[$i] ne "G" and $temp[$i] ne "T" and $temp[$i] ne "M" and $temp[$i] ne "R" and $temp[$i] ne "W" and $temp[$i] ne "S" and $temp[$i] ne "Y" and $temp[$i] ne "K")
					{
						$seqData{$individual}{$alleles[$k]}[$i] = "N";
					}
				}
			}
		}
		}
	}
}
##########################################################################################
##########################################################################################


##########################################################################################
#----------------------Count Alleles and Score Biallelic Sites---------------------------#
##########################################################################################
####
#Count alleles to check popmap is biallelic
####
my @biallelicSites = ();
for my $i (0..($site-1))
{
	#print OUT1 "$coords{chrom}{$i} $coords{position}{$i}";
	my @counts = ();
	for my $j (0..3)
	{
		$counts[$j] = 0;
	}
	foreach my $pop (sort keys %popList)
	{
		#print "Site $i in Population $pop\n";
		for my $j (0..(scalar(@{$populations{$pop}})-1))
		{
			#print "Counting alleles in individual $populations{$pop}[$j]\n";
			for my $k (0..1)
			{
				#print "Allele = $seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i]\n";
				if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "A")
				{
					$counts[0] = $counts[0] + 1;
				}
				if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "C")
				{
					$counts[1] = $counts[1] + 1;
				}
				if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "G")
				{
					$counts[2] = $counts[2] + 1;
				}
				if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "T")
				{
					$counts[3] = $counts[3] + 1;
				}	
			}
		}
		#print OUT1 " $counts[0],$counts[1],$counts[2],$counts[3]";
	}
	#print OUT1 "\n";
	my $biallelicCheck = 0;
	my $macCheck = 0;
	my $nsites = 0;
	for my $j (0..3)
	{
		if ($counts[$j] >= $passedArgs{mac})
		{
			$macCheck++;
		}
		if ($counts[$j] > 0)
		{
			$biallelicCheck++;
		}
		$nsites = $nsites + $counts[$j];
	}
	if ($biallelicCheck == 2 && $nsites >= $passedArgs{minsites} && $macCheck == 2)
	{
		$biallelicSites[$i] = 1;
	}
	elsif ($biallelicCheck != 2 || $nsites < $passedArgs{minsites} || $macCheck != 2 )
	{
		$biallelicSites[$i] = 0;
	}
}

my $countedSNPs = 0;
for my $i (0..($site - 1))
{
	my $locusLabelPrefix = $coords{chrom}{$i};
	$locusLabelPrefix =~ s/\_//g;
#	$locusLabelPrefix =~ s/\-//g;
	$locusLabelPrefix =~ s/\.//g;
	my $locusLabel = "$locusLabelPrefix" . "__" . "$coords{position}{$i}";
	if ($countedSNPs == 0 && $biallelicSites[$i] == 1)
	{
		print OUT1 "$locusLabel";
		$countedSNPs++;
	}
	elsif ($i > 0 && $biallelicSites[$i] == 1)
	{
		print OUT1 "\t$locusLabel";
		$countedSNPs++;
	}
}
print OUT1 "\n";
print "$countedSNPs SNPs retained afer population filtering\n";

my $popNumber = 0;
foreach my $pop (sort keys %popList)
{
	#print "Site $i in Population $pop\n";
	$popNumber++;
	for my $j (0..(scalar(@{$populations{$pop}})-1))
	{
		#print "Counting alleles in individual $populations{$pop}[$j]\n";
		for my $k (0..1)
		{
			print OUT1 "$populations{$pop}[$j]\t$popNumber";
			for my $i (0..($site - 1))
			{
				if ($biallelicSites[$i] == 1)
				{
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "A")
					{
						print OUT1 "\t1";
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "C")
					{
						print OUT1 "\t2";
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "G")
					{
						print OUT1 "\t3";
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "T")
					{
						print OUT1 "\t4";
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "N")
					{
						print OUT1 "\t-9";
					}
				}
			}
			print OUT1 "\n";
		}
	}
}
close OUT1;
}
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#arbitrary ploidy case (up to 8) assumes haploid sequences
if ($ploidyflag == 0)
{
	my %individualList = ();
	my %seqData = ();
	my $individual = "";
	my @alleles = ("allele1","allele2","allele3","allele4","allele5","allele6","allele7","allele8");
	my %haploidTracker = ();
	my $maxPloidy = 0;
	my %minorAlleles = ();
	
	open FH1,'<',"$inFile";
	while(<FH1>)
	{
		if (/^>(\S+)/)
		{
			$individual = $1;
#			print "$individual present in fasta\n";
			foreach my $tax (sort keys %individuals)
			{
#				print "\tchecking agaist $tax in pop map file\n";
				if (index($individual, $tax) >= 0)
				{
					if (! exists $haploidTracker{$tax})
					{
						$haploidTracker{$tax} = 1;
#						print "$individual from fasta found in pop map $tax\n";
					}
					elsif (exists $haploidTracker{$tax})
					{
						$haploidTracker{$tax} = $haploidTracker{$tax} + 1;
#						print "\tAND again\n";
					}
					
					if ($haploidTracker{$tax} > $maxPloidy)
					{
						$maxPloidy = $haploidTracker{$tax};
					}
					if (! exists $individualList{$individual})
					{
						$individualList{$individual} = $tax;
					}
				}
			}
		}
		elsif (/(\S+)/)
		{
			my $seq = $1;
			if (exists $individualList{$individual})
			{
			my @temp = ();
			@temp = split(//,$seq);
			#print "$individual\n";
			for my $i (0..(scalar(@temp)-1))
			{

				if ($temp[$i] eq "A")
				{
					$seqData{$individualList{$individual}}{$alleles[($haploidTracker{$individualList{$individual}} - 1)]}[$i] = "A";
				}
				if ($temp[$i] eq "C")
				{
					$seqData{$individualList{$individual}}{$alleles[($haploidTracker{$individualList{$individual}} - 1)]}[$i] = "C";
				}
				if ($temp[$i] eq "G")
				{
					$seqData{$individualList{$individual}}{$alleles[($haploidTracker{$individualList{$individual}} - 1)]}[$i] = "G";
				}
				if ($temp[$i] eq "T")
				{
					$seqData{$individualList{$individual}}{$alleles[($haploidTracker{$individualList{$individual}} - 1)]}[$i] = "T";
				}
				elsif ($temp[$i] ne "A" and $temp[$i] ne "C" and $temp[$i] ne "G" and $temp[$i] ne "T")
				{
					$seqData{$individualList{$individual}}{$alleles[($haploidTracker{$individualList{$individual}} - 1)]}[$i] = "N";
				}
			}
			}
		}
	}

##########################################################################################
##########################################################################################


##########################################################################################
#----------------------Count Alleles and Score Biallelic Sites---------------------------#
##########################################################################################
####
#Count alleles to check popmap is biallelic
####
my @biallelicSites = ();
for my $i (0..($site-1))
{
	#print OUT1 "$coords{chrom}{$i} $coords{position}{$i}";
	my @counts = ();
	for my $j (0..3)
	{
		$counts[$j] = 0;
	}
	foreach my $pop (sort keys %popList)
	{
		#print "Site $i in Population $pop\n";
		for my $j (0..(scalar(@{$populations{$pop}})-1))
		{
			#print "Counting alleles in individual $populations{$pop}[$j]\n";
			for my $k (0..($ploidylevels{$populations{$pop}[$j]} - 1))
			{
				#print "Allele = $seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i]\n";
				if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "A")
				{
					$counts[0] = $counts[0] + 1;
				}
				if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "C")
				{
					$counts[1] = $counts[1] + 1;
				}
				if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "G")
				{
					$counts[2] = $counts[2] + 1;
				}
				if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "T")
				{
					$counts[3] = $counts[3] + 1;
				}	
			}
		}
		#print OUT1 " $counts[0],$counts[1],$counts[2],$counts[3]";
	}
	#print OUT1 "\n";
	my $biallelicCheck = 0;
	my $macCheck = 0;
	my $minorAllele = 0;
	my $majorAllele = 0;
	my $nsites = 0;
	for my $j (0..3)
	{
		if ($counts[$j] >= $passedArgs{mac})
		{
			$macCheck++;
			if ($counts[$j] > $majorAllele)
			{
				if ($majorAllele != 0)
				{
					$minorAlleles{minor}[$i] = $minorAlleles{major}[$i];
				}
				$minorAlleles{major}[$i] = $j;
				$majorAllele = $counts[$j];
			}
			elsif ($counts[$j] <= $majorAllele)
			{
				$minorAlleles{minor}[$i] = $j;
			}
			else
			{
				print "Serious problem, allele counts not correct!\n"
			}	
		}
		if ($counts[$j] > 0)
		{
			$biallelicCheck++;
		}
		$nsites = $nsites + $counts[$j];
	}
	if ($biallelicCheck == 2 && $nsites >= $passedArgs{minsites} && $macCheck == 2)
	{
		$biallelicSites[$i] = 1;
	}
	elsif ($biallelicCheck != 2 || $nsites < $passedArgs{minsites} || $macCheck != 2)
	{
		$biallelicSites[$i] = 0;
	}
}

my $countedSNPs = 0;
print OUT2 "Individual\t";
print OUT3 "Individual\t";
print OUT4 "Individual\t";
print OUT5 "Individual\t";
for my $i (0..($site - 1))
{
	my $locusLabelPrefix = $coords{chrom}{$i};
	$locusLabelPrefix =~ s/\_//g;
#	$locusLabelPrefix =~ s/\-//g;
	$locusLabelPrefix =~ s/\.//g;
	my $locusLabel = "$locusLabelPrefix" . "__" . "$coords{position}{$i}";
	if ($countedSNPs == 0 && $biallelicSites[$i] == 1)
	{
		print OUT1 "$locusLabel";
		print OUT2 "$locusLabel";
		print OUT3 "$locusLabel";
		print OUT4 "$locusLabel";
		print OUT5 "$locusLabel";
		$countedSNPs++;
	}
	elsif ($i > 0 && $biallelicSites[$i] == 1)
	{
		print OUT1 "\t$locusLabel";
		print OUT2 "\t$locusLabel";
		print OUT3 "\t$locusLabel";
		print OUT4 "\t$locusLabel";
		print OUT5 "\t$locusLabel";
		$countedSNPs++;
	}
}
print OUT1 "\n";
print OUT2 "\n";
print OUT3 "\n";
print OUT4 "\n";
print OUT5 "\n";
print "$countedSNPs SNPs retained afer population filtering\n";

my $popNumber = 0;
foreach my $pop (sort keys %popList)
{
	#print "Site $i in Population $pop\n";
	$popNumber++;
	for my $j (0..(scalar(@{$populations{$pop}})-1))
	{
		#print "Counting alleles in individual $populations{$pop}[$j]\n";
		print OUT2 "$populations{$pop}[$j]";
		print OUT3 "$populations{$pop}[$j]";
		print OUT4 "$populations{$pop}[$j]";
		print OUT5 "$populations{$pop}[$j]";
################
#Output structure format
################
		for my $k (0..($ploidylevels{$populations{$pop}[$j]} - 1))
		{
			print OUT1 "$populations{$pop}[$j]\t$popNumber";
			for my $i (0..($site - 1))
			{
				if ($biallelicSites[$i] == 1)
				{
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "A")
					{
						print OUT1 "\t1";
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "C")
					{
						print OUT1 "\t2";
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "G")
					{
						print OUT1 "\t3";
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "T")
					{
						print OUT1 "\t4";
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "N")
					{
						print OUT1 "\t-9";
					}
				}
			}
			print OUT1 "\n";
		}
################
#Output pca matrices
################
#Missing data would need to be imputed for downstream analyses. There are some good EM methods to do this in R or just let adegenet handle it.
		for my $i (0..($site - 1))
		{
			if ($biallelicSites[$i] == 1)
			{
				my $pcaRaw = 0;
				my $pcaCompressed5 = 0;
				my $pcaCompressed3 = 0;
				my $pcaCompressed2 = 0;
				my $ismissing = 0;
				for my $k (0..($ploidylevels{$populations{$pop}[$j]} - 1))
				{	
#					print "$pop\t$populations{$pop}[$j]\t$i\t$k\t$minorAlleles{minor}[$i]\t:\t$pcaRaw ==>";
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "A")
					{
						if ($minorAlleles{minor}[$i] == 0 && $ismissing == 0)
						{
							$pcaRaw++;
						}
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "C")
					{
						if ($minorAlleles{minor}[$i] == 1 && $ismissing == 0)
						{
							$pcaRaw++;
						}
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "G")
					{
						if ($minorAlleles{minor}[$i] == 2 && $ismissing == 0)
						{
							$pcaRaw++;
						}
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "T")
					{
						if ($minorAlleles{minor}[$i] == 3 && $ismissing == 0)
						{
							$pcaRaw++;
						}
					}
					if ($seqData{$populations{$pop}[$j]}{$alleles[$k]}[$i] eq "N")
					{
						$ismissing = 1;
						$pcaRaw = "NA";
						$pcaCompressed5 = "NA";
						$pcaCompressed3 = "NA";
						$pcaCompressed2 = "NA";
					}
#					print "\t$pcaRaw\n";
				}

		
				if ($ismissing == 0)
				{
					if ($pcaRaw == 0)
					{
						$pcaCompressed5 = 0;
						$pcaCompressed3 = 0;
						$pcaCompressed2 = 0;
					}
					elsif (($pcaRaw/$ploidylevels{$populations{$pop}[$j]}) == 1)
					{
						$pcaCompressed5 = 4;
						$pcaCompressed3 = 2;
						$pcaCompressed2 = 1;
					}
					elsif (($pcaRaw/$ploidylevels{$populations{$pop}[$j]}) == 0.5)
					{
						$pcaCompressed5 = 2;
						$pcaCompressed3 = 1;
						$pcaCompressed2 = 1;
					}
					elsif ((($pcaRaw/$ploidylevels{$populations{$pop}[$j]}) > 0.5 && ($pcaRaw/$ploidylevels{$populations{$pop}[$j]}) != 1))
					{
						$pcaCompressed5 = 3;
						$pcaCompressed3 = 1;
						$pcaCompressed2 = 1;
					}
					elsif ((($pcaRaw/$ploidylevels{$populations{$pop}[$j]}) < 0.5 && ($pcaRaw/$ploidylevels{$populations{$pop}[$j]}) != 0))
					{
						$pcaCompressed5 = 1;
						$pcaCompressed3 = 1;
						$pcaCompressed2 = 1;
					}
				}
				elsif ($ismissing == 1)
				{
#					print "missing site $i for $populations{$pop}[$j] - santiy check\n";
				}
				print OUT2 "\t$pcaRaw";
				print OUT3 "\t$pcaCompressed5";
				print OUT4 "\t$pcaCompressed3";
				print OUT5 "\t$pcaCompressed2";
			}	
		}
		print OUT2 "\n";
		print OUT3 "\n";
		print OUT4 "\n";
		print OUT5 "\n";
	}
}

}
close OUT1;
close OUT2;
close OUT3;
close OUT4;
close OUT5;

print "Note that some locus labels may have changed in output files to avoid conflict with downstream interpretation of special characters\nPlease confirm labeling has not caused a loss of information.\n";

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
exit;