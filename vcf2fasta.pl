#!/usr/bin/perl -w
use strict;
#----------------------------------------------------------------------------------------#
#George P. Tiley
#3 February 2021
#contact: george.tiley@duke.edu
#Convert VCF file to a fasta file
#----------------------------------------------------------------------------------------#


##########################################################################################
#--------------------------------Input Flag Handling-------------------------------------#
##########################################################################################
my @checkArgs = ("i","o","t","d","b");
my %passedArgs = ();
if (scalar(@ARGV) == 0)
{
	die "/*--------INPUT PARAMETERS--------*/\n
	-i STRING <input vcf file>
	-o STRING <output fasta file>
	-t INT <flag for SNP thinning strategy> 
	-d INT <distance for thinning SNPS>
	-b INT <biallelic SNPS only? 0 = no || 1 = yes>
	\n/*--------EXAMPLE COMMAND--------*/\n
	./vcf2fasta.pl -i murinusClade.keepBad1.mac3.FS6.vcf -o Mmur_all.fasta -t 2 -d 10000 -b 1\n";
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
		if ($ARGV[$i] eq "-t")
		{
			$passedArgs{t} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-d")
		{
			$passedArgs{d} = $ARGV[$i+1];
		}
		if ($ARGV[$i] eq "-b")
		{
			$passedArgs{b} = $ARGV[$i+1];
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
my $outFile = $passedArgs{o};
my $thinning = $passedArgs{t};
my $interval = $passedArgs{d};
my $biallelic = $passedArgs{b};
##########################################################################################
##########################################################################################

##########################################################################################
#----------------------------------VCF Processesing--------------------------------------#
##########################################################################################
open FH1, '<', "$inFile";

my %metaData = ();
my %snps = ();
my %taxList = ();

while (<FH1>)
{
	my $line = $_;
	chomp $line;
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mgan001	mgan002	mgan004	mgan005	mgan006	mgan007	mgan008	mgan009	mgan010	mgan011	mgan014	mgan016	mgan017	mgan018	mgan019	mgan021	mgan022	mman003	mman004	mman005	mman006	mmur001	mmur002	mmur003	mmur004	mmur006	mmur009	mmur010	mmur011	mmur012	mmur013	mmur014	mmur016	mmur019	mmur030	mmur031	mmur032	mmur033	mmur037	mmur038	mmur039	mmur040	mmur041	mmur042	mmur043	mmur044	mmur045	mmur046	mmur047	mmur048	mmur049	mmur050	mmur051	mmur052	mmur053	mmur055	mmur056	mmur057	mmur058	mmur061	mmur063	mmur064	mmur065	mmur066	mmur068	mmur069	mmur070	mmur071	mmur072	mruf007
#NC_033660.1	3539628	.	G	A	2717.49	PASS	ABHet=0.57;ABHom=1;AC=9;AF=0.051;AN=176;BaseQRankSum=0;DP=2071;ExcessHet=0.0267;FS=0;InbreedingCoeff=0.3845;MLEAC=11;MLEAF=0.063;MQ=60;MQRankSum=0;OND=0;QD=12.02;ReadPosRankSum=0;SOR=0.732	GT:AD:DP:GQ:PL	0/0:20,0:20:60:0,60,630	0/0:37,0:37:99:0,102,1530	0/0:38,0:38:99:0,102,1292	0/0:34,0:34:99:0,102,1021	0/0:20,0:20:43:0,43,481	0/0:13,0:13:39:0,39,469	0/0:35,0:35:99:0,102,1150	0/0:18,0:18:54:0,54,682	0/0:39,0:39:99:0,101,1291	0/0:39,0:39:83:0,83,1290	0/0:26,0:26:61:0,61,856	0/0:20,0:20:43:0,43,537	0/0:18,0:18:54:0,54,605	0/0:20,0:20:60:0,60,609	0/0:31,0:31:93:0,93,1008	0/0:35,0:35:96:0,96,1440	0/0:15,0:15:45:0,45,396	0/0:49,0:49:99:0,111,1784	0/1:25,5:30:54:54,0,619	0/0:27,0:27:64:0,64,883	0/0:32,0:32:93:0,93,1063	0/0:28,0:28:84:0,84,947	0/0:24,0:24:72:0,72,840	0/0:32,0:32:90:0,90,1149	0/0:22,0:22:66:0,66,674	0/0:24,0:24:55:0,55,638	0/0:35,0:35:99:0,105,1204	0/0:12,0:12:36:0,36,451	0/0:63,0:63:99:0,120,1800	0/0:20,0:20:60:0,60,659	0/0:10,0:10:30:0,30,269	0/0:37,0:37:99:0,106,1300	0/0:25,0:25:72:0,72,793	0/0:41,0:41:99:0,106,1348	./.:2,0:2:6:0,6,49	0/0:14,0:14:42:0,42,565	0/0:8,0:8:24:0,24,347	./.:0,0:0:.:0,0,0	0/0:47,0:47:99:0,120,1800	0/1:17,18:35:99:385,0,402	0/0:26,0:26:78:0,78,945	1/1:0,23:23:69:742,69,0	0/0:16,0:16:48:0,48,587	0/0:7,0:7:21:0,21,236	0/0:16,0:16:48:0,48,423	0/0:41,0:41:99:0,116,1582	0/0:23,0:23:69:0,69,646	0/0:25,0:25:75:0,75,854	0/0:52,0:52:99:0,120,1800	0/0:14,0:14:42:0,42,380	0/0:25,0:25:72:0,72,1080	0/0:28,0:28:84:0,84,926	0/0:24,0:24:72:0,72,669	0/1:11,35:46:99:865,0,194	1/1:0,10:10:30:299,30,0	0/0:20,0:20:60:0,60,576	0/0:14,0:14:42:0,42,548	0/0:25,0:25:72:0,72,982	0/0:8,0:8:24:0,24,190	0/0:29,0:29:84:0,84,1260	0/0:29,0:29:87:0,87,871	0/0:6,0:6:18:0,18,148	0/0:25,0:25:72:0,72,1080	0/0:14,0:14:42:0,42,338	0/1:30,9:39:99:213,0,904	0/0:28,0:28:84:0,84,964	0/0:30,0:30:90:0,90,1105	0/0:7,0:7:21:0,21,277	0/0:35,0:35:99:0,105,1467	0/0:52,0:52:99:0,105,1648

#skip all of the comment lines 
	if ($line !~ m/^#.+/)
	{
		my @temp = ();
		@temp = split(/\t/, $line);
		my $allele1 = $temp[3];
		my $allele2 = $temp[4];
		
		my $chrom = $temp[0];
		my $position = $temp[1];
		
		my $alleleCount1 = 0;
		my $alleleCount2 = 0;
		
		
		for my $i (9..(scalar(@temp)-1))
		{
			my $consensus = "N";
			if ($temp[$i] =~ m/(\d+)\/(\d+)\:\S+/)
			{
				my $a1 = $1;
				my $a2 = $2;
				if ($a1 == 0 && $a2 == 0)
				{
					$consensus = $allele1;
					$alleleCount1 = $alleleCount1 + 2;
				}
				elsif ($a1 == 1 && $a2 == 1)
				{
					$consensus = $allele2;
					$alleleCount2 = $alleleCount2 + 2;
				}
				elsif (($a1 == 0 && $a2 == 1) || ($a1 == 1 && $a2 == 0))
				{
					$consensus = resolveHeterozygote($allele1,$allele2);
					$alleleCount1 = $alleleCount1 + 1;
					$alleleCount2 = $alleleCount2 + 1;
				}
				else
				{
					$consensus = "N";
				}
				
			}
			elsif ($temp[$i] =~ m/\.\/\.\:\S+/)
			{
				$consensus = "N";
			}
			push @{$snps{$taxList{$i}}}, $consensus;
		}
		push @{$metaData{chromosome}}, $chrom;
		push @{$metaData{position}}, $position;
		if ($alleleCount1 >= 1 && $alleleCount2 >= 1)
		{
			push @{$metaData{biallelic}}, 1;
		}
		elsif ($alleleCount1 == 0 || $alleleCount2 == 0)
		{
			push @{$metaData{biallelic}}, 0;
		}
	}
	#grab the header line
	elsif ($line =~ m/^#CHROM.+/)
	{
		my @temp = ();
		@temp = split(/\t/, $line);
		for my $i (9..(scalar(@temp)-1))
		{
			if (! exists $taxList{$i})
			{
#				$snps{$temp[$i]} = 1;
				$taxList{$i} = $temp[$i];
			}
		}
	}
}
####
#Dump SNPS out to fasta
#thinning option is greedy and makes SNPS be some interval apart
open OUT1,'>',"$outFile";
open OUT2,'>',"$outFile.pos";
my $posfile = 0;


if ($thinning == 0)
{
	foreach my $taxonNumber (sort keys %taxList)
	{
		print OUT1 ">$taxList{$taxonNumber}\n";
		for my $i (0..(scalar(@{$snps{$taxList{$taxonNumber}}})-1))
		{
			if ($biallelic == 0)
			{
				print OUT1 "$snps{$taxList{$taxonNumber}}[$i]";
				if ($posfile == 0)
				{
					print OUT2 "$metaData{chromosome}[$i]\t$metaData{position}[$i]\n";
				}
				
			}
			elsif ($biallelic == 1)
			{
				if ($metaData{biallelic}[$i] == 1)
				{
					print OUT1 "$snps{$taxList{$taxonNumber}}[$i]";
					if ($posfile == 0)
					{
						print OUT2 "$metaData{chromosome}[$i]\t$metaData{position}[$i]\n";
					}
				}
			}
		}
		print OUT1 "\n";
		$posfile = 1;
	}
}

elsif ($thinning == 1)
{
	my $thinnedSNPS = 0;
	my %snpChecker = ();
	foreach my $taxonNumber (sort keys %taxList)
	{
		print OUT1 ">$taxList{$taxonNumber}\n";
		if ($biallelic == 0)
		{
			for my $i (0..(scalar(@{$snps{$taxList{$taxonNumber}}})-1))
			{
				my $check = 0;
				if ($i > 0)
				{
					if ($metaData{chromosome}[$i] eq $metaData{chromosome}[$i-1])
					{
						print "$metaData{chromosome}[$i]\t$metaData{chromosome}[$i-1]\n";
						my $distance = $metaData{position}[$i] - $metaData{position}[$i-1];
						if ($distance < $interval)
						{
							$check = 1;
							print "==> linked sites at $distance bp distance\n";
						}	
					}
				}
				if ($check == 0)
				{
					print OUT1 "$snps{$taxList{$taxonNumber}}[$i]";
					if (! exists $snpChecker{$i})
					{
						$thinnedSNPS++;
						$snpChecker{$i} = 1;
					}
					if ($posfile == 0)
					{
						print OUT2 "$metaData{chromosome}[$i]\t$metaData{position}[$i]\n";
					}
				}
			}
		}
		elsif ($biallelic == 1)
		{
			my $skipped = 1;
			for my $i (0..(scalar(@{$snps{$taxList{$taxonNumber}}})-1))
			{
				my $check = 0;
				if ($metaData{biallelic}[$i] == 0)
				{	
					$skipped++;
				}
				if ($metaData{biallelic}[$i] == 1)
				{
					if ($i > 0)
					{
						if ($metaData{chromosome}[$i] eq $metaData{chromosome}[$i-$skipped])
						{
							print "$metaData{chromosome}[$i]\t$metaData{chromosome}[$i-$skipped]\n";
							my $distance = $metaData{position}[$i] - $metaData{position}[$i-$skipped];
							if ($distance < $interval)
							{
								$check = 1;
								print "==> linked sites at $distance bp distance\n";
							}	
						}
					}
					if ($check == 0)
					{
						print OUT1 "$snps{$taxList{$taxonNumber}}[$i]";
						if (! exists $snpChecker{$i})
						{
							$thinnedSNPS++;
							$snpChecker{$i} = 1;
						}
						if ($posfile == 0)
						{
							print OUT2 "$metaData{chromosome}[$i]\t$metaData{position}[$i]\n";
						}
					}
					$skipped = 1;
				}
			}
		}
		print OUT1 "\n";
		$posfile = 1;
	}
	print "$thinnedSNPS SNPs recovered after thinning at greedy $interval bp intervals\n";
}

elsif ($thinning == 2)
{
	my $thinnedSNPS = 0;
	my %snpChecker = ();
	my $ntax = 0;
	my $currentSNP = 0;
	
	foreach my $taxonNumber (sort keys %taxList)
	{
		#print OUT1 ">$taxList{$taxonNumber}\n";
		if ($biallelic == 0)
		{
			if ($ntax == 0)
			{
				for my $i (0..(scalar(@{$snps{$taxList{$taxonNumber}}})-1))
				{				
					if ($currentSNP != $i)
					{
						if ($metaData{chromosome}[$currentSNP] eq $metaData{chromosome}[$i])
						{		
							print "$metaData{chromosome}[$i]\t$metaData{chromosome}[$i-1]\n";
							my $distance = $metaData{position}[$i] - $metaData{position}[$currentSNP];
							if ($distance < $interval)
							{
								my $pSurvival = rand(1);
								if ($pSurvival < 0.5)
								{
									$currentSNP = $currentSNP;
								}
								elsif ($pSurvival >= 0.5)
								{
									$currentSNP = $i;
								}
								print "==> linked sites at $distance bp distance\n";
							}
							elsif ($distance >= $interval)
							{
								$snpChecker{$currentSNP} = 1;
								$currentSNP = $i;
								$thinnedSNPS++;
								
								if ($i == (scalar(@{$snps{$taxList{$taxonNumber}}})-1))
								{
									$snpChecker{$currentSNP} = 1;
									$thinnedSNPS++;
								}
								print "==> unlinked sites at $distance bp distance\n";
							}
						}
						elsif ($metaData{chromosome}[$currentSNP] ne $metaData{chromosome}[$i])
						{
							print "$metaData{chromosome}[$i]\t$metaData{chromosome}[$i-1]\n";
							print "==> unlinked sites on different chromosomes\n";
							$snpChecker{$currentSNP} = 1;
							$currentSNP = $i;
							$thinnedSNPS++;
							if ($i == (scalar(@{$snps{$taxList{$taxonNumber}}})-1))
							{
								$snpChecker{$currentSNP} = 1;
								$thinnedSNPS++;
							}
						}
					}
				}
				$ntax++;
			}
		}
		elsif ($biallelic == 1)
		{	
			if ($ntax == 0)
			{
				my $startSNP = 0;
				my $checkStart = 0;
				while ($checkStart == 0)
				{
					if ($metaData{biallelic}[$startSNP] == 0)
					{
						$startSNP++;
					}
					elsif ($metaData{biallelic}[$startSNP] == 1)
					{
						$checkStart = 1;
						$currentSNP = $startSNP;
					}
				}
				
				for my $i ($startSNP..(scalar(@{$snps{$taxList{$taxonNumber}}})-1))
				{	
					if ($metaData{biallelic}[$i] == 1)
					{		
						if ($currentSNP != $i)
						{
							if ($metaData{chromosome}[$currentSNP] eq $metaData{chromosome}[$i])
							{			
								#print "$metaData{chromosome}[$i]\t$metaData{chromosome}[$i-1]\n";
								my $distance = $metaData{position}[$i] - $metaData{position}[$currentSNP];
								if ($distance < $interval)
								{
									my $pSurvival = rand(1);
									if ($pSurvival < 0.5)
									{
										$currentSNP = $currentSNP;
									}
									elsif ($pSurvival >= 0.5)
									{
										$currentSNP = $i;
									}
									print "==> linked sites at $distance bp distance\n";
								}
								elsif ($distance >= $interval)
								{
									$snpChecker{$currentSNP} = 1;
									$currentSNP = $i;
									$thinnedSNPS++;
									if ($i == (scalar(@{$snps{$taxList{$taxonNumber}}})-1))
									{
										if ($metaData{biallelic}[$i] == 1)
										{
											$snpChecker{$currentSNP} = 1;
											$thinnedSNPS++;
										}
									}
									print "==> unlinked sites at $distance bp distance\n";
								}
							}
							elsif ($metaData{chromosome}[$currentSNP] ne $metaData{chromosome}[$i])
							{
								print "$metaData{chromosome}[$i]\t$metaData{chromosome}[$i-1]\n";
								print "==> unlinked sites on different chromosomes\n";
								$snpChecker{$currentSNP} = 1;
								$currentSNP = $i;
								$thinnedSNPS++;
								if ($i == (scalar(@{$snps{$taxList{$taxonNumber}}})-1))
								{
									if ($metaData{biallelic}[$i] == 1)
									{
										$snpChecker{$currentSNP} = 1;
										$thinnedSNPS++;
									}
								}
							}
						}
					}
					$ntax++;
				}
			}
		}
	}

	foreach my $taxonNumber (sort keys %taxList)
	{	
		print OUT1 ">$taxList{$taxonNumber}\n";
		foreach my $unlinkedSNP (sort {$a <=> $b} keys %snpChecker)
		{
			print OUT1 "$snps{$taxList{$taxonNumber}}[$unlinkedSNP]";
			if ($posfile == 0)
			{
				print OUT2 "$metaData{chromosome}[$unlinkedSNP]\t$metaData{position}[$unlinkedSNP]\n";
			}
		}
		print OUT1 "\n";
		$posfile = 1;
	}
	print "$thinnedSNPS SNPs recovered after thinning at $interval bp intervals with random SNP selection\n";
}
close OUT2;
##########################################################################################
##########################################################################################

##########################################################################################
#-------------------------------------Subroutines----------------------------------------#
##########################################################################################

sub resolveHeterozygote
{
	my $bp1 = $_[0];
	my $bp2 = $_[1];
	#initialize to N in case something weird happens
	my $resolvedConsensus = "N";
	
	#AC AG AT CG CT GT
	
	if (($bp1 eq "A" and $bp2 eq "C") or ($bp1 eq "C" and $bp2 eq "A"))
	{
		$resolvedConsensus = "M";
	}
	elsif (($bp1 eq "A" and $bp2 eq "G") or ($bp1 eq "G" and $bp2 eq "A"))
	{
		$resolvedConsensus = "R";
	}
	if (($bp1 eq "A" and $bp2 eq "T") or ($bp1 eq "T" and $bp2 eq "A"))
	{
		$resolvedConsensus = "W";
	}
	if (($bp1 eq "C" and $bp2 eq "G") or ($bp1 eq "G" and $bp2 eq "C"))
	{
		$resolvedConsensus = "S";
	}
	if (($bp1 eq "C" and $bp2 eq "T") or ($bp1 eq "T" and $bp2 eq "C"))
	{
		$resolvedConsensus = "Y";
	}
	if (($bp1 eq "G" and $bp2 eq "T") or ($bp1 eq "T" and $bp2 eq "G"))
	{
		$resolvedConsensus = "K";
	}

	return $resolvedConsensus;
}
##########################################################################################
##########################################################################################

exit;