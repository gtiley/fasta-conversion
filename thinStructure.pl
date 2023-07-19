#!/usr/bin/perl -w

$prefix = $ARGV[0];

@snps = ();
%snpBlocks = ();
%chrList = ();
$header = 1;

open FH1,'<',"$prefix.str";
open OUT1,'>',"$prefix.str.thin";
while (<FH1>)
{
	$line = $_;
	chomp $line;
	@temp = ();
	@temp = split(/\t/,$line);
	if ($header==1)
	{
		for $i (0..(scalar(@temp)-1))
		{
			if ($temp[$i] =~ m/(\S+)\_\_\d+/)
			{
				
				$chr = $1;
#				$thisPos = $2;
				push @{$snpBlocks{$chr}}, $i;
				if (! exists $chrList{$chr})
				{
					$chrList{$chr} = 1;
				}
			}
		}
		
		foreach $chr (keys %chrList)
		{
			$r = int(rand(scalar(@{$snpBlocks{$chr}})));
			push @snps, $snpBlocks{$chr}[$r];
		}
		
		$header = 0;
		
		for $i (0..(scalar(@snps)-1))
		{
			if ($i == 0)
			{
				print OUT1 "$temp[$snps[$i]]";
			}
			elsif ($i > 0)
			{
				print OUT1 "\t$temp[$snps[$i]]";
			}
		}
		print OUT1 "\n";
		$nsnps = scalar(@snps);
		print "$nsnps snps surviving after thinning\n"
	}
	elsif ($header == 0)
	{
		print OUT1 "$temp[0]\t$temp[1]";
		for $i (0..(scalar(@snps)-1))
		{
			print OUT1 "\t$temp[$snps[$i] + 2]"
		}
		print OUT1 "\n"
	}
}
close OUT1;
close FH1;
exit;
