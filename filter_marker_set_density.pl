#! /usr/bin/perl

use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Std;
use IO::File;
my % options=();
use IO::Handle;

my$opt_string='h:o:m:n:';
getopts("$opt_string",\%options) or usage();
usage()if $options{h};


my $marker_set="";
my $min_number_of_markers=100;
my $output_file="";

init();
main();
sub usage
{
	my $usage= "$0  -m parental_marker_set -n minimum_marker_set -o output\n\n" ;
	print $usage;
	exit();
}



sub init
{
	if($options{m})
	{
		$marker_set = $options{m};
	}
	else
	{
		die usage();
	}
	
	if($options{n})
	{
		$min_number_of_markers = $options{n};
	}
	
	if($options{o})
	{
		$output_file=$options{o};
	}
}



sub main 
{
	my %sample = ();
	
	my $inputFile = new IO::File($marker_set, "r") or die "could not open $marker_set: $!\n";
		
	while(my $line=$inputFile->getline)
	{
		chomp($line);
		my @a = split " ",$line;
		
		if(!exists $sample{$a[0]})
		{
			my @temp =($line);
			$sample{$a[0]} = [@temp];
		}
		else
		{
			push(@{$sample{$a[0]}},$line);	
		}
	}
	my $file_writer=IO::File->new();
	$file_writer->open(">".$output_file) or die "Could not open ".$output_file.".txt: $!\n" ;
	
	
	my @sorted_chr = sort {$a <=> $b} keys %sample;
	
	foreach my $chr (@sorted_chr)
	{
		if(@{$sample{$chr}} >= $min_number_of_markers)
		{
			foreach my $line (@{$sample{$chr}})
			{
				$file_writer->say($line);
			}			
		}
	}
}


