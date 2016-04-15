#! /usr/bin/perl

use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Std;
use IO::File;
use Time::HiRes;
#use Constant_Math_lib;
my % options=();
use IO::Handle;


my $break_file="";
my $output_file="output.txt";
my %POS=();

my$opt_string='hb:o:';
getopts("$opt_string",\%options) or usage();
usage()if $options{h};




init();
main();
sub usage
{
	my $usage = " usage:perl prepare_data_for_overview_plot.pl -b breakfile -o output \n";
	print $usage;
	exit();
}

sub init
{
	if($options{b})
	{
		$break_file= $options{b};
	}
	else
	{
		die usage();
	}
	if($options{o})
	{
		$output_file=$options{o};
	}
}

sub main
{
	my $input_file = new IO::File($break_file, "r") or die "could not open $break_file: $!\n";
	
	my $sample="";
	my $chr=0;
	my $begin="";
	my $end="";
	my $genotype="";
	my $file_writer=IO::File->new();
	$file_writer->open(">".$output_file) or die "Could not open ".$output_file." $!\n" ;
	
	while (my $line = $input_file->getline) 
	{
		chomp($line);
		#print " read line ".$line ."\n";
     	my   @a = split " ", $line;
        if($chr != $a[1])
        {
        	#print
        	if($chr !=0)
        	{
        		$file_writer->print($sample."\t".$chr."\t".$begin."\t".$end."\t".$genotype."\n");
        	}
        	$sample=$a[0];
		    $chr=$a[1];
			$begin=$a[2];
			$end=$a[3];
			$genotype=$a[4];
        }
        else
        {
        	if($genotype eq $a[4])
        	{
        		$end=$a[3];
        	}
        	else
        	{
        		#print $sample."\t".$a[1]."\t".$a[2]."\t".$a[3]."\t".$a[4]." current line".$line ."\n";
        		$file_writer->print($sample."\t".$chr."\t".$begin."\t".$end."\t".$genotype."\n");
        		$chr=$a[1];
				$begin=$a[2];
				$end=$a[3];
				$genotype=$a[4];
        	}
        }
	}	
	$file_writer->print($sample."\t".$chr."\t".$begin."\t".$end."\t".$genotype."\n");
	
	
}
