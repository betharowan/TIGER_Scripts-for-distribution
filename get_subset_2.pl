#! /usr/bin/perl

use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Std;
use IO::File;
my % options=();
use IO::Handle;

my$opt_string='h:o:c:p:';
getopts("$opt_string",\%options) or usage();
usage()if $options{h};


my $complet_set="";
my $ref_set="";
my $output_file="";

init();
main();
sub usage
{
	my $usage= "\n$0  -c complete_marker_set_file  -p ref_marker_set -o output\n\n" ;
	print $usage;
	exit();
}



sub init
{
	if($options{c} && $options{p})
	{
		$complet_set = $options{c};
		$ref_set = $options{p};
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


sub read_completSet
{
	
	my ($complet_set_file,$Pos_hash_ref)=@_;
	my $inputFile = new IO::File($complet_set_file, "r") or die "could not open $complet_set_file: $!\n";
		
	while(my $line=$inputFile->getline)
	{
		chomp($line);
		
		my @a= split " ",$line;
		
		$$Pos_hash_ref{$a[0]."\t".$a[1]} = $line;
	}
}



sub getSelectedMarkers
{
	my ($refSet_file,$Pos_hash_ref)=@_;
	
	my $inputFile = new IO::File($refSet_file, "r") or die "could not open $refSet_file: $!\n";
	
	my $file_writer=IO::File->new();
	$file_writer->open(">".$output_file) or die "Could not open ".$output_file.".txt: $!\n" ;
	
	
	while(my $line=$inputFile->getline)
	{
		chomp($line);
		my @a= split " ",$line;
		
		if(exists $$Pos_hash_ref{$a[0]."\t".$a[1]})
		{
			$file_writer->say($$Pos_hash_ref{$a[0]."\t".$a[1]});
		}
		else
		{
			if(@a==2)
			{
				$file_writer->say($a[0]."\t".$a[1]."\t\.\t0\t\.\t0");
			}
			else
			{
				$file_writer->say($a[0]."\t".$a[1]."\t".$a[3]."\t0".$a[4]."\t0");
			}
		}
	}
}

sub main 
{
	my %POS = ();
	read_completSet($complet_set,\%POS);
	getSelectedMarkers($ref_set,\%POS);	
}


