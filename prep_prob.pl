use strict;
use warnings;
use Class::Struct;
use Getopt::Std;
use IO::File;
use Time::HiRes;
use Cwd;
my % options=();
use IO::Handle;


my$opt_string='hb:m:o:c:s:';

getopts("$opt_string",\%options) or usage();
usage()if $options{h};

my $out_file_handler = IO::File->new();
my $basecalled="";
my $marker_file="";
my $chrsizes_file="";
my $sampleID="";
my $minimum_qualitiy_value=0;
my $output_file="output.txt";
my %parent_container=();
my %chrsizes=();
my $support=0;
init();
main();
print "done \n";



sub usage
{
	print "perl $0 -m allele_count -b basecalled -c chrsizes.txt -s sampleID -o output\n";
}

sub init
{
	if($options{b} && $options{m} && $options{c} && $options{s})
	{
		$basecalled=$options{b};
		$marker_file=$options{m};
		$chrsizes_file=$options{c};
		$sampleID=$options{s};
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
#chr,pos,base call,hmm out, hmm out masked,count a, count c, count g, count t

sub main
{

	open FILE, $chrsizes_file or die "cannot open file $chrsizes_file\n";
	while (<FILE>) {
		my @a = split " ";
		$chrsizes{$a[0]} = $a[1];
	}

	my $input_file = new IO::File($marker_file, "r") or die "could not open $marker_file: $!\n";	
	my %marker_container=();
	
	my $file_writer=IO::File->new();
	my $file_writer2=IO::File->new();
	
	$file_writer->open(">".$output_file) or die "Could not open ".$output_file." $!\n" ;
	if (substr($output_file, length($output_file)-4, 4) eq ".txt") {
		$output_file = substr($output_file, 0, length($output_file)-4);
	}
#	$file_writer2->open(">".$output_file.".breaks.txt") or die "Could not open output gfile\n" ;
	
	while(my $l = $input_file->getline)
	{
		if(substr($l,0,1)eq "#")
		{
			next;
		}
		my @a = split " ",$l;
		if(!exists $marker_container{$a[0]})
		{
			my @temp=($l);
			$marker_container{$a[0]}=[@temp];
			
		}
		else
		{
			push @{ $marker_container{$a[0]}},$l;
		}
				
	}
	#my @chr_key = sort {$a<=>$b} keys %marker_container;
	$input_file = new IO::File($basecalled, "r") or die "could not open $basecalled: $!\n";
	while(my $l = $input_file->getline)
	{
		if(substr($l,0,1)eq "#")
		{
			#print $l."\n";
			chomp($l);
			my @a = split ":",$l;
			my $chr = $a[2];
			#print $chr."\n";
			my @basecaller=split " ",$input_file->getline; #basecaller
			#my @hmm=split " ",$input_file->getline; #hmm
			#my @hmm_masked=split " ",$input_file->getline; #masked
			
			my $background = "start";
			my $last_start = 1;
			my $last_pos = -1;
			for(my $i=0; $i<@basecaller;++$i)
			{
				#print @{$marker_container{1}}."\n";
				my @internal_line=split " ",@{$marker_container{$chr}}[$i];
				#print @internal_line."\n";
				#last;
				$file_writer->print($sampleID."\t".$chr."\t".$internal_line[1]."\t".$basecaller[$i]."\t".$internal_line[2]."\t".$internal_line[3]."\t".$internal_line[4]."\t".$internal_line[5]."\n");			
			}
		}
	}	
}



