#get_subset.pl
# by Joerg Hagmann

#! /usr/bin/perl
use strict;

my $usage= "\n$0  <complete_set_file>  <list_of_columns(1-based)>  <subset_file>  <list_of_subset_columns(1-based)> <0 for intersection, 1 for completeset part, 2 for subset part, 3 for intersection with both lines printed in one, 4 for union>\n\n" ;

my $file = shift or die $usage;
my $cols_f = shift or die $usage;
my $subset = shift or die $usage;
my $cols_s = shift or die $usage;
my $rev = shift;

$rev = 0 if (!defined $rev);

my @cols_f = split/,/,$cols_f;
my @cols_s = split/,/,$cols_s;

my %POS = ();

open S, $subset or die $usage;
while (my $line = <S>) {
	chomp $line;
        my @a = split" ", $line;
	my $id = "";
	foreach (@cols_s) { $id .= $a[$_-1]."#"; }
	$POS{$id} = $line;
}
close S;

open FILE, $file or die $usage;
my $c=0;
while (my $line = <FILE>) {
	chomp $line;
        my @a = split" ", $line;
	my $id="";
	foreach (@cols_f) { $id .= $a[$_-1]."#"; }
	if (exists($POS{$id})) {
		if ($rev == 0 || $rev == 3 || $rev == 4) {
			#chomp $line if ($rev==3);
			print $line;
			print "\n" if ($rev!=3);
			print "\t".$POS{$id}."\n" if ($rev==3);
			$c++;
			delete $POS{$id} if ($rev == 4);
			#last if ($c>=scalar(keys %POS));
		}
		elsif ($rev == 2) {
			delete $POS{$id};
		}
	}
	elsif ($rev == 1 || $rev == 4) {
		print $line."\n";
	}
}
close FILE;

if ($rev == 2 || $rev == 4) {
	foreach (keys %POS) {
		print $POS{$_}."\n";
	}
}

