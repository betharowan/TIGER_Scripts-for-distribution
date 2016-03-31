#! /usr/bin/perl
use strict;

my $usage = "$0 allele_count_refinement breaks_file\n";

my $ac_file = shift or die $usage;
my $bf_file = shift or die $usage;

my %P2ID = ();
my %GENO_POS = ();
my %GENO_P1 = ();
my %GENO_P1C = ();
my %GENO_P2 = ();
my %GENO_P2C = ();

my $MIN_ISLAND_SIZE = 400000;
#my $MIN_ISLAND_SIZE = 0;
#my $MIN_ISLAND_SIZE=0;

####################################################################
## Read in genotyping at an extended set of marker:

open FILE, $ac_file or die "cannot open file\n";
my $c = 0;
my $cc = -1;
while (<FILE>) {
	my @a = split " ";
	$c = 0 if ($cc != $a[0]);
	$P2ID{$a[0]}{$a[1]} = $c;
	push @{$GENO_POS{$a[0]}}, $a[1];
	push @{$GENO_P1{$a[0]}}, $a[2];
	push @{$GENO_P1C{$a[0]}}, $a[3];
	push @{$GENO_P2{$a[0]}}, $a[4];
	push @{$GENO_P2C{$a[0]}}, $a[5];
	$cc = $a[0];
	$c++;
}
close FILE;



####################################################################
# Read in breaks file 

open FILE, $bf_file or die "cannot open file\n";

my $SAMPLE = "";
my %BREAK_STARTS = ();
my %BREAK_ENDS = ();
my %BREAK_GENOTYPES = ();

while (<FILE>) {
        my @a = split " ";

	push @{$BREAK_STARTS{$a[1]}}, $a[2];
	push @{$BREAK_ENDS{$a[1]}}, $a[3];
	push @{$BREAK_GENOTYPES{$a[1]}}, $a[4];
}
close FILE;



####################################################################
# Parse breaks and remove short interspersed double 
# recombination (USE WITH CARE!)

my %BREAK_REF_STARTS = ();
my %BREAK_REF_ENDS = ();
my %BREAK_REF_GENOTYPES = ();

foreach my $chr (sort {$a <=> $b} keys %BREAK_STARTS) {
	if ((@{$BREAK_STARTS{$chr}}+0)>2) { # otherwise there is not double recombination possible
                $BREAK_REF_STARTS{$chr}[0] = $BREAK_STARTS{$chr}[0];
                $BREAK_REF_ENDS{$chr}[0] = $BREAK_ENDS{$chr}[0];
                $BREAK_REF_GENOTYPES{$chr}[0] = $BREAK_GENOTYPES{$chr}[0];
		
		my $c = 1;
		my $last_valid = 0;

		for (my $i = 1; $i < (@{$BREAK_STARTS{$chr}}+0)-1; $i++) 
		{
			if ($BREAK_GENOTYPES{$chr}[$last_valid] eq $BREAK_GENOTYPES{$chr}[$i+1] and # there needs to the same genotype on both sides of an island 
				$BREAK_ENDS{$chr}[$i] - $BREAK_STARTS{$chr}[$i] +1 < $MIN_ISLAND_SIZE) { # and the region needs to be short
				print "found $chr $BREAK_ENDS{$chr}[$i] \n";
				$BREAK_REF_ENDS{$chr}[$c-1] = $BREAK_ENDS{$chr}[$i+1];
				$last_valid = $i+1;
				$i++;
			}
			else {
				$BREAK_REF_STARTS{$chr}[$c] = $BREAK_STARTS{$chr}[$i];
	                	$BREAK_REF_ENDS{$chr}[$c] = $BREAK_ENDS{$chr}[$i];
                		$BREAK_REF_GENOTYPES{$chr}[$c] = $BREAK_GENOTYPES{$chr}[$i];	
				$last_valid = $i;
				$c++;
			}
		}


		if ($last_valid == (@{$BREAK_STARTS{$chr}}-2) ) {
			my $le = (@{$BREAK_STARTS{$chr}})-1;
	                $BREAK_REF_STARTS{$chr}[$c] = $BREAK_STARTS{$chr}[$le];
        	        $BREAK_REF_ENDS{$chr}[$c] = $BREAK_ENDS{$chr}[$le];
                	$BREAK_REF_GENOTYPES{$chr}[$c] = $BREAK_GENOTYPES{$chr}[$le];
		}

	}
	else {
		$BREAK_REF_STARTS{$chr} = $BREAK_STARTS{$chr};
		$BREAK_REF_ENDS{$chr} = $BREAK_ENDS{$chr};
		$BREAK_REF_GENOTYPES{$chr} = $BREAK_GENOTYPES{$chr};
	}


#print STDERR "$chr start:", join(",", @{$BREAK_REF_STARTS{$chr}}), "\n";
#print STDERR "$chr ends:", join(",", @{$BREAK_REF_ENDS{$chr}}), "\n";
#print STDERR "$chr geno:", join(",", @{$BREAK_REF_GENOTYPES{$chr}}), "\n";

}


####################################################################
# Now perform the actual refinement based on the filter recombination

my $fn = $bf_file;
if (substr($fn, length($fn)-11, 11) eq ".breaks.txt") {
        $fn = substr($fn, 0, length($fn)-11);
}
open OUT2, ">".$fn.".refined.recomb.txt";
open OUT3, ">".$fn.".recomb.txt";
open OUT, ">".$fn.".refined.breaks.txt";

foreach my $chr (sort {$a<=>$b} keys %BREAK_REF_STARTS) {

	print OUT "$SAMPLE\t$chr\t".$BREAK_REF_STARTS{$chr}[0]."\t";

	my $start = $BREAK_REF_ENDS{$chr}[0];
	my $geno_up = $BREAK_REF_GENOTYPES{$chr}[0];

	for (my $i = 1; $i < (@{$BREAK_REF_STARTS{$chr}}+0); $i++) { 
		my $refined_start = $start;
		my $end = $BREAK_REF_STARTS{$chr}[$i];
		my $refined_end = $end;
		my $geno_down = $BREAK_REF_GENOTYPES{$chr}[$i];

		my $p = $P2ID{$chr}{$start} + 1;
		my $pos = $GENO_POS{$chr}[$p];

#print STDERR "Check recomb: $geno_up $start  ---  $end $geno_down\n";

		if (is_homo($geno_up)) { # if upstream genotype is homozygous try to refine
			UPREFINE: while ($pos < $end) {
				if (new_allele_present($geno_up, $GENO_P1C{$chr}[$p], $GENO_P2C{$chr}[$p])) {
					last UPREFINE;
				}
				if (old_allele_present($geno_up, $GENO_P1C{$chr}[$p], $GENO_P2C{$chr}[$p])) {
					$refined_start = $pos;
				}

				## increase position
				$p++;
				$pos = ${$GENO_POS{$chr}}[$p];
			
			}
		}
		else {
			UPREFINE: while ($pos < $end) {
                                if (non_bg_allele_present($geno_down, $GENO_P1C{$chr}[$p], $GENO_P2C{$chr}[$p])) {
					$refined_start = $pos;
                                }

                                ## increase position
                                $p++;
                                $pos = ${$GENO_POS{$chr}}[$p];
                        }
		}


		$p = $P2ID{$chr}{$end} - 1;
                $pos = $GENO_POS{$chr}[$p];

		if (is_homo($geno_down)) { # if upstream genotype is homozygous try to refine

                        DOWNREFINE: while ($pos > $refined_start) {
                                if (new_allele_present($geno_up, $GENO_P1C{$chr}[$p], $GENO_P2C{$chr}[$p])) {
                                        last DOWNREFINE;
                                }
                                if (old_allele_present($geno_up, $GENO_P1C{$chr}[$p], $GENO_P2C{$chr}[$p])) {
                                        $refined_end = $pos;
                                }

                                ## decrease position
                                $p--;
                                $pos = ${$GENO_POS{$chr}}[$p];
                        }
                }
		else {
                        DOWNREFINE: while ($pos > $refined_start) {
                                if (non_bg_allele_present($geno_down, $GENO_P1C{$chr}[$p], $GENO_P2C{$chr}[$p])) {
                                        $refined_start = $pos;
                                }

                                ## decrease position
                                $p--;
                                $pos = ${$GENO_POS{$chr}}[$p];
                        }
                }

		print OUT2 "$SAMPLE\t$chr\t$refined_start\t$refined_end\t",($refined_end-$refined_start-1),"\t",($refined_start+int(($refined_end-$refined_start)/2)),"\n"; 
		print OUT3 "$SAMPLE\t$chr\t$start\t$end\t",($end-$start-1),"\t",($refined_start+int(($end-$start)/2)),"\n"; 

		print OUT "$refined_start\t$geno_up\n";
		print OUT "$SAMPLE\t$chr\t".$refined_end."\t";

		# reset for the next recombination:
	        $start = $BREAK_REF_ENDS{$chr}[$i];
	        $geno_up = $BREAK_REF_GENOTYPES{$chr}[$i];

	}
	
#print STDERR "finish genome   <$start>   <$geno_up>\n";
	print OUT "$start\t$geno_up\n";

}


##################################################################################################
##################################################################################################
##################################################################################################

sub is_homo {
	my ($geno) = @_;

	if ($geno eq "CC" or $geno eq "LL") {
		return 1;
	}
	return 0;
}


sub non_bg_allele_present {
        my ($geno, $ac1, $ac2) = @_;

        if ($geno eq "CC") {
                if ($ac2 != 0) {
                        return 1;
                }
                else {
                        return 0;
                }
        }
        elsif ($geno eq "LL") {
                if ($ac1 != 0) {
                        return 1;
                }
                else {
                        return 0;
                }
        }
        else {
                if ($ac1 != 0 or $ac2 != 0) {
                        return 1;
                }
                else {
                        return 0;
                }
        }
}


sub old_allele_present {
	my ($old_geno, $ac1, $ac2) = @_;
	
	if ($old_geno eq "CC") {
                if ($ac1 != 0) {
                        return 1;
                }
                else {
                        return 0;
                }
        }
	elsif ($old_geno eq "LL") {
                if ($ac2 != 0) {
                        return 1;
                }
                else {
                        return 0;
                }
        }
	else {
                if ($ac1 != 0 or $ac2 != 0) {
                        return 1;
                }
                else {
                        return 0;
                }
        }
}


sub new_allele_present {
	my ($old_geno, $ac1, $ac2) = @_;

	if ($old_geno eq "CC") {
		if ($ac2 != 0) {
			return 1;
		}
		else {
			return 0;
		}
	}
	elsif ($old_geno eq "LL") {
		if ($ac1 != 0) {
                        return 1;
                }
                else {
                        return 0;
                }
	}
	else {
                return 0;
	}

}

