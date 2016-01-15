use strict;
use warnings;
use Class::Struct;
use Getopt::Std;
use IO::File;
use Time::HiRes;
use Cwd;
my % options=();
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
use IO::Handle;

my$opt_string='hs:p:o:a:c:';

getopts("$opt_string",\%options) or usage();
usage()if $options{h};
#change here if you have a gzip file
my $GZIP=1;

my $REPORT_ALL = 0;

my $sliding_window_file="";
my $hmm_converted_file="";
my $out_file_handler = IO::File->new();
my $p1_min_val=-50;
my $p2_min_val=50;
my $output_file_name="output.txt";
my %POS=();
my %genotypes=();
my %transistion=();
my @starting_transition=(0,0,0);
my $chr_sizes_file="";
my %chr_sizes=();
my $allel_border_file="";
init();

main();

sub usage
{
	print "perl $0 -s sliding_window_outputFile -o output -p hmm_converted_outputFile -a allele_borders -c chr_sizes\n";
}

sub init
{
	if($options{s} && $options{p} && $options{c} && $options{p})
	{
		$sliding_window_file=$options{s};
		$hmm_converted_file=$options{p};
		$chr_sizes_file=$options{c};
		$allel_border_file=$options{a}
		
	}
	else
	{
		die usage();
	}
	if($options{o})
	{
		$output_file_name=$options{o};
	}
}


sub main
{
	if(-e $allel_border_file)
	{
		my $input_file_1 =new IO::File($allel_border_file, "r") or die "could not open $allel_border_file: $!\n";
		my $line = $input_file_1->getline;	
	
		chomp($line);
	
		$p1_min_val = $line;
		$line = $input_file_1->getline;
		$p2_min_val=$line;	
	}
	else
	{
		system("R --slave --vanilla --args $sliding_window_file 1 < beta_mixture_model.R"."\n");
		
		if(-e $allel_border_file)
		{
			my $input_file_1 =new IO::File($allel_border_file, "r") or die "could not open $allel_border_file: $!\n";
			my $line = $input_file_1->getline;	
		
			chomp($line);
		
			$p1_min_val = $line;
			$line = $input_file_1->getline;
			$p2_min_val=$line;	
		}
		else
		{	
			#default.. if somethings goes wrong
			$p1_min_val=-50;
			$p2_min_val=50;
		}
	}	
	
	my $input_file_1 =new IO::File($chr_sizes_file, "r") or die "could not open $chr_sizes_file: $!\n";
	
	while(my $line = $input_file_1->getline)
	{
		my @a=split " ",$line;
		my $chr = $a[0];
		$chr=~ s/\D//g;
		$chr_sizes{$chr}=$a[1];
	}
	
	
	$input_file_1 = new IO::File($sliding_window_file, "r") or die "could not open $sliding_window_file: $!\n";
	my $input_file_2 = new IO::File($hmm_converted_file, "r") or die "could not open $hmm_converted_file: $!\n";
	
	my $file_writer=IO::File->new();
#	$file_writer->open(">".$output_file_name."_sliding_window") or die "Could not open ".$output_file_name."_sliding_window $!\n" ;
	my $file_writer2=IO::File->new();
	$file_writer2->open(">".$output_file_name."_sliding_window.breaks.txt") or die "Could not open ".$output_file_name."_sliding_window $!\n" ;
	my $prev="";
	
	my $chr=0;
	my $start_pos=1;
	my $sample_id="";

	initTransition("CC");
	initTransition("LL");
	initTransition("CL");
	while (my $line_file_1 = $input_file_1->getline) 
	{
		my $line_file_2=$input_file_2->getline;
		my @a_sliding_window_file = split " ",$line_file_1;
		my @a_hmm_converted_file = split " ",$line_file_2;
		#100     1       154     CU      CL      CL      C       3       T       0
		#1       154     -31.643356643356647
		
		if($chr != $a_hmm_converted_file[1] || $sample_id ne $a_hmm_converted_file[0])
		{
			if($chr!=0)
			{
				$file_writer2->print($sample_id."\t".$chr."\t".$start_pos."\t".$chr_sizes{$chr}."\t".$prev."\n");
			}
			$sample_id=$a_hmm_converted_file[0];
			$prev="";
			$chr=$a_hmm_converted_file[1];
		}
		if($prev eq "")
		{
			$start_pos=1;
		}
		my $prev_temp="";
		$prev_temp=$prev;		
		
		if($a_sliding_window_file[2]<=$p1_min_val)
		{
			add_genotype("LL",$a_hmm_converted_file[3]);
			addTransition(\$prev,"LL");
			#$file_writer->print($a_hmm_converted_file[0]."\t".$a_hmm_converted_file[1]."\t".$a_hmm_converted_file[2]."\t".$a_hmm_converted_file[3]."\t".$a_hmm_converted_file[4]."\tLL"."\t".$a_hmm_converted_file[6]."\t".$a_hmm_converted_file[7]."\t".$a_hmm_converted_file[8]."\t".$a_hmm_converted_file[9]."\n")
		}
		elsif($a_sliding_window_file[2]>=$p2_min_val)
		{
			add_genotype("CC",$a_hmm_converted_file[3]);
			addTransition(\$prev,"CC");
			#$file_writer->print($a_hmm_converted_file[0]."\t".$a_hmm_converted_file[1]."\t".$a_hmm_converted_file[2]."\t".$a_hmm_converted_file[3]."\t".$a_hmm_converted_file[4]."\tCC"."\t".$a_hmm_converted_file[6]."\t".$a_hmm_converted_file[7]."\t".$a_hmm_converted_file[8]."\t".$a_hmm_converted_file[9]."\n")
		}
		else
		{
			add_genotype("CL",$a_hmm_converted_file[3]);
			addTransition(\$prev,"CL");
			#$file_writer->print($a_hmm_converted_file[0]."\t".$a_hmm_converted_file[1]."\t".$a_hmm_converted_file[2]."\t".$a_hmm_converted_file[3]."\t".$a_hmm_converted_file[4]."\tCL"."\t".$a_hmm_converted_file[6]."\t".$a_hmm_converted_file[7]."\t".$a_hmm_converted_file[8]."\t".$a_hmm_converted_file[9]."\n")
		}
		if($prev ne $prev_temp && $prev_temp ne "")
		{
			$file_writer2->print($sample_id."\t".$chr."\t".$start_pos."\t".$a_hmm_converted_file[2]."\t".$prev_temp."\n");
			$start_pos=$a_hmm_converted_file[2];
		}				
	}
	$file_writer2->print($sample_id."\t".$chr."\t".$start_pos."\t".$chr_sizes{$chr}."\t".$prev."\n");
	$file_writer2->flush();
	$file_writer2->close();
	normalise();
	#$file_writer->flush();
	#$file_writer->close();
	$file_writer->open(">".$output_file_name."_hmm_model") or die "Could not open ".$output_file_name."_hmm_model $!\n" ;
	############################################
	#
	#write f2 starting prob 25 25 50
	#
	$file_writer->print("starting probs:\n");
	if($starting_transition[0]!=0)
	{
		$file_writer->print("CC ". ($starting_transition[0]/($starting_transition[0]+$starting_transition[1]+$starting_transition[2]))."\n");
	}
	else
	{
		$file_writer->print("CC 0.25\n");
	}
	if($starting_transition[1]!=0)
	{
		$file_writer->print("LL ". ($starting_transition[1]/($starting_transition[0]+$starting_transition[1]+$starting_transition[2]))."\n");
	}
	else
	{
		$file_writer->print("LL 0.25\n");
	}
	if($starting_transition[2]!=0)
	{
		$file_writer->print("CL ". ($starting_transition[2]/($starting_transition[0]+$starting_transition[1]+$starting_transition[2]))."\n");
	}
	else
	{
		$file_writer->print("CL 0.50\n");
	}
	
	$file_writer->print("transition:\n");
	if(!exists $transistion{"CC"})
	{
		$file_writer->print("CC -> CC "."0"."\n");
		$file_writer->print("CC -> LL "."0"."\n");
		$file_writer->print("CC -> CL "."0"."\n");
	}
	else
	{
		$file_writer->print("CC -> CC ".@{$transistion{"CC"}}[0]."\n");
		$file_writer->print("CC -> LL ".@{$transistion{"CC"}}[1]."\n");
		$file_writer->print("CC -> CL ".@{$transistion{"CC"}}[2]."\n");
	}
	$file_writer->print("\n");
	if(!exists $transistion{"LL"})
	{
		$file_writer->print("LL -> CC "."0"."\n");
		$file_writer->print("LL -> LL "."0"."\n");
		$file_writer->print("LL -> CL "."0"."\n");
	}
	else
	{
		$file_writer->print("LL -> CC ".@{$transistion{"LL"}}[0]."\n");
		$file_writer->print("LL -> LL ".@{$transistion{"LL"}}[1]."\n");
		$file_writer->print("LL -> CL ".@{$transistion{"LL"}}[2]."\n");	
	}
	
	$file_writer->print("\n");
	#print length(@{$transistion{"CL"}})."\n";
	if(!exists $transistion{"CL"})
	{
		$file_writer->print("CL -> CC "."0"."\n");
		$file_writer->print("CL -> LL "."0"."\n");
		$file_writer->print("CL -> CL "."0"."\n");
	}
	else
	{
		$file_writer->print("CL -> CC ".@{$transistion{"CL"}}[0]."\n");
		$file_writer->print("CL -> LL ".@{$transistion{"CL"}}[1]."\n");
		$file_writer->print("CL -> CL ".@{$transistion{"CL"}}[2]."\n");
	}
	$file_writer->print("\n");
	$file_writer->print("emission:\n");
	if(!exists $genotypes{"CC"})
	{
		$file_writer->print("CC -> Discrete distribution --- CC "."0".", LL "."0".", CL "."0".", UN "."0".", CU "."0".", LU "."0".", AA 0, AB 0, AC 0, AL 0, BB 0, BC 0, BL 0, AP 0, BP 0, CP 0, LP 0, OP 0, OC 0, OL 0, OK 0, MP 0, MB 0, ML 0, MJ 0, NP 0, NB 0, NC 0, NI 0, IP 0, IA 0, IL 0, JP 0, JA 0, JC 0, KP 0, KA 0, KB 0, QP 0, QL 0, RP 0, RC 0, TP 0, TA 0, SP 0, SB 0\n");
	
	}
	else
	{
		$file_writer->print("CC -> Discrete distribution --- CC ".@{$genotypes{"CC"}}[0].", LL ".@{$genotypes{"CC"}}[1].", CL ".@{$genotypes{"CC"}}[2].", UN ".@{$genotypes{"CC"}}[5].", CU ".@{$genotypes{"CC"}}[3].", LU ".@{$genotypes{"CC"}}[4].", AA 0, AB 0, AC 0, AL 0, BB 0, BC 0, BL 0, AP 0, BP 0, CP 0, LP 0, OP 0, OC 0, OL 0, OK 0, MP 0, MB 0, ML 0, MJ 0, NP 0, NB 0, NC 0, NI 0, IP 0, IA 0, IL 0, JP 0, JA 0, JC 0, KP 0, KA 0, KB 0, QP 0, QL 0, RP 0, RC 0, TP 0, TA 0, SP 0, SB 0\n");
	}
	
	if(!exists $genotypes{"LL"})
	{

		$file_writer->print("LL -> Discrete distribution --- CC "."0".", LL "."0".", CL "."0".", UN "."0".", CU "."0".", LU "."0".", AA 0, AB 0, AC 0, AL 0, BB 0, BC 0, BL 0, AP 0, BP 0, CP 0, LP 0, OP 0, OC 0, OL 0, OK 0, MP 0, MB 0, ML 0, MJ 0, NP 0, NB 0, NC 0, NI 0, IP 0, IA 0, IL 0, JP 0, JA 0, JC 0, KP 0, KA 0, KB 0, QP 0, QL 0, RP 0, RC 0, TP 0, TA 0, SP 0, SB 0\n");
	
	}
	else
	{
		$file_writer->print("LL -> Discrete distribution --- CC ".@{$genotypes{"LL"}}[0].", LL ".@{$genotypes{"LL"}}[1].", CL ".@{$genotypes{"LL"}}[2].", UN ".@{$genotypes{"LL"}}[5].", CU ".@{$genotypes{"LL"}}[3].", LU ".@{$genotypes{"LL"}}[4].", AA 0, AB 0, AC 0, AL 0, BB 0, BC 0, BL 0, AP 0, BP 0, CP 0, LP 0, OP 0, OC 0, OL 0, OK 0, MP 0, MB 0, ML 0, MJ 0, NP 0, NB 0, NC 0, NI 0, IP 0, IA 0, IL 0, JP 0, JA 0, JC 0, KP 0, KA 0, KB 0, QP 0, QL 0, RP 0, RC 0, TP 0, TA 0, SP 0, SB 0\n");
	}
	
	if(!exists $genotypes{"CL"})
	{

		$file_writer->print("CL -> Discrete distribution --- CC "."0".", LL "."0".", CL "."0".", UN "."0".", CU "."0".", LU "."0".", AA 0, AB 0, AC 0, AL 0, BB 0, BC 0, BL 0, AP 0, BP 0, CP 0, LP 0, OP 0, OC 0, OL 0, OK 0, MP 0, MB 0, ML 0, MJ 0, NP 0, NB 0, NC 0, NI 0, IP 0, IA 0, IL 0, JP 0, JA 0, JC 0, KP 0, KA 0, KB 0, QP 0, QL 0, RP 0, RC 0, TP 0, TA 0, SP 0, SB 0\n");
	}
	else
	{
		$file_writer->print("CL -> Discrete distribution --- CC ".@{$genotypes{"CL"}}[0].", LL ".@{$genotypes{"CL"}}[1].", CL ".@{$genotypes{"CL"}}[2].", UN ".@{$genotypes{"CL"}}[5].", CU ".@{$genotypes{"CL"}}[3].", LU ".@{$genotypes{"CL"}}[4].", AA 0, AB 0, AC 0, AL 0, BB 0, BC 0, BL 0, AP 0, BP 0, CP 0, LP 0, OP 0, OC 0, OL 0, OK 0, MP 0, MB 0, ML 0, MJ 0, NP 0, NB 0, NC 0, NI 0, IP 0, IA 0, IL 0, JP 0, JA 0, JC 0, KP 0, KA 0, KB 0, QP 0, QL 0, RP 0, RC 0, TP 0, TA 0, SP 0, SB 0");
	}
		
	$file_writer->flush();
	$file_writer->close();
}

sub initTransition
{
	my($prev_state_ref)=@_;
	if(!exists $transistion{$prev_state_ref})
	{
		$transistion{$prev_state_ref}=[0,0,0];
	}
}


sub addTransition
{
	my($prev_state_ref,$current_state)=@_;
	
	if($$prev_state_ref eq "")
	{
		if($current_state eq "CC")
		{
			$starting_transition[0]+=1;
		}
		elsif($current_state eq "LL")
		{
			$starting_transition[1]+=1;
		}
		else
		{
			$starting_transition[2]+=1;
		}
		$$prev_state_ref=$current_state;
	}
	else
	{
		if(!exists $transistion{$$prev_state_ref})
		{
			$transistion{$$prev_state_ref}=[0,0,0];
		}
		my @temp=@{$transistion{$$prev_state_ref}};
		
		if($current_state eq "CC")
		{
			$temp[0]+=1;
			$transistion{$$prev_state_ref}=[@temp];
		}
		elsif($current_state eq "LL")
		{
			$temp[1]+=1;
			$transistion{$$prev_state_ref}=[@temp];
		}
		else
		{
			$temp[2]+=1;
			$transistion{$$prev_state_ref}=[@temp];
		}
		$$prev_state_ref=$current_state;
	}
}


sub normalise
{
	my @hash_values = keys %genotypes;
	
	foreach my $val (@hash_values)
	{
		my $sum=0;
		for(my $i=0;$i<@{$genotypes{$val}};++$i)
		{
			$sum+=@{$genotypes{$val}}[$i];
		}
		
		for(my $i=0;$i<@{$genotypes{$val}};++$i)
		{
			if(@{$genotypes{$val}}[$i]>0)
			{
				@{$genotypes{$val}}[$i] =@{$genotypes{$val}}[$i]/$sum;
			} 
		}		
	}
	
	@hash_values = keys %transistion;
	
	foreach my $val (@hash_values)
	{
		my $sum=0;
		for(my $i=0;$i<@{$transistion{$val}};++$i)
		{
			$sum+=@{$transistion{$val}}[$i];
		}
		
		for(my $i=0;$i<@{$transistion{$val}};++$i)
		{
			if(@{$transistion{$val}}[$i]>0)
			{
				@{$transistion{$val}}[$i] =@{$transistion{$val}}[$i]/$sum;
			}
		}
	}	
}


sub add_genotype
{
	my($genotype_sliding_window,$genotype_hmm)=@_;
	if(!exists $genotypes{$genotype_sliding_window})
	{
		$genotypes{$genotype_sliding_window}=[0,0,0,0,0,0];
	}
	
	my @temp=@{$genotypes{$genotype_sliding_window}};
	if($genotype_hmm eq "CC")
	{
		$temp[0]+=1;
		$genotypes{$genotype_sliding_window}=[@temp];
	}
	elsif($genotype_hmm eq "LL")
	{
		$temp[1]+=1;
		$genotypes{$genotype_sliding_window}=[@temp];
		
	}
	elsif($genotype_hmm eq "CL")
	{
		$temp[2]+=1;
		$genotypes{$genotype_sliding_window}=[@temp];
		
	}
	elsif($genotype_hmm eq "CU")
	{
		$temp[3]+=1;
		$genotypes{$genotype_sliding_window}=[@temp];
		
	}
	elsif($genotype_hmm eq "LU")
	{
		$temp[4]+=1;
		$genotypes{$genotype_sliding_window}=[@temp];
		
	}
	elsif($genotype_hmm eq "UN")
	{
		$temp[5]+=1;
		$genotypes{$genotype_sliding_window}=[@temp];
		
	}
}

