#!/bin/bash
#samtotiger.sh
#uses sam files to generate vcf files and files for TIGER input
#22 September 2015
#By Beth Rowan

#modified 13 April 2016
#By Vipul Patel

#usage
if [ $# -eq 0 ]
  then
echo "bamtotiger.sh <path to sam tools> <reference genome> <path to bcf tools> <path to tiger scripts> <complete_marker_file> <corrected_marker_file> <sampleno>"
exit 1
fi

#variables
SAMDIR=$1
REF=$2
BCFDIR=$3
TIGER=$4
MARKCOMP=$5
MARKCORR=$6
NUM=$7

#Convert sam file to sorted bam file
#from each sample directory
#cd into each directory
for i in $(eval echo {1..$NUM\});
do
    cd INDEX$i

         #convert to bam
#$SAMDIR/samtools view -bS INDEX$i.fastq.sam > INDEX$i.bam

#sort
$SAMDIR/samtools sort INDEX$i.bam INDEX$i.sorted

#convert sorted bam file to vcf
#$SAMDIR/samtools mpileup -vu -f $REF INDEX$i.sorted.bam | $BCFDIR/bcftools call -m -Ov >  INDEX$i.vcf
$SAMDIR/samtools mpileup -E --skip-indels --BCF --output-tags DP,DV,SP,DP4 -f $REF INDEX$i.sorted.bam | $BCFDIR/bcftools call --output-type v --multiallelic-caller --skip-variants indels -R $MARKCOMP  > INDEX$i.vcf
#extract DP4 lines
$BCFDIR/bcftools query -f '%CHROM %POS %REF %ALT %QUAL [ %INDEL %DP %DP4]\n' INDEX$i.vcf -o INDEX$i.comma.txt

#replace comma separators with tabs
tr ',' '\t' < INDEX$i.comma.txt > INDEX$i.tabbed.txt

#get rid of mito and cp reads; create columns for read counts for ref allele and alt allele by adding the first two of the DP4 fields and the second two of the DP4 fields. ##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">

awk '{if ($1 !="ChrC" && $1!="ChrM") print $1 "\t" $2 "\t" $3 "\t" $8+$9 "\t"  $4 "\t"  $10+$11}' INDEX$i.tabbed.txt > INDEX$i.input.temp

#call bash script to change chromosome column to numbers only (to match TIGER input)
bash $TIGER/chr_convert.sh INDEX$i.input.temp

#filter indels only (to generate "complete" file for TIGER)
#perl $TIGER/get_subset.pl input.out 1,2 $MARKCOMP 1,2 0 > input_complete$i.txt
perl $TIGER/get_subset_2.pl -c input.out -p $MARKCOMP -o input_complete$i.txt
#use stricter filters (to generate "corrected" file for TIGER)
#perl $TIGER/get_subset.pl input.out 1,2 $MARKCORR 1,2 0 > input_corrected$i.txt

perl $TIGER/get_subset_2.pl -c input.out -p $MARKCORR -o input_complete$i.txt


rm *tabbed*
rm *temp
rm *comma*

cd ../

done
