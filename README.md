# TIGER_Scripts-for-distribution

This repository contains the scripts and documentation for the TIGER pipeline described in Rowan, Patel et al. 2015 in G3 (doi: 10.1534/g3.114.016501) for reconstructing recombination genomes from low-coverage, short-read sequencing data.

#TIGER Documentation

The TIGER Scripts were written by Vipul Patel and Korbinian Schneeberger.

Scripts for generating TIGER input files were written by Beth Rowan.

The get_subset.pl script was written by Joerg Hagmann.

This documentation file was written by Beth Rowan.

Questions about these scripts can be directed to Beth Rowan (beth.rowan@tuebingen.mpg.de)

Last updated 15 January 2016

#Getting started

1.Create input files and get set up

Convert SHORE alignment output (map.list) or BWA alignment file (.bam) into an input file. The input file has a format like this:
<#chr><#pos><#parent1allele><#reads supportingparent1allele><parent2allele><#readssupportingparent2allele>

Make one file with all markers with only indels filtered out (this one is called “complete”)

In our example, the file is called input.complete24.txt

Make a second file with all filters applied (this one is called “corrected”)

In our example, the file is called input.corrected24.txt 

Put these input files into a separate directory (called input).

You will also need a tab-delimited file showing the chromosome numbers and their lengths.

We have scripts for generating these input files from either sam (samtotiger.sh) or bam (bamtotiger.sh) files. You will need to provide your own filtered and unfiltered marker files. All that is required for the marker files is that the chromosome number is in the first column and the position (in bp) is in the second column.

Finally, you will need to install the Java 7 environment to run the java scripts.

What follows is a list of the individual commands for the TIGER pipeline. I have also provided an example of how to run all of the steps in parallel for our compute cluster (tigerclusterrun.sh).

2. Run the base caller on the corrected input file.

General command:

java -jar base_caller.jar -r $CORRECTED_INPUT_FILE -o $BASE_CALL_OUTPUT -n bi

Notes:
-r specifies the corrected marker counts input file
-o designates the output file
-n specifies that it is biparental

In the provided example, this was the command I ran:
java -jar base_caller.jar -r ./input/allele_count.corrected.txt -o allele_count.base_call.txt -n bi

This output file is basically a single line for each chromosome with all of the positions assigned a base caller genotype. Note that the genotype notation uses C and L for the two different parental alleles instead of A and B as described in Rowan, Patel et al. 2015.

3. Run the allele frequency estimator

General command:

java -jar allele_freq_estimator.jar -r $CORRECTED_INPUT_FILE -o $ALLELE_FREQUENCIES_FOR_BMM -n bi -w 1000

-r specifies the input file (this is again the corrected marker counts input file)
-0 specifies an output file
-n specifies that it’s biparental
-w specifies window size

In the provided example, this was the command I ran:
java -jar allele_freq_estimator.jar -r ./input/allele_count.corrected.txt -o frequencies_for_bmm.txt -n bi -w 1000

The output file has this format:
<chr><pos><number specifying read ratio distribution>

4. Apply the beta mixture model

General command:

R --slave --vanilla --args $ALLELE_FREQUENCIES_FOR_BMM $BMM_OUTPUT < beta_mixture_model.R

First argument: This is the output of the allele_freq_estimator.jar command
Second argument: Specifies an output file.

The output file will contain two numbers that specify the intersections of the curves

In the example, this was the command that I ran:

R --slave --vanilla --args frequencies_for_bmm.txt bmm.intersections.txt < beta_mixture_model.R

Please note that this step will probably take several hours to process. It takes about three hours for me.

5.  Prepare files for HMM probability estimation using the BASECALLER output and the output of the beta mixture model.

General command:

perl prep_prob.pl -s $LABEL -m $CORRECTED_INPUT_FILE -b $BASE_CALL_OUTPUT -c $CHRSIZES -o $FILE_FOR_PROB

-s specifies sample label (24 in our example case)
-m corrected marker input file
-b file that was the output from the base caller
-c chromosome sizes file
-o specifies output file

The output file looks like this:
<sample><chr><pos><basecaller><parent1><reads for parent1><parent2><reads for parent2>

In the example, here is the command that I ran:
perl prep_prob.pl -s 24 -m ./input/allele_count.corrected.txt -b allele_count.base_call.txt -c TAIR10_chrSize.txt -o file_for_probabilities.txt

6.  Calculate transmission and emission probabilities for the HMM
General command:
perl hmm_prob.pl -s $ALLELE_FREQUENCIES_FOR_BMM -p $FILE_FOR_PROB -o sample -a $BMM_OUTPUT -c $CHRSIZES

-s output from the sliding window (allele frequency estimator)
-p output file from previous script (prep_prob.pl)
-o sample (gives a prefix for the output files)
-a output from beta mixture model
-c chromosome sizes file

This gives two output files:
sample_hmm_model (probabilities for the HMM)
sample_sliding_window (genotyping only based on the sliding window)

In the provided example, this is the command I used:

perl hmm_prob.pl -s frequencies_for_bmm.txt -p file_for_probabilities.txt -o sample -a bmm.intersections.txt -c TAIR10_chrSize.txt

7.  Run the HMM

General command:

java -jar hmm_play.jar -r $BASE_CALL_OUTPUT -o $HMM_OUTPUT -t bi -z sample_hmm_model

-r output of the base caller file
-o specify output file
-t bi (for biparental)
-z output from last script (probabilities contained in file sample_hmm_model)

Example script:
java -jar hmm_play.jar -r allele_count.base_call.txt -o hmm.out.txt -t bi -z sample_hmm_model

The output file has a single line of base caller genotypes for each chromosome, then two lines of genotypes inferred from HMM, and a fourth line with just the probabilities.

8. Get rough estimate of recombination breakpoint positions

General command:

perl prepare_break.pl -s $LABEL -m $CORRECTED_INPUT_FILE -b hmm.out.txt -c $CHRSIZES -o $ROUGH_CO


-s sample label
-m marker file with corrected markers
-b output file from previous script (hmm_play.jar)
-c chromosome sizes file
-o output file

Example script:
perl prepare_break.pl -s 24 -m input/allele_count.corrected.txt  -b hmm.out.txt -c TAIR10_chrSize.txt -o rough_COs.txt

Two output files:
$ROUGH_CO.txt
<sample label><chr><pos><basecaller genotype><hmm_inferred_genotype1><hmm_inferred_genotype2><parent1alelle><countsfor parent1><parent2allele><countsforparent2>

$ROUGH_CO.breaks.txt
<sample number><chr><startpos><endpos><genotype>

9. Refine recombination breaks

General command:
perl refine_recombination_break.pl $COMPLETE_INPUT_FILE $ROUGH_CO.breaks.txt

First argument: marker input file with complete data
Second argument: “breaks” output file from previous script

Example script:
perl refine_recombination_break.pl input/allele_count.complete.txt rough_COs.breaks.txt

Gives output
$ROUGH_CO.recomb.txt
$ROUGH_CO.refined.breaks.txt
$ROUGH_CO.refined.recomb.txt

10. Smooth out breaks

General command:

perl breaks_smoother.pl -b $ROUGH_CO.refined.breaks.txt -o $SMOOTH_CO

-b “refined breaks” output from previous script
-o specify output file

Example script:
perl breaks_smoother.pl -b rough_COs.refined.breaks.txt -o corrected.refined.breaks.txt

output
Like the “breaks” output, but with corrected breakpoints based on the markers that had been filtered out.


11. Visualize output

General output

R --slave --vanilla --args  plot_genotyping.R
sample_id $VISUALIZATION.pdf $ROUGH_CO.breaks.txt $ROUGH_CO.refined.breaks.txt corrected.refined.breaks.txt $ALLELE_FREQUENCIES_FOR_BMM sample_sliding_window.breaks.txt < plot_genotyping.R

usage: <sample label><output file><unrefined breaks file><refined breaks file><corrected refined breaks file><allele frequency estimator output><file with breaks based only on sliding window output>

Example script:
R --slave --vanilla --args 24 visual_out.pdf rough_COs.breaks.txt rough_COs.refined.breaks.txt corrected.refined.breaks.txt frequencies_for_bmm.txt sample_sliding_window.breaks.txt < plot_genotyping.R
