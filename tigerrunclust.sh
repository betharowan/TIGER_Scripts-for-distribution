#!/bin/bash
#tigerrun.sh
#runs TIGER pipeline on cluster
#23 September 2015
#By Beth Rowan

#specify path to Java installation and java memory
export JAVA_HOME=/ebio/abt6_projects9/natural_variation_recombination/code/jre1.7.0_80/bin/
export _JAVA_OPTIONS=-Xmx1024M

#specify variables
TIGER=$1
DATA=$2
NUM=$SGE_TASK_ID

#run base caller from TIGER directory
cd $TIGER
java -jar base_caller.jar -r $DATA/INDEX$NUM/input_corrected$NUM.txt -o $DATA/INDEX$NUM/allele_count_base_call$NUM.txt -n bi

#run allele frequency estimator TIGER directory
java -jar allele_freq_estimator.jar -r $DATA/INDEX$NUM/input_corrected$NUM.txt -o $DATA/INDEX$NUM/frequencies_for_bmm$NUM.txt -n bi -w 1000

#apply the beta mixture model
R --slave --vanilla --args $DATA/INDEX$NUM/frequencies_for_bmm$NUM.txt $DATA/INDEX$NUM/bmm.intersections$NUM.txt < beta_mixture_model.R

#Prepare files for HMM probability estimation using the BASECALLER output
perl prep_prob.pl -s $NUM -m $DATA/INDEX$NUM/input_corrected$NUM.txt -b $DATA/INDEX$NUM/allele_count_base_call$NUM.txt -c $TIGER/TAIR10_chrSize.txt -o $DATA/INDEX$NUM/file_for_probabilities$NUM.txt

#Calculate transmission and emission probabilities for the HMM and breakpoint estimates based on the sliding window
cd $DATA/INDEX$NUM

perl $TIGER/hmm_prob.pl -s $DATA/INDEX$NUM/frequencies_for_bmm$NUM.txt -p $DATA/INDEX$NUM/file_for_probabilities$NUM.txt -o INDEX$NUM -a $DATA/INDEX$NUM/bmm.intersections$NUM.txt -c $TIGER/TAIR10_chrSize.txt

#Run the HMM
cd $TIGER
java -jar hmm_play.jar -r $DATA/INDEX$NUM/allele_count_base_call$NUM.txt -o $DATA/INDEX$NUM/hmm.out$NUM.txt -t bi -z $DATA/INDEX$NUM/INDEX$NUM\_hmm_model 

#Get rough estimate of recombination breakpoint positions
perl prepare_break.pl -s $NUM -m $DATA/INDEX$NUM/input_corrected$NUM.txt -b $DATA/INDEX$NUM/hmm.out$NUM.txt -c $TIGER/TAIR10_chrSize.txt -o $DATA/INDEX$NUM/rough_COs$NUM.txt

#refine breaks by using nearby markers that had been previously filtered out
perl $TIGER/refine_recombination_break.pl $DATA/INDEX$NUM/input_complete$NUM.txt $DATA/INDEX$NUM/rough_COs$NUM.breaks.txt

#smooth out breaks
perl $TIGER/breaks_smoother.pl -b $DATA/INDEX$NUM/rough_COs$NUM.refined.breaks.txt -o $DATA/INDEX$NUM/$NUM\refined.corrected.breaks

#visualize output
cd $DATA/INDEX$NUM
R --slave --vanilla --args $NUM sample_$NUM.visual_out.pdf $DATA/INDEX$NUM/rough_COs$NUM.breaks.txt $DATA/INDEX$NUM/rough_COs$NUM.refined.breaks.txt $DATA/INDEX$NUM/$NUM\refined.corrected.breaks $DATA/INDEX$NUM/frequencies_for_bmm$NUM.txt $DATA/INDEX$NUM/INDEX$NUM\_sliding_window.breaks.txt < $TIGER/plot_genotyping.R
