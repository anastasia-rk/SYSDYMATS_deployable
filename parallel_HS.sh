#!/bin/bash

#$ -P acsehpc
#$ -q acsehpc.q
#$ -pe smp 12
#$ -l rmem=4G
#$ -l h_rt=150:00:00
#$ -cwd 
#$ -m bea                            
#$ -M a.kadochnikova@sheffield.ac.uk
#$ -o errors_HS.txt
#$ -j y

module load apps/matlab/2019b
matlab -nodesktop -nojvm -nosplash -r parallel_main_HS