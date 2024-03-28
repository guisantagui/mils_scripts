#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=NEAT_enrich
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH -c 1
#SBATCH --time=00-24:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/cureMILS/proteomics/scripts/output_NEAT_enrich.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/cureMILS/proteomics/scripts/error_NEAT_enrich.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal


module load lang/R/4.0.5-foss-2020b

# Variables for the pipeline
########################################################################################################################

DE_file="/home/users/gsantamaria/projects/cureMILS/proteomics/results/DE/proteom_log2_sick_UT_vs_d_almSign_oplsda_res_sign.csv"
outNameNEAT="proteom_log2_sick_UT_vs_d_neat.csv"

Rscript NEAT_enrich.R --DE $DE_file --outName $outNameNEAT
