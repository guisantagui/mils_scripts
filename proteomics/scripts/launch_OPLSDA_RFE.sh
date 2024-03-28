#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=OPLS_RFE
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH -c 1
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/cureMILS/proteomics/scripts/output_OPLSDA_RFE.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/cureMILS/proteomics/scripts/error_OPLSDA_RFE.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal


module load lang/R/4.0.5-foss-2020b

# Variables for the pipeline
########################################################################################################################

input="/home/users/gsantamaria/projects/cureMILS/proteomics/results/preprocessing/proteom_log2.csv"
sampInfo="/home/users/gsantamaria/projects/cureMILS/proteomics/results/sample_info.csv"
respVar='treatment'
outDir="/home/users/gsantamaria/projects/cureMILS/proteomics/results/sup_analysis/OPLSDA_RFE/"
#splitTrainTest='0.66'
splitTrainTest='no'
filtSampsByDisPhen='sick'

Rscript OPLSDA_RFE.R $input --sampMetadata $sampInfo --respVar $respVar --outDir $outDir --splitTrainTest $splitTrainTest --filtSampsByDisPhen $filtSampsByDisPhen