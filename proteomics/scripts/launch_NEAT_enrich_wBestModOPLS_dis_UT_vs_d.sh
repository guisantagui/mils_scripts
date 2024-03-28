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
#SBATCH --time=00-03:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/cureMILS/proteomics/scripts/output_NEAT_enrich_dis_UT_vs_d_bestMod.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/cureMILS/proteomics/scripts/error_NEAT_enrich_dis_UT_vs_d_bestMod.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal


module load lang/R/4.0.5-foss-2020b

# Variables for the pipeline
########################################################################################################################
DE_file="/home/users/gsantamaria/projects/cureMILS/proteomics/results/DE/proteom_log2_OPLS_dis_UT_vs_d_bestMod_resSign.csv"
outNameNEAT="proteom_log2_sick_UT_vs_d_oplsBestMod_neat.csv"
conv2ENSEMBL=true
FCFile="~/projects/cureMILS/proteomics/data/FunCoup/FC5.0_H.sapiens_full"

# If the datasets and dataformat are not installed this script will install them
bash install_ncbi_CLTools.sh

# Create the dataframe with ENSEMBL IDs
if [ $conv2ENSEMBL = true ]
then
    Rscript symbs_2_IDs.R $DE_file --whatOutID ENSEMBL --remResDups
    csvExt=".csv"
    DE4NEAT="${DE_file/$csvExt/}"
    DE4NEAT=$DE4NEAT"_ENSEMBL_IDs.csv"
else
    DE4NEAT=$DE_file
fi

# Run NEAT
Rscript NEAT_enrich.R --DE $DE4NEAT --outName $outNameNEAT --FCFile $FCFile
