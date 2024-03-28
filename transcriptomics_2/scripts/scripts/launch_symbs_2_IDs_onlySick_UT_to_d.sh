#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=convSymb2Entr
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100G
#SBATCH -c 1
#SBATCH --time=00-00:10:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/cureMILS/transcriptomics_2/scripts/output_symbs2IDs.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/cureMILS/transcriptomics_2/scripts/error_symbs2IDs.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal


module load lang/R/4.0.5-foss-2020b

# Variables for the pipeline
########################################################################################################################
DE_file="/home/users/gsantamaria/projects/cureMILS/transcriptomics_2/results/DE/DEG_no_drug_to_drug_just_mutation.xlsx"
#outNameNEAT="DEG_no_drug_to_drug_05.xlsx_neat.csv"
conv2ENSEMBL=true
#FCFile="~/projects/cureMILS/proteomics/data/FunCoup/FC5.0_H.sapiens_full"

# If the datasets and dataformat are not installed this script will install them
bash install_ncbi_CLTools.sh

# Create the dataframe with ENSEMBL IDs
if [ $conv2ENSEMBL = true ]
then
    Rscript symbs_2_IDs.R $DE_file --whatOutID ENTREZID --remResDups --saveIDsDF
    csvExt=".csv"
    xlsxExt=".xlsx"
    DE4NEAT="${DE_file/$csvExt/}"
    DE4NEAT="${DE_file/$xlsxExt/}"
    DE4NEAT=$DE4NEAT"_ENSEMBL_IDs.csv"
else
    DE4NEAT=$DE_file
fi