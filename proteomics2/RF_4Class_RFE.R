################################################################################
# CureMILS project:                                                            #
# PROTEOMICS: Random Forest with RFE doing a 4 category classification, with.  #
# the purpose of doing variable selection                                      #
################################################################################
if(!require("caret", quietly = T)){
        install.packages("caret",
                         repos='http://cran.us.r-project.org')
}
library(caret)
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)

# Terminal argument parser
################################################################################
parser <- arg_parser("Performs RFE with Random Forest with the aim of doing feature selection for the classification in 4 categories: sick UT, sick d, ctrl UT and ctrl d. Saves result as an .Rdata file")

parser <- add_argument(parser = parser,
                       arg = c("input", "--sampInfo", "--outFile"),
                       help = c("Input omic dataset.",
                                "Dataframe with sample information of the phenotypes and treatments. Needs to have at least sample, disease and treatment columns.",
                                "Name of the output .RData file"),
                       flag = c(F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
inFile <- parsed$input
sampInfoFile <- parsed$sampInfo
outFile <- parsed$outFile

outDir <- paste0(dirname(outFile), "/")

if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}


# Load functions and data
################################################################################

# Load the proteomics table
proteom <- read.csv(inFile,
                    header = T,
                    row.names = 1)

# Load the dataframe of sample information
sample_info <- read.csv(sampInfoFile,
                        header = T,
                        row.names = 1)

# Do the RF
################################################################################
groups <- paste(sample_info$disease,
                sample_info$treatment, sep = "_")[match(rownames(proteom),
                                                        sample_info$sample)]

proteom4RF <- proteom
proteom4RF$groups <- as.factor(groups)
set.seed(12345)
inTrain <- createDataPartition(y = proteom4RF$groups, p = .66, list = FALSE)
training <- proteom4RF[inTrain, ]
testing <- proteom4RF[-inTrain, ]

sample_info[match(rownames(proteom4RF)[inTrain[, 1]], sample_info$sample), ]

rfe_ctrl <- rfeControl(functions = rfFuncs,
                       method = "repeatedcv",
                       number = 5,
                       repeats = 2,
                       allowParallel = T)


# Create a vector of sizes to test, as testing every size would 
# take ages to run
nFeats <- ncol(training) - 1

sizes2Tst <- c(1:((nFeats - (nFeats %% 10))/ 10),
               seq(from = ((nFeats - (nFeats %% 10))/ 10) + 10,
                   to = ((nFeats - (nFeats %% 10))/ 10) + 1010,
                   by = 10),
               seq(from = ((nFeats - (nFeats %% 10))/ 10) + 1110,
                   to = nFeats - (nFeats %% 100),
                   by = 100),
                   nFeats)

print("Sizes to be tested in RFE:")
print(sizes2Tst)

print("Running RFE...")
set.seed("1273627")
rfe_RF <- rfe(groups ~ .,
              data = training,
              metric = "Accuracy",
              sizes = sizes2Tst,
              rfeControl = rfe_ctrl)

save(rfe_RF, file = outFile)

print(sprintf("%s saved at %s",
              basename(outFile),
              outDir))