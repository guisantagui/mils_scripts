if(!require("ropls", quietly = T)) BiocManager::install("ropls", update = F)
library(ropls)
if(!require(argparser, quietly = T)){
        install.packages("argparser", repos='http://cran.us.r-project.org')
}
library(argparser)

################################################################################
# CureMILS project:                                                            #
# PROTEOMICS: OPLS-DA RFE. Given a matrix with rows being samples and columns  #
# features, and the trait to classify, performs Recursive Feature Elimination  #
# to find the set of features that maximizes the Q2 of the model.              #
################################################################################

# Terminal argument parser
################################################################################
parser <- arg_parser("RFE for omic data classification")

parser <- add_argument(parser = parser,
                       arg = c("input",
                               "--sampMetadata",
                               "--respVar",
                               "--outDir",
                               "--splitTrainTest",
                               "--filtSampsByDisPhen",
                               "--filtSampsByTreatment"),
                       help = c("Input omic dataset.",
                                "Metadata with information of each sample (disease, treatment, etc).",
                                "Response variable to be used in the classification (disease or treatment).",
                                "Directory where all the produced files 
                                will be saved.",
                                "Training:testing proportion. If not used the model will be fitted in all the dataset.",
                                "Filter the samples to keep either ctrl or sick ones. If not used no filtering will be done.",
                                "Filter the samples to keep either UT or d ones. If not used no filtering will be done."),
                       flag = c(F,
                                F, F, F, F, F, F),
                       default = list("input" = "proteom_imp.csv",
                                      "--respVar" = "treatment",
                                      "--sampMetadata" = "/Users/guillem.santamaria/Documents/postdoc/comput/cureMILS/proteomics/sample_info.csv",
                                      "--outDir" = "/Users/guillem.santamaria/Documents/postdoc/comput/cureMILS/proteomics/results/sup_analysis/OPLSDA_RFE/",
                                      "--splitTrainTest" = "no",
                                      "--filtSampsByDisPhen" = "no",
                                      "--filtSampsByTreatment" = "no"))

parsed <- parse_args(parser)


# Functions
################################################################################

# Computes several iteration of OPLS-DA removing the worse variable for the 
# prediction, evaluated by VIP, and stores R2Y and Q2Y for each iteration in a 
# dataframe, which is the output
opls_RFE <- function(DF, y, orthoI = 3){
        oplsFit <- opls(DF,
                        y = y,
                        orthoI = orthoI,
                        permI = 2,
                        fig.pdfC = "none")
        r2y <- oplsFit@summaryDF$`R2Y(cum)`
        q2y <- oplsFit@summaryDF$`Q2(cum)`
        
        r2y_vec <- r2y
        q2y_vec <- q2y
        print(dim(DF))
        nFeats <- ncol(DF)
        nFeatsVec <- nFeats
        featsInMod <- paste(colnames(DF), collapse = ", ")
        featsInModVec <- featsInMod
        while(nFeats >= (orthoI + 4)){
                worseFeat <- sort(oplsFit@vipVn, decreasing = F)[1]
                print(worseFeat)
                DF <- DF[, !grepl(names(worseFeat), colnames(DF))]
                print(dim(DF))
                oplsFit <- opls(DF,
                                y = y,
                                orthoI = orthoI,
                                permI = 2,
                                fig.pdfC = "none")
                r2y <- oplsFit@summaryDF$`R2Y(cum)`
                q2y <- oplsFit@summaryDF$`Q2(cum)`
                r2y_vec <- c(r2y_vec, r2y)
                q2y_vec <- c(q2y_vec, q2y)
                nFeats <- ncol(DF)
                nFeatsVec <- c(nFeatsVec, nFeats)
                featsInMod <- paste(colnames(DF), collapse = ", ")
                featsInModVec <- c(featsInModVec, featsInMod)
        }
        statsDF <- data.frame(nFeats = nFeatsVec,
                              r2y = r2y_vec,
                              q2y = q2y_vec,
                              featsInMod = featsInModVec)
        return(statsDF)
}

# Given the RFE dataframe generated with previous function, obtains the best
# model
getBestOplsMod <- function(RFE_DF, protDF, orthoI = 3){
        bestProts <- RFE_DF$featsInMod[which.max(RFE_DF$q2y)]
        print(bestProts)
        bestProts <- strsplit(bestProts, split = ", ")[[1]]
        protDF_bestProts <- protDF[, bestProts]
        bestOplsMod <- opls(protDF_bestProts,
                            factor(sample_info$treatment[match(rownames(protDF_bestProts),
                                                               sample_info$sample)],
                                   levels = sort(unique(sample_info$treatment),
                                                 decreasing = T)),
                            orthoI = orthoI,
                            predI = 1,
                            permI = 200)
        return(bestOplsMod)
}

# Directory stuff
################################################################################
inFile <- parsed$input
metDatFile <- parsed$sampMetadata
respVar <- parsed$respVar
outDir <- parsed$outDir
splitTrain <- parsed$splitTrainTest
filtByDis <- parsed$filtSampsByDisPhen
filtByTreat <- parsed$filtSampsByTreatment

inName <- sub("\\.\\w+$", "", basename(inFile))
filtDisName <- gsub("no", "", filtByDis)
filtTreatName <- gsub("no", "", filtByTreat)
trainSplitName <- gsub("no", "", as.character(splitTrain))
if(nchar(trainSplitName) > 0){
        trainSplitName <- paste0(trainSplitName,
                                 "_inTrain")
}
respVarName <- paste(respVar, "class", sep = "_")

outName <- sprintf("%s%s_%s%s_%s_%s",
                   outDir,
                   inName,
                   filtByDis,
                   filtByTreat,
                   trainSplitName,
                   respVarName)


if(!dir.exists(outDir)){
        dir.create(outDir, recursive = T)
}

# Load data and filter it if required
################################################################################

DF <- read.csv(inFile, row.names = 1)
sample_info <- read.csv(metDatFile, row.names = 1)

if(filtByDis == "sick"){
        DF <- DF[rownames(DF) %in% sample_info$sample[sample_info$disease ==
                                                              "sick"], ]
}else if(filtByDis == "ctrl"){
        DF <- DF[rownames(DF) %in% sample_info$sample[sample_info$disease == 
                                                              "ctrl"], ]
}

if(filtByTreat == "UT"){
        DF <- DF[rownames(DF) %in% sample_info$sample[sample_info$treatment ==
                                                              "UT"], ]
}else if(filtByTreat == "d"){
        DF <- DF[rownames(DF) %in% sample_info$sample[sample_info$treatment ==
                                                              "d"], ]
}

samp_resp <- sample_info[match(rownames(DF),
                               sample_info$sample),
                         respVar]

# Analysis
################################################################################

print(dim(DF))
print(rownames(DF))

# If split train proportion is specified, split the dataset in training and testing.
if(splitTrain != "no"){
        splitTrain <- as.numeric(splitTrain)
        nToSamp <- floor(table(samp_resp) * splitTrain)
        print(nToSamp)
        classes <- names(nToSamp)
        print(classes)
        inTrain <- c()
        for(i in seq_along(classes)){
                cl <- classes[i]
                numCl <- nToSamp[i]
                print(cl)
                print(numCl)
                sampIdxs <- which(rownames(DF) %in% sample_info$sample[sample_info[, respVar] == cl])
                set.seed(122 + i)
                sampIdxs_4Train <- sample(sampIdxs, size = numCl)
                inTrain <- c(inTrain, sampIdxs_4Train)
        }
        print(inTrain)
        DF_train <- DF[inTrain, ]
        DF_test <- DF[-inTrain, ]
        print("Samples used for training:")
        print(rownames(DF_train))
        samp_resp <- sample_info[match(rownames(DF_train),
                                       sample_info$sample),
                                 respVar]
        DF <- list(train = DF_train,
                   test = DF_test)
}else{
        DF <- list(train = DF)
}

# Order the factors of the response variable in order to have them
# either ctrl-->sick or UT-->d depending on what was chosen as 
# response variable.
if(respVar == "disease"){
        samp_resp <- factor(samp_resp,
                            levels = sort(unique(samp_resp)))
}else if(respVar == "treatment"){
        samp_resp <- factor(samp_resp,
                            levels = sort(unique(samp_resp),
                                          decreasing = T))
}

rfeStatsFile <- sprintf("%s_statsDF.csv", outName)
if(!file.exists(rfeStatsFile)){
        RFE_DF <- opls_RFE(DF$train,
                   y = samp_resp)
        
        write.csv(RFE_DF, rfeStatsFile)
        print(sprintf("%s saved at %s.",
              basename(rfeStatsFile),
              dirname(rfeStatsFile)))
}else{
        RFE_DF <- read.csv(rfeStatsFile, row.names = 1)
        print(sprintf("%s loaded from %s.",
              basename(rfeStatsFile),
              dirname(rfeStatsFile)))
}


pdf(file = sprintf("%s.pdf", outName))
bestMod <- getBestOplsMod(RFE_DF, DF$train)
dev.off()

print(sprintf("%s_bestMod.pdf saved at %s.",
              basename(outName),
              dirname(outName)))

save(bestMod,
     file = sprintf("%s_bestMod.RData",
                    outName))

print(sprintf("%s_bestMod.RData saved at %s.",
              basename(outName),
              dirname(outName)))               


if(length(DF) == 2){
        print("Samples used for testing:")
        print(rownames(DF$test))
        test_classes <- sample_info[match(rownames(DF$test),
                                          sample_info$sample),
                                    respVar]

        DF_test <- DF$test
        DF_test <- DF_test[, colnames(DF_test) %in% names(bestMod@vipVn)]
        pred <- predict(bestMod, DF_test)
        print(pred)
        print(test_classes)
        save(DF,
             file = sprintf("%s%s_trainTest.RData",
                            outDir,
                            as.character(splitTrain)))

        # Compute NIR
        nir <- max(table(test_classes)/length(test_classes))
        
        # Get Confusion Matrix
        confMat <- table(data.frame(predicted = pred,
                                    real = test_classes))

        print(confMat)
        
        # Calculate p-value with successes and errors, compared to
        # NIR
        binomPVal = binom.test(c(confMat[1, 1] + confMat[2, 2],
                                 confMat[2, 1] + confMat[1, 2]),
                               alternative = "greater",
                               p = nir)$p.value

        print(sprintf("NIR: %s", as.character(nir)))
        print(sprintf("p-value: %s", as.character(binomPVal)))
}