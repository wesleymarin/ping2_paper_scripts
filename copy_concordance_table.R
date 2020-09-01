

setwd('/home/LAB_PROJECTS/PING2_PAPER/PING2')
source('Resources/general_functions.R')
source('Resources/extractor_functions.R')
source('ping_gc.R')

# Initialization variables ------------------------------------------------
rawFastqDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_script_results/extractedFastq/'
#sequenceDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/sequence_data/'
fastqPattern <- '_KIR_'
threads <- 12
resultsDirectory <- '/home/LAB_PROJECTS/PING2_PAPER/3_script_results/'


# Preparation -------------------------------------------------------------
# Build up a list of sample objects
sampleList <- general.paired_sample_objects(rawFastqDirectory, fastqPattern, resultsDirectory)


# PING2 extractor ---------------------------------------------------------
cat('\n\n----- Moving to PING2 KIR extraction -----')
# Define the extracted fastq directory
extractedFastqDirectory <- file.path(resultsDirectory,'extractedFastq')
# Run PING2 extractor
sampleList <- extractor.run(sampleList,threads,extractedFastqDirectory,forceRun=F)


# PING2 gene content and copy number --------------------------------------
cat('\n\n----- Moving to PING2 gene content and copy determination -----')
sampleList <- ping_gc(sampleList=sampleList,threads=threads,resultsDirectory=resultsDirectory,predictCopy=F)



# copy comparison ---------------------------------------------------------

validCopyDF <- read.csv('../validation_data/INDIGO_controls_copyNumber.csv',stringsAsFactors=F,row.names=1,check.names=F)
validCopyDF$KIR2DL5 <- validCopyDF$KIR2DL5A + validCopyDF$KIR2DL5B

ping2CopyDF <- read.csv('../3_script_results/manualCopyNumberFrame.csv',stringsAsFactors=F,row.names=1,check.names=F)
rownames(ping2CopyDF) <- tstrsplit(rownames(ping2CopyDF),'_',fixed=T)[[1]]

ping2RatioDF <- read.csv('../3_script_results/locusRatioFrame.csv',stringsAsFactors=F,row.names=1,check.names=F)
rownames(ping2RatioDF) <- tstrsplit(rownames(ping2RatioDF),'_',fixed=T)[[1]]

inBothRowVect <- intersect(rownames(validCopyDF), rownames(ping2CopyDF))
validCopyDF <- validCopyDF[inBothRowVect,]
ping2CopyDF <- ping2CopyDF[inBothRowVect,]

inBothColVect <- intersect(colnames(validCopyDF), colnames(ping2CopyDF))
validCopyDF <- validCopyDF[,inBothColVect]
ping2CopyDF <- ping2CopyDF[,inBothColVect]


for(currentLocus in colnames(ping2CopyDF)){
  cat('\n',currentLocus)
  copyLevelVect <- sort(unique(ping2CopyDF[inBothRowVect,currentLocus]))
  
  for(copyLevel in copyLevelVect){
    cat('\n\t',copyLevel)
    copyLevelPingIDVect <- inBothRowVect[ping2CopyDF[inBothRowVect,currentLocus] == copyLevel]
    copyLevelValidIDVect <- inBothRowVect[validCopyDF[inBothRowVect,currentLocus] == copyLevel]
    
    # Identify recombinant samples in the validation dataset
    recombValidIDVect <- inBothRowVect[validCopyDF[inBothRowVect,currentLocus] == 'R']
    
    # Remove recomb samples from ping2 results
    copyLevelPingIDVect <- setdiff(copyLevelPingIDVect, recombValidIDVect)
    
    matchIDVect <- intersect(copyLevelPingIDVect, copyLevelValidIDVect)
    mismatchIDVect <- setdiff(copyLevelPingIDVect, copyLevelValidIDVect)
    
    totalVectLen <- length(copyLevelPingIDVect)
    matchLen <- length(matchIDVect)
    mismatchLen <- length(mismatchIDVect)
    
    cat('\tMatch:',matchLen,'/',totalVectLen)
    
    if(mismatchLen > 0){
      cat('\tMismatched IDs:',paste0(mismatchIDVect,collapse=' '))
    }
  }
}






