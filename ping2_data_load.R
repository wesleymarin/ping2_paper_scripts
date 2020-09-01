library(data.table,lib.loc = '~/R/x86_64-redhat-linux-gnu-library/3.4/')
library(ggplot2)
library(stringr)
library(methods)
library(plotly,lib.loc='~/R/x86_64-redhat-linux-gnu-library/3.4/')
library(gtools)


kirLocusFeatureNameList <- list()
kirLocusFeatureNameList[['KIR2DL1']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL2']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL3']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL4']] <- c("5UTR","E1","I1","E2","I2","E3","I3/4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DL5']] <- c("5UTR","E1","I1","E2","I2","E3","I3/4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DP1']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS1']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS2']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS3']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS4']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR2DS5']] <- c("5UTR","E1","I1","E2","I2","PE3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DL1']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DL2']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DL3']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5/6","E7","I7","E8","I8","E9","3UTR")
kirLocusFeatureNameList[['KIR3DP1']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","3UTR")
kirLocusFeatureNameList[['KIR3DS1']] <- c("5UTR","E1","I1","E2","I2","E3","I3","E4","I4","E5","I5","E6","I6","E7","I7","E8","I8","E9","3UTR")

## -- PING workflow read in
ping.1.typeDF <- read.csv('/home/wmarin/PING_projects/PING2_scripts/PING2_paper_results/ping2_results/pingAlleleCalls.final.csv',
                          stringsAsFactors=F,check.names=F,row.names=1)
ping.1.typeDF <- ping.1.typeDF[ apply(ping.1.typeDF, 1, function(x){ !all(is.na(x)) } ), ]


## -- ITER workflow read in
iter.2.typeDF <- read.csv('/home/wmarin/PING_projects/PING2_scripts/PING2_paper_results/ping2_results/iterAlleleCalls.final.csv', 
                          stringsAsFactors = F, check.names = F, row.names=1)

iter.2.typeDF <- iter.2.typeDF[ apply(iter.2.typeDF, 1, function(x){ !all(is.na(x)) } ), ]

#iter.2.typeDF$'KIR2DL23' <- ''
iter.2.typeDF$'KIR2DL23' <- unlist(apply(iter.2.typeDF[,c('KIR2DL2','KIR2DL3','KIR2DL23')],1,function(x){
  kir2DL2Call <- x[1]
  kir2DL3Call <- x[2]
  kir2DL23Call <- x[3]
  
  if( is.na(kir2DL2Call) & is.na(kir2DL3Call) & is.na(kir2DL23Call) ){
    return(NA)
  }else if( is.na(kir2DL2Call) ){
    return( kir2DL3Call )
  }else if( is.na(kir2DL3Call) ){
    return( kir2DL2Call )
  }else if( kir2DL2Call == 'new' | kir2DL3Call == 'new' & is.na(kir2DL23Call) ){
    return('new')
  }else{
    return(kir2DL23Call)
  }
  
  #else{
  #  kir2DL2Vect <- unlist(tstrsplit(kir2DL2Call,' '))
  #  kir2DL3Vect <- unlist(tstrsplit(kir2DL3Call,' '))
  #  
  #  combinationDF <- expand.grid(kir2DL2Vect,kir2DL3Vect)
  #  kir2DL23Call <- paste0(apply(combinationDF,1,paste0,collapse='+'), collapse=' ')
  #  return(kir2DL23Call)
  #}
}))


## -- FILTER workflow read in
filter.2.typeDF <- read.csv('/home/wmarin/PING_projects/PING2_scripts/PING2_paper_results/ping2_results/filterAlleleCalls.final.csv',
                            stringsAsFactors = F, check.names=F, row.names=1 )

filter.2.typeDF <- filter.2.typeDF[ apply(filter.2.typeDF, 1, function(x){ !all(is.na(x)) } ), ]



## -- KFF results
#kff.2.countDF <- read.csv('/home/wmarin/PING_projects/PING2_scripts/PING2_paper_results/kff_ambiguity_reduction_data/kffCountFrame.csv',
#                          stringsAsFactors = F, check.names=F, row.names=1)


## -- PING/FILTER 2DS3/5 allele split
## KIR2DS5 isolation in PING typings


post.filter.addKIR2DS5 <- function( typeDF ){
  for(i in 1:nrow(typeDF)){
    
    index2DS5Int <- which(grepl('KIR2DS5',typeDF[i,],fixed=T))
    
    if( length(index2DS5Int) == 0 ){
      newBool <- grepl('unresolved',typeDF[i,'KIR2DS35'],fixed=T)
      
      if(!newBool){
        next
      }
      type2DS5Str <- 'unresolved'
    }else{
      type2DS5Str <- typeDF[i,'KIR2DS35'][[1]]
    }
    
    vectType2DS5 <- tstrsplit(type2DS5Str,' ',fixed=T)
    listType2DS5 <- list()
    j <- 1
    for(sub2DS5Str in vectType2DS5){
      
      ## KIR2DS3 split
      if(any(grepl('KIR2DS3', sub2DS5Str, fixed=T))){
        splitVect <- strsplit(sub2DS5Str,'+',fixed=T)[[1]]
        sub2DS5Str <- splitVect[grepl('KIR2DS5',splitVect,fixed=T)]
        listType2DS5[[j]] <- sub2DS5Str
        j <- j+1
      }else{
        listType2DS5[[j]] <- sub2DS5Str
        j <- j+1
      }
      
    }
    
    typeDF[i,'KIR2DS5'] <- paste0(listType2DS5, collapse=' ')
  }
  return(typeDF)
}

post.filter.addKIR2DS3 <- function( typeDF ){
  
  ## KIR2DS3 isolation in PING typings
  for(i in 1:nrow(typeDF)){
    
    index2DS3IntVect <- which(grepl('KIR2DS3',typeDF[i,],fixed=T))
    
    if( length(index2DS3IntVect) == 0 ){
      
      newBool <- any(grepl('unresolved',typeDF[i,c('KIR2DS3','KIR2DS35')],fixed=T))
      
      if( !newBool ){
        next
      }else{
        type2DS3Str <- 'unresolved'
      }
    }else{
      type2DS3Str <- typeDF[i,index2DS3IntVect,drop=T][[1]]
    }
    
    if(length(index2DS3IntVect) > 1){
      cat('\n',type2DS3Str)
    }
    typeDF[i,'KIR2DS3'] <- type2DS3Str
  }
  return(typeDF)
}

ping.1.typeDF$KIR2DS5 <- ''
ping.1.typeDF <- post.filter.addKIR2DS5( ping.1.typeDF )
ping.1.typeDF <- post.filter.addKIR2DS3( ping.1.typeDF )


filter.2.typeDF$KIR2DS5 <- ''
filter.2.typeDF <- post.filter.addKIR2DS5( filter.2.typeDF )
filter.2.typeDF <- post.filter.addKIR2DS3( filter.2.typeDF )


## Creating the sample object class
sample <- setRefClass("sample",
                      fields=list(name='character',
                                  haploResult='list',
                                  trimHaploResult='list',
                                  pingResult='list',
                                  trimPingResult='list',
                                  validResult='list',
                                  trimValidResult='list',
                                  copyResult='list',
                                  interResult='list'))

valid.ping2.buildSampleList <- function(iter.2.typeDF, filter.2.typeDF, validTypeDF, copyDF){
  haploSampleVect <- rownames(iter.2.typeDF)
  pingSampleVect <- rownames(filter.2.typeDF)
  validSampleVect <- rownames(validTypeDF)
  copySampleVect <- rownames(copyDF)
  
  sampleNameVect <- unique(c(haploSampleVect, pingSampleVect, validSampleVect, copySampleVect))
  
  output.sampleList <- list()
  for(sampleNameStr in sampleNameVect){
    cat('\n',sampleNameStr)
    
    ## Checking haplo results
    if(sampleNameStr %in% haploSampleVect){
      ## Pull out typings
      haploType <- iter.2.typeDF[sampleNameStr,]
      
      ## If oll results for this sample are NA, then return NA
      con1 <- all(is.na(haploType))
      
      if(con1){
        
        sampleHaploType <- list()
        
      }else{
        
        sampleHaploType <- list()
        ## Otherwise iterate over each typed locus
        for(currentLocus in names(haploType)){
          cat('\t',currentLocus)
          ## Pull out the typing for the current locus
          currentLocusType <- haploType[[currentLocus]][[1]]
          
          ## If the current locus typing is NA, then record NA
          if(is.na(currentLocusType)){
            sampleHaploType[[currentLocus]] <- NA
            next
          }
          
          if( grepl('unresolved', currentLocusType) ){
            currentLocusType <- 'new'
          }
          
          ## Split ambiguous typings
          currentLocusTypeVect <- strsplit(currentLocusType,' ',fixed=T)[[1]]
          
          currentLocusTypeVect <- unique(currentLocusTypeVect)
          
          ## Expanding typings that need to be expanded
          if(all(!grepl('+',currentLocusTypeVect,fixed=T)) & length(currentLocusTypeVect)>1){
            
            combMat <- combinations(length(currentLocusTypeVect),2,currentLocusTypeVect,repeats.allowed = T)
            currentLocusTypeVect <- apply(combMat, 1, paste0, collapse='+')
          }
          
          ## Initialize a dataframe for storing allele typings, 1 row for each ambiguity
          alleleTypeFrame <- data.frame(matrix('',nrow=length(currentLocusTypeVect),ncol=2),
                                        stringsAsFactors = F)
          colnames(alleleTypeFrame) <- c('allele1','allele2')
          
          i <- 0
          ## For each ambiguity, record the typings in the dataframe
          for(singleType in currentLocusTypeVect){
            i <- i+1
            con2 <- grepl('+',singleType,fixed=T)
            
            if(con2){
              hetTypeVect <- unlist(tstrsplit(singleType,'+',fixed=T))
              
              if(length(hetTypeVect) > 2){
                alleleTypeFrame <- NA
              }else{
                alleleTypeFrame[as.character(i),] <- hetTypeVect
              }
            }else{
              alleleTypeFrame[as.character(i),] <- singleType
            }
          }
          
          ## Save the dataframe to the type list
          sampleHaploType[[currentLocus]] <- alleleTypeFrame
        }
      }
    }else{
      sampleHaploType <- list()
    }
    
    ## Checking ping results
    if(sampleNameStr %in% pingSampleVect){
      ## Pull out typings
      pingType <- filter.2.typeDF[sampleNameStr,]
      
      ## Check for bad sample
      con1 <- all(is.na(pingType) | pingType == '')
      
      if(con1){
        samplePingType <- list()
      }else{
        
        samplePingType <- list()
        ## Otherwise iterate over each typed locus
        for(currentLocus in names(pingType)){
          
          ## Pull out the typing for the current locus
          currentLocusType <- pingType[[currentLocus]][[1]]
          
          ## If the current locus typing is NA, then record NA
          if(is.na(currentLocusType) | currentLocusType == ''){
            samplePingType[[currentLocus]] <- NA
            next
          }
          
          if( grepl('unresolved', currentLocusType) ){
            currentLocusType <- 'new'
          }
          
          ## Split ambiguous typings
          currentLocusTypeVect <- strsplit(currentLocusType,' ',fixed=T)[[1]]
          
          ## Expanding typings that need to be expanded
          if(all(!grepl('+',currentLocusTypeVect,fixed=T)) & length(currentLocusTypeVect)>1){
            
            combMat <- combinations(length(currentLocusTypeVect),2,currentLocusTypeVect,repeats.allowed = T)
            currentLocusTypeVect <- apply(combMat, 1, paste0, collapse='+')
          }
          
          ## Initialize a dataframe for storing allele typings, 1 row for each ambiguity
          alleleTypeFrame <- data.frame(matrix('',nrow=length(currentLocusTypeVect),ncol=2),
                                        stringsAsFactors = F)
          colnames(alleleTypeFrame) <- c('allele1','allele2')
          
          i <- 0
          ## For each ambiguity, record the typings in the dataframe
          for(singleType in currentLocusTypeVect){
            i <- i+1
            con2 <- grepl('+',singleType,fixed=T)
            
            if(con2){
              hetTypeVect <- unlist(tstrsplit(singleType,'+',fixed=T))
              
              if(length(hetTypeVect) > 2){
                stop('ping')
              }
              alleleTypeFrame[as.character(i),] <- hetTypeVect
            }else{
              alleleTypeFrame[as.character(i),] <- singleType
            }
          }
          
          ## Save the dataframe to the type list
          samplePingType[[currentLocus]] <- alleleTypeFrame
        }
      }
    }else{
      samplePingType <- list()
    }
    
    ## Checking valid results
    if(sampleNameStr %in% validSampleVect){
      ## Pull out typings
      validType <- validTypeDF[sampleNameStr,]
      
      ## Check for bad sample
      con1 <- all(is.na(validType) | validType == '')
      
      if(con1){
        sampleValidType <- list()
      }else{
        
        sampleValidType <- list()
        ## Otherwise iterate over each typed locus
        for(currentLocus in names(validType)){
          
          ## Pull out the typing for the current locus
          currentLocusType <- validType[[currentLocus]][[1]]
          
          ## If the current locus typing is NA, then record NA
          if(is.na(currentLocusType) | 
             currentLocusType == '' | 
             any(grepl('something',currentLocusType)) | 
             any(grepl('?',currentLocusType,fixed=T))){
            sampleValidType[[currentLocus]] <- NA
            next
          }
          
          ## Split ambiguous typings
          currentLocusTypeVect <- strsplit(currentLocusType,' ',fixed=T)[[1]]
          
          ## Initialize a dataframe for storing allele typings, 1 row for each ambiguity
          alleleTypeFrame <- data.frame(matrix('',nrow=length(currentLocusTypeVect),ncol=2),
                                        stringsAsFactors = F)
          colnames(alleleTypeFrame) <- c('allele1','allele2')
          
          i <- 0
          ## For each ambiguity, record the typings in the dataframe
          for(singleType in currentLocusTypeVect){
            i <- i+1
            con2 <- grepl('+',singleType,fixed=T)
            
            if(con2){
              hetTypeVect <- unlist(tstrsplit(singleType,'+',fixed=T))
              
              ## Excluding copy 3 results
              if(length(hetTypeVect) > 2){
                alleleTypeFrame <- NA
              }else{
                alleleTypeFrame[as.character(i),] <- hetTypeVect
              }
            }else{
              alleleTypeFrame[as.character(i),] <- singleType
            }
          }
          
          ## Save the dataframe to the type list
          sampleValidType[[currentLocus]] <- alleleTypeFrame
        }
      }
    }else{
      sampleValidType <- list()
    }
    
    ## Checking copy results
    if(sampleNameStr %in% copySampleVect){
      ## Pull out typings
      copyType <- copyDF[sampleNameStr,]
      
      ## Check for bad sample
      con1 <- all(is.na(copyType) | copyType == '')
      
      if(con1){
        sampleCopyType <- list()
      }else{
        sampleCopyType <- list()
        for(currentLocus in names(copyType)){
          currentCopy <- copyType[[currentLocus]]
          
          if(is.numeric(currentCopy)){
            sampleCopyType[[currentLocus]] <- currentCopy
          }else{
            sampleCopyType[[currentLocus]] <- NA
          }
        }
      }
    }else{
      sampleCopyType <- list()
    }
    
    output.sampleList[[sampleNameStr]] <- sample(name=sampleNameStr,
                                                 haploResult=sampleHaploType,
                                                 pingResult=samplePingType,
                                                 validResult=sampleValidType,
                                                 copyResult=sampleCopyType)
  }
  
  return(output.sampleList)
}


## Build up a list of sample objects that contain all results
#ping2.sampleList <- ping2.buildSampleList( iter.2.typeDF, filter.2.typeDF, validTypeDF, copyDF )

# sample list for ping2_allele_concordance.R
ping2.sampleList <- valid.ping2.buildSampleList( iter.2.typeDF, ping.1.typeDF, validTypeDF, copyDF )

