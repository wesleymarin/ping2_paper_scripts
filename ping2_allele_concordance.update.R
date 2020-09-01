library(data.table,lib.loc = '~/R/x86_64-redhat-linux-gnu-library/3.4/')
library(ggplot2)
library(stringr)
library(methods)
library(plotly,lib.loc='~/R/x86_64-redhat-linux-gnu-library/3.4/')
library(gtools)

setwd('/home/wmarin/PING_projects/PING2_scripts/')

source('/home/wmarin/PING_projects/PING2/Resources/gc_functions.R')
#source('ping2_paper_INDIGO_results_read.R')

resultsDirectory <- '/home/wmarin/PING_projects/PING2_scripts/PING2_paper_results/'

kirLocusList <- c('KIR3DP1','KIR2DS5','KIR2DL3','KIR2DP1',
                  'KIR2DS3','KIR2DS2','KIR2DL4','KIR3DL3',
                  'KIR3DL1','KIR3DS1','KIR2DL2','KIR3DL2','KIR2DS4','KIR2DL1', 'KIR2DS1', 'KIR2DL5')

gcResourceDirectory <- normalizePath('/home/wmarin/PING_projects/PING2/Resources/gc_resources', mustWork = T)
kirReferenceAlignedFasta <- normalizePath(file.path(gcResourceDirectory, 'filled_kir_reference', 'KIR_gen_onelines_filled.fasta'), mustWork=T)

## Read in the aligned fasta and convert it into a list of dataframes (a dataframe for each locus)
cat('\nLoading SNP dfs...\t')
kirAlleleDFList <- list()
for(locusName in kirLocusList){
  cat('',locusName)
  # Load the SNP df
  locusSnpDF <- read.csv( normalizePath(
    file.path(gcResourceDirectory,
              paste0('ping2_alleles/', locusName, '_alleleSNPs.csv')
    )), check.names=F, row.names=1, stringsAsFactors = F)
  
  # Pull out the exon labels
  exonFeatVect <- grep('E',kirLocusFeatureNameList[[locusName]], value=T, fixed=T)
  exonFeatVect <- grep('PE',exonFeatVect,value=T, invert=T, fixed=T)
  
  # Subset the SNP df by exonic positions
  locusSnpDF <- locusSnpDF[, tstrsplit( colnames(locusSnpDF), '_', fixed=T )[[1]] %in% exonFeatVect ]
  kirAlleleDFList[[locusName]] <- locusSnpDF
}
#cat('\nReading in known KIR alleles...')
#kirAlleleDFList <- read.kir_allele_dataframe_from_reference_fasta(kirReferenceAlignedFasta, kirLocusList)
cat('\nFinished.')


## Set up the copy number results path
#haploTypePath <- '/home/common_arse/INDIGO/4_post_analysis/all_indigo_type_frame.csv'
copyPath <- '/home/common_arse/INDIGO/4_post_analysis/INDIGO_copyNumberFrame_2018dec12.csv'
#haploDistancePath <- '/home/common_arse/INDIGO/4_post_analysis/all_indigo_distance_frame.csv'


## Read in data
#haploTypeDF <- read.csv(haploTypePath, stringsAsFactors = F, check.names = F)
copyDF <- read.csv(copyPath, stringsAsFactors = F, check.names= F, row.names = 1)
#distanceDF <- read.csv(haploDistancePath,stringsAsFactors=F,check.names=F,row.names=1)
probeRatioDF <- read.csv('3DL2RatioDF.csv',stringsAsFactors = F,check.names=F,row.names=1)
rownames(probeRatioDF) <- tstrsplit(rownames(probeRatioDF),'_',fixed=T)[[1]]
probeConvDF <- read.csv('kffResolver.v2.csv',stringsAsFactors = F,check.names=F,row.names=1)

probeConvList <- list()
probeConvList[['ratioA']] <- unlist(tstrsplit(probeConvDF['probe4','allele_matches'],' ',fixed=T))
probeConvList[['ratioB']] <- unlist(tstrsplit(probeConvDF['probe7','allele_matches'],' ',fixed=T))
probeConvList[['ratioA']] <- unique(unlist(lapply(probeConvList$ratioA,kir.allele_resolution, 5)))
probeConvList[['ratioB']] <- unique(unlist(lapply(probeConvList$ratioB,kir.allele_resolution, 5)))

## Clean up conflicting copy results
copyDF <- copyDF[apply(copyDF,1,function(x){!any(is.na(x))}),]
copyDF$'KIR2DL23' <- 2

## Matching row and column naming
row.names(copyDF) <- tstrsplit(row.names(copyDF), '_')[[1]]


## ONLY VALIDATION SAMPLES
#validSampleVect <- intersect(rownames(validTypeDF), rownames(haploTypeDF))
#validSampleVect <- intersect(rownames(validTypeDF), rownames(iter.2.typeDF))
#validSampleVect <- intersect(validSampleVect, rownames(filter.2.typeDF))
#validTypeDF <- validTypeDF[validSampleVect,]
#iter.2.typeDF <- iter.2.typeDF[validSampleVect,]
#filter.2.typeDF <- filter.2.typeDF[validSampleVect,]

## Build up a list of sample objects that contain all results
# sample list for ping2_allele_concordance.R
#ping2.sampleList <- valid.ping2.buildSampleList( iter.2.typeDF, filter.2.typeDF, validTypeDF, copyDF ) <- this is now done in ping2_results_read_in.R

#ping2.sampleList <- ping2.sampleList[names(ping2.sampleList) %in% validSampleVect]

#### Check for validation consistency
for(currentSample in ping2.sampleList){
  haploResultBool <- length(currentSample$haploResult) != 0
  pingResultBool <- length(currentSample$pingResult) != 0
  validResultBool <- length(currentSample$validResult) != 0
  copyResultBool <- length(currentSample$copyResult) != 0
  
  ## If type results exist for all 3 datasets
  if(haploResultBool & pingResultBool & validResultBool){
    
    cat('\n',currentSample$name)
    
    ## Transfer the result lists
    haploResultList <- currentSample$haploResult
    pingResultList <- currentSample$pingResult
    validResultList <- currentSample$validResult
    
    ## Identify which loci all datasets share
    sharedLociVect <- intersect(intersect(names(haploResultList), names(pingResultList)), names(validResultList))
    
    ## Iterate through each shared locus
    for(currentLocus in sharedLociVect){
      
      if(currentLocus == 'KIR3DL1' | currentLocus == 'KIR3DS1'){
        if((currentSample$copyResult[['KIR3DL1']] + currentSample$copyResult[['KIR3DS1']]) > 2){
          cat('\nToo many 3DL1/S1 copies')
          next
        }
      }else if(currentSample$copyResult[[currentLocus]] > 2){
        cat('\nToo many',currentLocus,'copies')
        next
      }
      
      ## Pull out the typing table from each dataset
      haploLocusResultDF <- haploResultList[[currentLocus]]
      pingLocusResultDF <- pingResultList[[currentLocus]]
      validLocusResultDF <- validResultList[[currentLocus]]
      
      
      ## Check for any NA results
      haploLocusResultBool <- is.na(haploLocusResultDF)
      pingLocusResultBool <- is.na(pingLocusResultDF)
      validLocusResultBool <- is.na(validLocusResultDF)
      
      ## If there are NA results, skip this result
      if(any(c(haploLocusResultBool, pingLocusResultBool, validLocusResultBool))){
        next
      }
      
      #### Protein lvl typing conversion section
      ## Convert the haplo allele result frame to res3
      haploProtLocusResultDF <- haploLocusResultDF
      haploProtLocusResultDF[] <- lapply(haploResultList[[currentLocus]], function(x){
        tempX <- strsplit(x,'*',fixed=T)
        return(unlist(lapply(tempX, function(y){
          if(y[1] == 'new'){
            return('new')
          }else{
            return(paste0(y[1],'*',substr(y[2],1,3)))
          }
        })))
      })
      haploProtLocusResultDF <- as.data.frame(t(apply(haploProtLocusResultDF,1,sort)),stringsAsFactors = F)
      colnames(haploProtLocusResultDF) <- c('allele1','allele2')
      haploProtLocusResultDF <- unique(haploProtLocusResultDF)
      currentSample$trimHaploResult[[currentLocus]] <- haploProtLocusResultDF
      
      ## Convert the ping allele result frame to res3
      pingProtLocusResultDF <- pingLocusResultDF
      pingProtLocusResultDF[] <- lapply(pingResultList[[currentLocus]], function(x){
        tempX <- strsplit(x,'*',fixed=T)
        return(unlist(lapply(tempX, function(y){
          if(y[1] == 'new'){
            return('new')
          }else{
            return(paste0(y[1],'*',substr(y[2],1,3)))
          }
        })))
      })
      pingProtLocusResultDF <- as.data.frame(t(apply(pingProtLocusResultDF,1,sort)),stringsAsFactors = F)
      colnames(pingProtLocusResultDF) <- c('allele1','allele2')
      pingProtLocusResultDF <- unique(pingProtLocusResultDF)
      currentSample$trimPingResult[[currentLocus]] <- pingProtLocusResultDF
      
      ## Convert the validation allele result frame to res3
      validProtLocusResultDF <- validLocusResultDF
      validProtLocusResultDF[] <- lapply(validResultList[[currentLocus]], function(x){
        tempX <- strsplit(x,'*',fixed=T)
        return(unlist(lapply(tempX, function(y){
          if(y[1] == 'new'){
            return('new')
          }else{
            return(paste0(y[1],'*',substr(y[2],1,3)))
          }
        })))
      })
      #### END SECTION
      validProtLocusResultDF <- as.data.frame(t(apply(validProtLocusResultDF,1,sort)),stringsAsFactors = F)
      colnames(validProtLocusResultDF) <- c('allele1','allele2')
      validProtLocusResultDF <- unique(validProtLocusResultDF)
      
      currentSample$trimValidResult[[currentLocus]] <- validProtLocusResultDF
      
      ###################### 19-11-21 NEW HOPEFULLY FINAL VALIDATION SECTION
      #### Comparing known allele priority integration results to valiattion data
      ## Adding new integration method prioritizing method with not 'new' call
      haploNewBool = F
      if(any(haploProtLocusResultDF$allele1 == 'new')){
        haploNewBool = T
      }
      pingNewBool = F
      if(any(pingProtLocusResultDF$allele1 == 'new')){
        pingNewBool = T
      }
      
      if(haploNewBool & !pingNewBool){
        #### Condition if HAPLO is 'new' allele and PING is known
        ## Set the inter result to be the PING result
        interProtLocusResultDF <- pingProtLocusResultDF
        
      }else if(pingNewBool & !haploNewBool){
        #### Condition if PING is 'new' allele and haplo is known
        ## Set the inter result to be the HAPLO result
        interProtLocusResultDF <- haploProtLocusResultDF
        
      }else{
        #### Condition if both PING and HAPLO are 'new' or known alleles
        ## Intersecting haplo results with ping results
        
        interProtBoolVectFW <- apply(haploProtLocusResultDF, 1, function(x){
          return(any(apply(pingProtLocusResultDF, 1, function(y){
            return(x[2] %in% y[2] & x[1] %in% y[1])
          })))
        })
        
        interProtBoolVectRW <- apply(haploProtLocusResultDF, 1, function(x){
          return(any(apply(pingProtLocusResultDF, 1, function(y){
            return(x[1] %in% y[2] & x[2] %in% y[1])
          })))
        })
        
        interProtBoolVect <- interProtBoolVectFW | interProtBoolVectRW
        
        interProtLocusResultDF <- haploProtLocusResultDF[interProtBoolVect,]
        
        ## PING call bias specific for KIR2DL1
        if( nrow(interProtLocusResultDF) == 0 & currentLocus == 'KIR2DL1' ){
          interProtLocusResultDF <- pingProtLocusResultDF
        }
        
        ## PING call bias specific for KIR2DS4
        if(currentLocus == 'KIR2DS4'){
          interProtLocusResultDF <- pingProtLocusResultDF
        }
        
        ## PING call bias specific for KIR2DS4
        if( currentLocus == 'KIR2DL5' ){
          interProtLocusResultDF <- pingProtLocusResultDF
        }
        
        if( nrow(interProtLocusResultDF) == 0 & currentLocus == 'KIR3DL2'){
          interProtLocusResultDF <- haploProtLocusResultDF
        }
        
        # ## PING call bias specific for KIR3DL2
        # if(currentLocus == 'KIR3DL2'){
        #   
        #   if(is.na(probeRatioDF[currentSample$name,'ratioA'])){
        #     interProtLocusResultDF <- haploProtLocusResultDF[interProtBoolVect,]
        #   }else{
        #   
        #     ratioABool <- probeRatioDF[currentSample$name,'ratioA'] >= 0.15
        #     ratioBBool <- probeRatioDF[currentSample$name,'ratioB'] >= 0.15
        #     
        #     con1 <- ratioABool & !ratioBBool
        #     con2 <- !ratioABool & ratioBBool
        #     con3 <- ratioABool & ratioBBool
        #     
        #     if(con1){
        #       
        #       pingDF <- currentSample$pingResult$KIR3DL2[apply(currentSample$pingResult$KIR3DL2,1,function(x){
        #         all(x %in% probeConvList$ratioA)
        #       }),]
        #       
        #       haploDF <- currentSample$haploResult$KIR3DL2[apply(currentSample$haploResult$KIR3DL2,1,function(x){
        #         all(x %in% probeConvList$ratioA)
        #       }),]
        #       
        #     }else if(con2){
        #       
        #       pingDF <- currentSample$pingResult$KIR3DL2[apply(currentSample$pingResult$KIR3DL2,1,function(x){
        #         all(x %in% probeConvList$ratioB)
        #       }),]
        #       
        #       haploDF <- currentSample$haploResult$KIR3DL2[apply(currentSample$haploResult$KIR3DL2,1,function(x){
        #         all(x %in% probeConvList$ratioB)
        #       }),]
        #       
        #     }else if(con3){
        #       
        #       pingDF <- currentSample$pingResult$KIR3DL2[apply(currentSample$pingResult$KIR3DL2,1,function(x){
        #         con1allele1Match <- x[1] %in% probeConvList$ratioA
        #         con1allele2Match  <- x[2] %in% probeConvList$ratioB
        #         
        #         con2allele1Match <- x[1] %in% probeConvList$ratioB
        #         con2allele2Match <- x[2] %in% probeConvList$ratioA
        #         
        #         if(con1allele1Match & con1allele2Match){
        #           return(TRUE)
        #         }else if(con2allele1Match & con2allele2Match){
        #           return(TRUE)
        #         }else{
        #           return(FALSE)
        #         }
        #       }),]
        #       
        #       haploDF <- currentSample$haploResult$KIR3DL2[apply(currentSample$haploResult$KIR3DL2,1,function(x){
        #         con1allele1Match <- x[1] %in% probeConvList$ratioA
        #         con1allele2Match  <- x[2] %in% probeConvList$ratioB
        #         
        #         con2allele1Match <- x[1] %in% probeConvList$ratioB
        #         con2allele2Match <- x[2] %in% probeConvList$ratioA
        #         
        #         if(con1allele1Match & con1allele2Match){
        #           return(TRUE)
        #         }else if(con2allele1Match & con2allele2Match){
        #           return(TRUE)
        #         }else{
        #           return(FALSE)
        #         }
        #       }),]
        #       
        #     }else{
        #       stop('no 3DL2 probe matches found')
        #       }
        #     
        #     # Call bias (currently haplo biased 8/10/2020)
        #     if(nrow(haploDF) > 0 & nrow(pingDF) > 0){
        #       interProtLocusResultDF <- haploProtLocusResultDF
        #     }else if(nrow(haploDF) > 0){
        #       interProtLocusResultDF <- haploDF
        #     }else if(nrow(pingDF) > 0){
        #       interProtLocusResultDF <- pingDF
        #     }else{
        #       interProtLocusResultDF <- haploProtLocusResultDF
        #     }
        #     
        #     interProtLocusResultDF[] <- lapply(interProtLocusResultDF, function(x){
        #       tempX <- strsplit(x,'*',fixed=T)
        #       return(unlist(lapply(tempX, function(y){
        #         if(y[1] == 'new'){
        #           return('new')
        #         }else{
        #           return(paste0(y[1],'*',substr(y[2],1,3)))
        #         }
        #       })))
        #     })
        #     
        #     interProtLocusResultDF <- as.data.frame(t(apply(interProtLocusResultDF,1,sort)),stringsAsFactors = F)
        #     colnames(interProtLocusResultDF) <- c('allele1','allele2')
        #     interProtLocusResultDF <- unique(interProtLocusResultDF)
        #   }
        # }
        
        ## PING call bias specific for KIR3DL3
        if(nrow(interProtLocusResultDF) == 0 & currentLocus == 'KIR3DL3'){
          interProtLocusResultDF <- pingProtLocusResultDF
        }
        
        ## PING call bias specific for KIR3DL3
        if(nrow(interProtLocusResultDF) == 0 & currentLocus == 'KIR2DL4'){
          interProtLocusResultDF <- haploProtLocusResultDF
        }
        
        ## PING call bias specific for KIR3DL3
        if(nrow(interProtLocusResultDF) == 0 & currentLocus == 'KIR3DL1'){
          interProtLocusResultDF <- haploProtLocusResultDF
        }
        
        ## PING call bias specific for KIR3DL3
        if(nrow(interProtLocusResultDF) == 0 & currentLocus == 'KIR2DL23'){
          interProtLocusResultDF <- haploProtLocusResultDF
        }
        
        ## PING call bias specific for KIR3DL3
        if(nrow(interProtLocusResultDF) == 0 & currentLocus == 'KIR2DS3'){
          interProtLocusResultDF <- haploProtLocusResultDF
        }
        
        ## PING call bias specific for KIR3DL3
        if(nrow(interProtLocusResultDF) == 0 & currentLocus == 'KIR2DS5'){
          interProtLocusResultDF <- haploProtLocusResultDF
        }
        
        ## PING call bias specific for KIR3DL3
        if( currentLocus == 'KIR3DS1' ){
          interProtLocusResultDF <- haploProtLocusResultDF
        }
        
        if( nrow(interProtLocusResultDF) == 0 & currentLocus == 'KIR2DP1' ){
          interProtLocusResultDF <- haploProtLocusResultDF
        }
        
      }
    
      currentSample$interResult[[currentLocus]] <- unique(interProtLocusResultDF)
    }
  }
}


## -- Allele call concordance
lociVect <- unique(unlist(sapply(ping2.sampleList, function(x)names(x[['interResult']]))))
#lociVect <- lociVect[lociVect != 'KIR2DP1']
sampleIDVect <- names(ping2.sampleList)[sapply(ping2.sampleList, function(x)length(x[['interResult']])>0)]
matchIDList <- list()
mismatchIDList <- list()

tempMismatchIDList <- list()

cat('\n\t\tPING\tPING2\tN')
for(currentLocus in lociVect){
  cat('\n',currentLocus)
  
  matchIDList[[currentLocus]] <- list()
  matchI <- 0
  tempMismatchIDList[[currentLocus]] <- list()
  tempMismatchI <- 0
  mismatchIDList[[currentLocus]] <- list()
  mismatchI <- 0
  
  pingMatchCount <- 0
  pingMismatchCount <- 0
  haploMatchCount <- 0
  haploMismatchCount <- 0
  
  for(sampleID in sampleIDVect){
    if(!currentLocus %in% names(ping2.sampleList[[sampleID]]$interResult)){
      next
    }
    
    currentSample <- ping2.sampleList[[sampleID]]
    
    checkDF <- currentSample$interResult[[currentLocus]]
    tempCheckDF <- checkDF
    truthDF <- currentSample$trimValidResult[[currentLocus]]
    
    ## Ignore cases where validation says 'new' allele
    if(any(truthDF == 'new')){
      next
    }
    
    fwCheckBoolVect <- apply(checkDF, 1, function(x){
      return(any(apply(truthDF, 1, function(y){
        return(x[2] %in% y[2] & x[1] %in% y[1])
      })))
    })
    
    rwCheckBoolVect <- apply(checkDF, 1, function(x){
      return(any(apply(truthDF, 1, function(y){
        return(x[1] %in% y[2] & x[2] %in% y[1])
      })))
    })
    
    checkBool <- any(fwCheckBoolVect | rwCheckBoolVect)
    
    if(checkBool){
      matchI <- matchI + 1
      matchIDList[[currentLocus]][[matchI]] <- list('sampleID'=sampleID,
                                                    'haploType'=currentSample$haploResult[[currentLocus]],
                                                    'pingType'=currentSample$pingResult[[currentLocus]],
                                                    'ping2Type'=checkDF,
                                                    'validType'=currentSample$validResult[[currentLocus]])
    }else{
      mismatchI <- mismatchI + 1
      mismatchIDList[[currentLocus]][[mismatchI]] <- list('sampleID'=sampleID,
                                                          'haploType'=currentSample$haploResult[[currentLocus]],
                                                          'pingType'=currentSample$pingResult[[currentLocus]],
                                                          'ping2Type'=checkDF,
                                                          'validType'=currentSample$validResult[[currentLocus]])
    }
    
    checkDF <- currentSample$trimPingResult[[currentLocus]]
    
    fwCheckBoolVect <- apply(checkDF, 1, function(x){
      return(any(apply(truthDF, 1, function(y){
        return(x[2] %in% y[2] & x[1] %in% y[1])
      })))
    })
    
    rwCheckBoolVect <- apply(checkDF, 1, function(x){
      return(any(apply(truthDF, 1, function(y){
        return(x[1] %in% y[2] & x[2] %in% y[1])
      })))
    })
    
    pingCheckBool <- any(fwCheckBoolVect | rwCheckBoolVect)
    
    if(pingCheckBool){
      pingMatchCount <- pingMatchCount + 1
      
      tempMismatchI <- tempMismatchI + 1
      tempMismatchIDList[[currentLocus]][[tempMismatchI]] <- list('sampleID'=sampleID,
                                                                  'haploType'=currentSample$haploResult[[currentLocus]],
                                                                  'pingType'=currentSample$pingResult[[currentLocus]],
                                                                  'ping2Type'=tempCheckDF,
                                                                  'validType'=currentSample$validResult[[currentLocus]])
    }else{
      pingMismatchCount <- pingMismatchCount + 1
    }
    
    checkDF <- currentSample$trimHaploResult[[currentLocus]]
    
    fwCheckBoolVect <- apply(checkDF, 1, function(x){
      return(any(apply(truthDF, 1, function(y){
        return(x[2] %in% y[2] & x[1] %in% y[1])
      })))
    })
    
    rwCheckBoolVect <- apply(checkDF, 1, function(x){
      return(any(apply(truthDF, 1, function(y){
        return(x[1] %in% y[2] & x[2] %in% y[1])
      })))
    })
    
    haploCheckBool <- any(fwCheckBoolVect | rwCheckBoolVect)
    
    if(haploCheckBool){
      haploMatchCount <- haploMatchCount + 1
    }else{
      haploMismatchCount <- haploMismatchCount + 1
    }
  }
  matchStr <- paste(length(matchIDList[[currentLocus]]),'/',(length(matchIDList[[currentLocus]]) + length(mismatchIDList[[currentLocus]])))
  matchPer <- round(length(matchIDList[[currentLocus]])/(length(matchIDList[[currentLocus]]) + length(mismatchIDList[[currentLocus]])),3)
  
  pingMatchPer <- round(pingMatchCount/(pingMatchCount+pingMismatchCount),3)
  
  haploMatchPer <- round(haploMatchCount/(haploMatchCount+haploMismatchCount),3)
  
  nVal <- (length(matchIDList[[currentLocus]]) + length(mismatchIDList[[currentLocus]]))
  cat('\t',pingMatchPer,'\t',haploMatchPer,'\t',matchPer,'\t',nVal,'\t',matchStr)
}

#mismatchIDList[['KIR2DL1']]

#make_unique_pos_frame( kirAlleleDFList[['KIR2DL1']][c('KIR2DL1*010','KIR2DL1*00306'),,drop=F] )

#kirAlleleDFList[['KIR3DL3']][c('KIR3DL3*00102','KIR3DL3*078','KIR3DL3*007'),]
                             
## -- Unresolved genotype frequency

apply(filter.2.typeDF[sampleIDVect,], 2, function(x){
  newSum <- sum( na.omit( grepl('unresolved',x,fixed=T) ) )
  
  knownSum <- sum( na.omit( grepl('*',x,fixed=T) ) )
  

  paste0(newSum,'/',(newSum+knownSum))
})


ping1.unresolvedSumList <- list()
ping1.resolvedSumList <- list()
ping2.unresolvedSumList <- list()
ping2.resolvedSumList <- list()

for(currentLocus in lociVect){
  ping1.unresolvedSumList[[currentLocus]] <- 0
  ping1.resolvedSumList[[currentLocus]] <- 0
  ping2.unresolvedSumList[[currentLocus]] <- 0
  ping2.resolvedSumList[[currentLocus]] <- 0
  for(sampleID in sampleIDVect){
    
    if(!currentLocus %in% names(ping2.sampleList[[sampleID]]$interResult)){
      next
    }
    
    currentSample <- ping2.sampleList[[sampleID]]
    
    checkDF <- currentSample$interResult[[currentLocus]]
    
    ## Ignore cases where validation says 'new' allele
    if(any(checkDF == 'new')){
      ping2.unresolvedSumList[[currentLocus]] <- ping2.unresolvedSumList[[currentLocus]] + 1
    }else{
      ping2.resolvedSumList[[currentLocus]] <- ping2.resolvedSumList[[currentLocus]] + 1
    }
    
    pingType <- ping.1.typeDF[sampleID ,currentLocus]
    
    if( is.na( pingType ) ){
      cat('\nNA found', sampleID, currentLocus)
    }
    
    if( grepl( 'unresolved', pingType, fixed=T) ){
      ping1.unresolvedSumList[[currentLocus]] <- ping1.unresolvedSumList[[currentLocus]] + 1
    }else{
      ping1.resolvedSumList[[currentLocus]] <- ping1.resolvedSumList[[currentLocus]] + 1
    }
    
  }
  ping2.totalN <- (ping2.unresolvedSumList[[currentLocus]] + ping2.resolvedSumList[[currentLocus]])
  ping2.unresN <- ping2.unresolvedSumList[[currentLocus]]
  ping2.unresPer <- round( ping2.unresN / ping2.totalN, 3 )
  
  ping1.totalN <- (ping1.unresolvedSumList[[currentLocus]] + ping1.resolvedSumList[[currentLocus]])
  ping1.unresN <- ping1.unresolvedSumList[[currentLocus]]
  ping1.unresPer <- round( ping1.unresN / ping1.totalN, 3 )
  
  cat('\n',currentLocus,'\t',ping1.unresPer,'\t',ping1.unresN,'\t',ping1.totalN,'\t',ping2.unresPer,'\t',ping2.unresN,'\t',ping2.totalN)
}

