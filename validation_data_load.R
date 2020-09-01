library(data.table,lib.loc = '~/R/x86_64-redhat-linux-gnu-library/3.4/')
library(ggplot2)
library(stringr)
library(methods)
library(plotly,lib.loc='~/R/x86_64-redhat-linux-gnu-library/3.4/')
library(gtools)


#### Validation results read in
validResultsDirectory <- file.path('/home/common_arse/INDIGO/Manually_typed_genos/')
validCenTypePath <- normalizePath(file.path(validResultsDirectory,'INDIGO_controls_first_380_KIRallelesPJN.cen.csv'),mustWork = T)
validTelTypePath <- normalizePath(file.path(validResultsDirectory,'INDIGO_controls_first_380_KIRallelesPJN.tel.csv'),mustWork = T)

validCenTypeDF <- read.csv(validCenTypePath, stringsAsFactors = F, check.names = F,row.names=1)
validCenTypeDF <- validCenTypeDF[,!unlist(lapply(validCenTypeDF, function(x){all(is.na(x))}))]

validTelTypeDF <- read.csv(validTelTypePath, stringsAsFactors = F, check.names = F)
validTelTypeDF <- validTelTypeDF[1:379,]
rownames(validTelTypeDF) <- validTelTypeDF[,1]
validTelTypeDF <- validTelTypeDF[,2:ncol(validTelTypeDF)]

validTelTypeDF <- validTelTypeDF[,!unlist(lapply(validTelTypeDF, function(x){all(is.na(x))}))]

rawValidTypeDF <- merge(validCenTypeDF, validTelTypeDF, by=0, all=TRUE)

### Sanity fixes
rawValidTypeDF['259','2DL23_Allele 1'] <- 'L2*00301/L2*009'
rawValidTypeDF['294','2DL23_Allele 1'] <- 'L2*00301/L2*009'
rawValidTypeDF['316','2DL23_Allele 2'] <- 'L3*00501/L3*017'
rawValidTypeDF['327','2DL23_Allele 1'] <- 'L2*00301/L2*009'
rawValidTypeDF['105',c("2DL5C_allele1","2DL5C_allele2")] <- c('B*00801','B*00601')
rawValidTypeDF['105','2DL5t_allele1'] <- ''
rawValidTypeDF['106',c('3DL3_2')] <- '?'


validColNameInterpretVect <- c('sampleName','KIR3DL3','KIR3DL3','extra','KIR2DS2','KIR2DS2','KIR2DL23','KIR2DL23','KIR2DL5','KIR2DL5','KIR2DS35','KIR2DS35',
                               'KIR2DP1','KIR2DP1','KIR2DL1','KIR2DL1','extra','KIR2DL4','KIR2DL4','KIR2DL4','KIR2DL4','KIR3DL1S1','KIR3DL1S1','KIR3DL1S1','KIR3DL1S1',
                               'KIR2DL5','KIR2DL5','KIR2DS35','KIR2DS35','KIR2DS1','KIR2DS1','KIR2DS4','KIR2DS4','KIR2DS4','KIR3DL2','KIR3DL2','extra')

## Initialize data frame for storing the final valid data
sampleNameVect <- rawValidTypeDF[nchar(rawValidTypeDF[,1]) == 8,1]
validTypeDF <- data.frame(matrix('',nrow=length(sampleNameVect),ncol=14),row.names=sampleNameVect,stringsAsFactors = F)
colnames(validTypeDF) <- c('KIR3DL3','KIR2DS2','KIR2DL23','KIR2DL5','KIR2DS3','KIR2DS5','KIR2DP1',
                           'KIR2DL1','KIR2DL4','KIR3DL1','KIR3DS1','KIR2DS1','KIR2DS4','KIR3DL2')

## Loop through the sample names and build up the final valid dataframe using the raw valid frame
for(sampleNameStr in sampleNameVect){
  rowIndexInt <- which(rawValidTypeDF[,1] == sampleNameStr)
  
  ### Start KIR3DL3 resluts
  genoResultVect <- unlist(rawValidTypeDF[rowIndexInt,c('3DL3_1','3DL3_2')])
  
  con1 <- any(grepl('?',genoResultVect, fixed=T))
  if(con1){
    validTypeDF[sampleNameStr,'KIR3DL3'] <- NA
  }else{
    
    splitGenoList <- strsplit(genoResultVect,'/',fixed=T)
    
    expandGenoDF <- expand.grid(splitGenoList[[1]],splitGenoList[[2]],stringsAsFactors = F)
    
    fullGenoVect <- c()
    for(i in rownames(expandGenoDF)){
      subGenoList <- tstrsplit(unlist(expandGenoDF[i,]),'*',fixed=T)
      subGenoVect <- subGenoList[[length(subGenoList)]]
      subGenoVect <- unlist(lapply(subGenoVect,function(x){
        chopInt <- min(c(nchar(x),5))
        return(substr(x,1,chopInt))
      }))
      
      subGenoVect <- unique(subGenoVect)
      subGenoStr <- paste0('KIR3DL3*',subGenoVect, collapse='+')
      fullGenoVect <- c(fullGenoVect,subGenoStr)
    }
    
    validTypeDF[sampleNameStr,'KIR3DL3'] <- paste0(fullGenoVect,collapse=' ')
  }
  ### Start KIR2DS2 results
  genoResultVect <- unlist(rawValidTypeDF[rowIndexInt,c('2DS2','2DS2.1')])
  
  if(all(lapply(genoResultVect,nchar) == 0)){
    validTypeDF[sampleNameStr,'KIR2DS2'] <- ''
  }else{
    
    if(any(lapply(genoResultVect,nchar) == 0)){
      nullIndexInt <- which(lapply(genoResultVect,nchar) == 0)
      fullIndexInt <- which(lapply(genoResultVect,nchar) != 0)
      
      genoResultVect[nullIndexInt] <- genoResultVect[fullIndexInt]
    }
    splitGenoList <- strsplit(genoResultVect,'/',fixed=T)
    
    expandGenoDF <- expand.grid(splitGenoList[[1]],splitGenoList[[2]],stringsAsFactors = F)
    
    fullGenoVect <- c()
    for(i in rownames(expandGenoDF)){
      subGenoList <- tstrsplit(unlist(expandGenoDF[i,]),'*',fixed=T)
      subGenoVect <- subGenoList[[length(subGenoList)]]
      subGenoVect <- unlist(lapply(subGenoVect,function(x){
        chopInt <- min(c(nchar(x),5))
        return(substr(x,1,chopInt))
      }))
      
      subGenoVect <- unique(subGenoVect)
      subGenoStr <- paste0('KIR2DS2*',subGenoVect, collapse='+')
      fullGenoVect <- c(fullGenoVect,subGenoStr)
    }
    
    validTypeDF[sampleNameStr,'KIR2DS2'] <- paste0(fullGenoVect,collapse=' ')
  }
  
  ### Start KIR2DL23 results
  genoResultVect <- unlist(rawValidTypeDF[rowIndexInt,c('2DL23_Allele 1','2DL23_Allele 2')])
  
  if(any(grepl('R',genoResultVect,fixed=T))){
    validTypeDF[sampleNameStr,'KIR2DL23'] <- 'new'
  }else if(any(nchar(genoResultVect) == 2)){
    validTypeDF[sampleNameStr,'KIR2DL23'] <- 'new'
  }else{
    splitGenoList <- strsplit(genoResultVect,'/',fixed=T)
    
    expandGenoDF <- expand.grid(splitGenoList[[1]],splitGenoList[[2]],stringsAsFactors = F)
    
    fullGenoVect <- c()
    for(i in rownames(expandGenoDF)){
      subGenoList <- tstrsplit(unlist(expandGenoDF[i,]),'*',fixed=T)
      locusNameVect <- subGenoList[[1]]
      subGenoVect <- subGenoList[[length(subGenoList)]]
      subGenoVect <- unlist(lapply(subGenoVect,function(x){
        chopInt <- min(c(nchar(x),5))
        return(substr(x,1,chopInt))
      }))
      
      subGenoVect <- unique(paste0('KIR2D',locusNameVect,'*',subGenoVect))
      subGenoStr <- paste0(subGenoVect,collapse='+')
      fullGenoVect <- c(fullGenoVect,subGenoStr)
    }
    
    validTypeDF[sampleNameStr,'KIR2DL23'] <- paste0(fullGenoVect,collapse=' ')
  }
  
  ### Start KIR2DL5 results
  genoResultVect <- unlist(rawValidTypeDF[rowIndexInt,c("2DL5C_allele1","2DL5C_allele2","2DL5t_allele1","2DL5t_allele2")])
  
  genoResultVect <- genoResultVect[nchar(genoResultVect) > 0]
  
  genoResultVect <- unique(genoResultVect)
  
  if(length(genoResultVect) > 0){
    validTypeDF[sampleNameStr,'KIR2DL5'] <- paste0('KIR2DL5',genoResultVect,collapse='+')
  }else{
    validTypeDF[sampleNameStr,'KIR2DL5'] <- ''
  }
  
  ### Start KIR2DS35 results
  genoResultVect <- rawValidTypeDF[rowIndexInt,c("2DS35  C allele 1","2DS35 C allele 2","2DS35 allele 1","2DS35 allele 2")]
  genoResultVect <- genoResultVect[nchar(genoResultVect) > 0]
  
  genoResultVect <- str_replace(genoResultVect,fixed('_'),'*')
  
  s3GenoResultVect <- unlist(genoResultVect[grepl('S3*',genoResultVect,fixed=T)])
  s5GenoResultVect <- unlist(genoResultVect[grepl('S5*',genoResultVect,fixed=T)])
  
  if(length(s3GenoResultVect) > 0){
    s3GenoResultVect <- unique(s3GenoResultVect)
    
    if(any(grepl('/',s3GenoResultVect,fixed=T))){
      s3GenoResultStr <-  paste0('KIR2D',unlist(strsplit(s3GenoResultVect,'/',fixed=T)), collapse=' ')
    }else{
      s3GenoResultStr <- paste0('KIR2D',s3GenoResultVect, collapse='+')
    }
    
    validTypeDF[sampleNameStr,'KIR2DS3'] <- s3GenoResultStr
  }
  if(length(s5GenoResultVect) > 0){
    s5GenoResultVect <- unique(s5GenoResultVect)
    s5GenoResultStr <- paste0('KIR2D',s5GenoResultVect,collapse='+')
    validTypeDF[sampleNameStr,'KIR2DS5'] <- s5GenoResultStr
  }
  
  ### Start KIR2DP1 results
  genoResultVect <- rawValidTypeDF[rowIndexInt,c("2DP1_Allele 1","2DP1_Allele 2")]
  
  side1 <- unlist(genoResultVect[1], use.names=F)
  side2 <- unlist(genoResultVect[2],use.names=F)
  
  if(side1 == 'neg'){
    side1 <- ''
  }else if(side1 == 'something'){
    side1 <- 'something'
  }else{
    chopInt <- min(c(nchar(side1), 6))
    side1 <- substr(side1,1,chopInt)
    
    if( grepl('*005', side1, fixed=T) ){
      side1 <- '*002'
    }
    
    if( grepl('*016b', side1, fixed=T) ){
      side1 <- '*016'
    }
    
    side1 <- paste0('KIR2DP1',side1)
  }
  if(side2 == 'neg'){
    side2 <- ''
  }else if(side2 == 'something'){
    side2 <- 'something'
  }else{
    side2 <- unlist(strsplit(side2,'/',fixed=T))
    side2List <- lapply(side2,function(x){
      chopInt <- min(c(nchar(x), 6))
      x <- substr(x,1,chopInt)
      
      if( grepl('*005', x, fixed=T) ){
        x <- '*002'
      }
      
      if( grepl('*016b', x, fixed=T) ){
        x <- '*016'
      }
      
      return(x)
    })
    
    side2 <- paste0('KIR2DP1',side2List,collapse=' ')
    
  }
  
  KIR2DP1TypeVect <- unique(c(side1, side2))
  KIR2DP1TypeVect <- KIR2DP1TypeVect[nchar(KIR2DP1TypeVect) > 0]
  
  kir2DP1Type <- paste0(KIR2DP1TypeVect,collapse='+')
  validTypeDF[sampleNameStr,'KIR2DP1'] <- kir2DP1Type
  
  ### Start KIR2DL1 results
  genoResultVect <- rawValidTypeDF[rowIndexInt,c("2DL1_Allele 1","2DL1_Allele 2")]
  
  con1 <- any(grepl('R',genoResultVect,fixed=T))
  con2 <- any(grepl('/',genoResultVect,fixed=T))
  
  if(con1){
    genoResultVect <- 'new'
    validTypeDF[sampleNameStr,'KIR2DL1'] <- genoResultVect
  }else if(con2){
    validTypeDF[sampleNameStr,'KIR2DL1'] <- NA
  }else{
    
    if(any(grepl('neg', genoResultVect,fixed=T))){
      genoResultVect <- genoResultVect[!grepl('neg', genoResultVect,fixed=T)]
    }
    
    genoResultVect <- unlist(lapply(genoResultVect,function(x){
      chopInt <- min(c(nchar(x), 6))
      x <- substr(x,1,chopInt)
      return(x)
    }))
    
    genoResultVect <- unique(genoResultVect)
    genoResultStr <- paste0('KIR2DL1',genoResultVect,collapse='+')
    
    if(length(genoResultVect) == 0){
      validTypeDF[sampleNameStr,'KIR2DL1'] <- ''
    }else{
      validTypeDF[sampleNameStr,'KIR2DL1'] <- genoResultStr
    }
  }
  
  ### Start KIR2DL4 results
  genoResultVect <- rawValidTypeDF[rowIndexInt,c("2DL4 allele 1","allele 2","allele 3","allele 4")]
  
  con1 <- any(grepl('R',genoResultVect,fixed=T))
  con2 <- any(grepl('?', genoResultVect,fixed=T))
  
  if(con1){
    genoResultVect <- 'new'
    validTypeDF[sampleNameStr,'KIR2DL4'] <- genoResultVect
  }else if(con2){
    validTypeDF[sampleNameStr,'KIR2DL4'] <- NA
  }else{
    genoResultVect <- genoResultVect[!grepl('neg', genoResultVect,fixed=T)]
    genoResultVect <- genoResultVect[genoResultVect != '']
    
    genoResultList <- strsplit(genoResultVect,'/',fixed=T)
    
    genoResultChopVect <- unlist(lapply(genoResultList,function(x){
      
      chopInt <- min(c(nchar(x), 6))
      x <- substr(x,1,chopInt)
      return(paste0('KIR2DL4',x,collapse='/'))
    }))
    
    if(any(grepl('/',genoResultChopVect,fixed=T))){
      genoList <- strsplit(genoResultChopVect,'/',fixed=T)
      
      expandGenoDF <- expand.grid(genoList[[1]],genoList[[2]],stringsAsFactors = F)
      
      fullGenoVect <- c()
      for(i in rownames(expandGenoDF)){
        fullGenoVect <- c(fullGenoVect, paste0(expandGenoDF[i,], collapse='+'))
      }
      
      genoResultChopVect <- paste0(fullGenoVect, collapse=' ')
    }
    
    genoResultChopVect <- unique(genoResultChopVect)
    
    validTypeDF[sampleNameStr,'KIR2DL4'] <- paste0(genoResultChopVect,collapse='+')
  }
  
  ### Start KIR3DL1/S1 results
  genoResultVect <- rawValidTypeDF[rowIndexInt,c("3DL1/S1 allele 1","allele 2.1","allele 3.1","allele 4.1")]
  genoResultVect <- genoResultVect[!grepl('neg',genoResultVect,fixed=T)]
  
  genoResultVect <- genoResultVect[genoResultVect != '']
  
  s1Geno <- genoResultVect[grepl('S1',genoResultVect,fixed=T)]
  l1Geno <- genoResultVect[!grepl('S1',genoResultVect,fixed=T)]
  
  con1 <- any(grepl('R',genoResultVect,fixed=T))
  con2 <- any(grepl('?',genoResultVect,fixed=T))
  
  if(con1){
    validTypeDF[sampleNameStr,'KIR3DL1'] <- 'new'
    validTypeDF[sampleNameStr,'KIR3DS1'] <- 'new'
  }else if(con2){
    validTypeDF[sampleNameStr,'KIR3DL1'] <- NA
    validTypeDF[sampleNameStr,'KIR3DS1'] <- NA
  }else{
    if(length(s1Geno) > 0){
      s1Geno <- tstrsplit(s1Geno,'*',fixed=T)[[2]]
      s1ChopGeno <- unlist(lapply(s1Geno,function(x){
        chopInt <- min(c(nchar(x), 5))
        x <- substr(x,1,chopInt)
        return(x)
      }))
      
      s1ChopGeno <- unique(s1ChopGeno)
      validTypeDF[sampleNameStr,'KIR3DS1'] <- paste0('KIR3DS1*',s1ChopGeno,collapse='+')
    }
    
    if(length(l1Geno) > 0){
      l1Geno <- tstrsplit(l1Geno,'*',fixed=T)[[2]]
      
      l1ChopGeno <- unlist(lapply(l1Geno,function(x){
        chopInt <- min(c(nchar(x), 5))
        x <- substr(x,1,chopInt)
        return(x)
      }))
      
      l1ChopGeno <- unique(l1ChopGeno)
      validTypeDF[sampleNameStr,'KIR3DL1'] <- paste0('KIR3DL1*',l1ChopGeno,collapse='+')
    }
  }
  
  ### Start KIR2DS1 results
  genoResultVect <- rawValidTypeDF[rowIndexInt,c("2DS1_allele1","2DS1_allele2")]
  genoResultVect <- genoResultVect[!grepl('neg',genoResultVect,fixed=T)]
  
  genoResultVect <- genoResultVect[genoResultVect != '']
  
  if(length(genoResultVect) > 0){
    genoResultVect <- unique(genoResultVect)
    validTypeDF[sampleNameStr,'KIR2DS1'] <- paste0('KIR2DS1',genoResultVect,collapse='+')
  }else{
    validTypeDF[sampleNameStr,'KIR2DS1'] <- ''
  }
  
  ### Start KIR2DS4 results
  genoResultVect <- rawValidTypeDF[rowIndexInt,c("2DS4_Allele 1","2DS4_Allele 2","2DS4_Allele 3")]
  genoResultVect <- genoResultVect[!grepl('neg',genoResultVect,fixed=T)]
  
  genoResultVect <- genoResultVect[genoResultVect != '']
  
  genoResultVect <- str_replace_all(genoResultVect,'KIR2DS4_',fixed('*'))
  
  if(length(genoResultVect) > 0){
    genoResultVect <- unique(genoResultVect)
    
    if(any(grepl('/',genoResultVect,fixed=T))){
      genoList <- strsplit(genoResultVect,'/',fixed=T)[[1]]
      
      genoResultVect <- paste0('KIR2DS4',genoList, collapse=' ')
      validTypeDF[sampleNameStr,'KIR2DS4'] <- genoResultVect
    }else{
      validTypeDF[sampleNameStr,'KIR2DS4'] <- paste0('KIR2DS4',genoResultVect,collapse='+')
    }
  }else{
    validTypeDF[sampleNameStr,'KIR2DS4'] <- ''
  }
  
  ### Start KIR3DL2 results
  genoResultVect <- rawValidTypeDF[rowIndexInt,c("3DL2 allele 1","3DL2 allele 2")]
  genoResultVect <- genoResultVect[!grepl('neg',genoResultVect,fixed=T)]
  
  genoResultVect <- genoResultVect[genoResultVect != '']
  
  con1 <- any(grepl('R',genoResultVect,fixed=T))
  con2 <- any(grepl('g',genoResultVect,fixed=T))
  con3 <- any(grepl('/',genoResultVect,fixed=T))
  
  if(con1){
    validTypeDF[sampleNameStr,'KIR3DL2'] <- 'new'
  }else{
    
    if(con2){
      genoResultVect <- str_replace_all(genoResultVect,'g..','')
    }
    
    if(con3){
      genoList <- strsplit(genoResultVect,'/',fixed=T)
      
      expandGenoDF <- expand.grid(genoList[[1]],genoList[[2]],stringsAsFactors = F)
      
      fullGenoVect <- c()
      for(i in rownames(expandGenoDF)){
        currentGeno <- expandGenoDF[i,]
        currentGeno <- unlist(lapply(currentGeno,function(x){
          chopInt <- min(c(nchar(x), 6))
          x <- substr(x,1,chopInt)
          return(x)
        }))
        fullGenoVect <- c(fullGenoVect, paste0('KIR3DL2',currentGeno, collapse='+'))
      }
      
      fullGenoVect <- unique(fullGenoVect)
      
      genoResultVect <- paste0(fullGenoVect, collapse=' ')
      validTypeDF[sampleNameStr,'KIR3DL2'] <- genoResultVect
    }else{
      fullGenoVect <- unlist(lapply(genoResultVect,function(x){
        chopInt <- min(c(nchar(x), 6))
        x <- substr(x,1,chopInt)
        return(x)
      }))
      
      fullGenoVect <- unique(fullGenoVect)
      validTypeDF[sampleNameStr,'KIR3DL2'] <- paste0('KIR3DL2',fullGenoVect, collapse='+')
    }
  }
}
## Some corrections
validTypeDF['IND00283',c('KIR2DS3','KIR2DS5')] <- 'new'

validTypeDF <- validTypeDF[rownames(validTypeDF) != 'IND00026',]
validTypeDF <- validTypeDF[rownames(validTypeDF) != 'IND00033',]
validTypeDF <- validTypeDF[rownames(validTypeDF) != 'IND00086',]

new3DL3Vect <- c('IND00006',
                 'IND00050',
                 'IND00088',
                 'IND00099',
                 'IND00142',
                 'IND00143',
                 'IND00167',
                 'IND00172',
                 'IND00179',
                 'IND00189',
                 'IND00201',
                 'IND00261',
                 'IND00279',
                 'IND00328',
                 'IND00340',
                 'IND00389')

validTypeDF[new3DL3Vect,'KIR3DL3'] <- 'new'

new2DL1Vect <- c('IND00143','IND00167','IND00279')
validTypeDF[new2DL1Vect,'KIR2DL1'] <- 'new'

new2DP1Vect <- c('IND00067', 'IND00159', 'IND00166', 'IND00167', 'IND00185', 'IND00345')
validTypeDF[new2DP1Vect,'KIR2DP1'] <- 'new'

new3DL1Vect <- c('IND00062', 'IND00143', 'IND00156', 'IND00178', 'IND00357')
validTypeDF[new3DL1Vect,'KIR3DL1'] <- 'new'

new3DL2Vect <- c('IND00021', 'IND00148', 'IND00178', 'IND00203', 'IND00204', 'IND00226', 'IND00253', 'IND00257', 'IND00270'
                 , 'IND00284', 'IND00321', 'IND00322', 'IND00353', 'IND00356', 'IND00363', 'IND00379', 'IND00390')
validTypeDF[new3DL2Vect,'KIR3DL2'] <- 'new'


new2DL23Vect <- c('IND00039', 'IND00042', 'IND00060', 'IND00118', 'IND00134', 'IND00143', 'IND00154', 'IND00161'
                  ,'IND00166', 'IND00167', 'IND00178', 'IND00205', 'IND00211', 'IND00255', 'IND00258')
validTypeDF[new2DL23Vect,'KIR2DL23'] <- 'new'

