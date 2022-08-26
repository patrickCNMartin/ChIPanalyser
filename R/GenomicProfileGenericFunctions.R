
### Computes the PWM matrix from the PFM matrix
### PFM: the PFM matrix from file
### pseudocount=1: the pseudocount used to compute the PWM
### BPFrequency=rep(0.25,4): a vector with the genome wide
###frequencies of the four nucleotides.
### naturalLog=FALSE: if this is TRUE the method uses the natural logarithm
####when computing the PWM (Berg and von Hippel, 1987),
####otherwise the method uses logarithm in base 2 (Stormo, 2000)
### noOfSites=NULL: the number of DNA sequences used to compute the PFM.
###If NULL it will be computed based on the PFM
### return: a 4 row matrix of the PWM, with rows: A, C, G, T
.computePWM <- function(PFM, pseudocount, BPFrequency,naturalLog,noOfSites){
    if(min(PFM) < 0){
    #the provided PFM is the actual PWM
        PWM <- PFM
    } else{
    #get the sum on each colum (it should be equal)
        if(noOfSites == 0){
            noOfSites <- max(apply(PFM,2,sum))
    }
    #normalise the pfm
    for(i in seq_len(ncol(PFM))){
        PFM[,i] <- PFM[,i]/max(apply(PFM,2,sum))
    }

    #generate the pwm based on the pfm
    PWM <- .computePWMFromPFM(PFM,max(noOfSites), pseudocount,
        BPFrequency,naturalLog)
    }


    return(PWM);
}



### Computed the PWM matrix using the PFM matrix
### PFM: the PFM matrix
### noOfSites: the total number of sites used to compute the PFM matrix
### pseudocount=1: the pseudocount used to compute the PWM
### BPFrequency=rep(0.25,4): a vector with the genome wide
###frequencies of the four nucleotides .
### naturalLog=FALSE: if this is TRUE the method uses the natural
###logarithm  when computing the PWM (Berg and von Hippel, 1987),
###otherwise logarithm  in base 2 is used (Stormo, 2000)
### return: a 4 row matrix of the PWM, with rows: A, C, G, T
.computePWMFromPFM <- function(PFM,noOfSites, pseudocount,
    BPFrequency,naturalLog){
    PWM <- matrix(0,nrow = nrow(PFM),ncol = nrow(PFM))
    PWM <- PFM
    for(i in seq_len(nrow(PFM))){
        for(j in seq_len(ncol(PFM))){
            PWM[i,j] <- (PFM[i,j]*noOfSites+BPFrequency[i]*pseudocount)/
                (noOfSites+pseudocount)
            if(naturalLog){
                PWM[i,j] <- log(PWM[i,j]/BPFrequency[i])
            } else{
                PWM[i,j] <- log2(PWM[i,j]/BPFrequency[i])
            }
        }
    }
    return(PWM)
}

### Attempts to parse a file and extract the PFM. It uses the
###following formats: RAW, TRANSFAC, Jaspar or constructs one from sequences.
###For the correct formats see
###http://www.benoslab.pitt.edu/stamp/help.html#input
### filename: the file storing the motif
### format="raw": the format of the PWM: raw, transfac, jaspar,
###seqs (list of binding sequences)
### return: a PFM matrix, with rows: A, C, G, T
.parsePFM <- function(filename, format){
    PFM <- NULL
    if(format == "raw"){
        PFM <- .parseRawPFM(filename)
    }
    if(format == "transfac"){
        PFM <- .parseTransfacPFM(filename)
    }
    if(format == "JASPAR"){
        PFM <- .parseJASPARPFM(filename)
    }
    if(format == "sequences"){
        PFM <- .parsePFMFromFile(filename)
    }
    #transpose the PFM/PWM if required
    if(nrow(PFM) != 4){
        if(ncol(PFM) == 4){
            PFM <- t(PFM);
        } else{
        stop("Could not load PFM/PWM of 4 rows or columns. Wrong PFM/PWM")
        }
    }
    return(PFM)
}


### Extracts the PFM matrix from a raw PSSM style file.
###It also works for MEME output
### filename: the file containing the motif
### return: the PFM matrix
.parseRawPFM <- function(filename){
    buffer <- read.table(filename,row.names = NULL,stringsAsFactors = FALSE)
    if(ncol(buffer) == 4){
        PFM <- t(buffer)
        rownames(PFM) <- c("A","C","G","T")
    } else{
        PFM <- NULL
    }
    return(PFM)
}


### Extracts the PFM matrix from a transfac style file.
### filename: the file containing the motif
### return: a 4 row PFM matrix, with rows: A, C, G, T
.parseTransfacPFM <- function(filename){
    buffer <- read.table(filename,row.names = NULL,stringsAsFactors = FALSE)
    if(ncol(buffer) > 4){
        PFM <- t(buffer[,2:5])
        if(is.character(PFM[,1])){
        PFM <- PFM[,2:ncol(PFM)]
        }
        rownames(PFM)=c("A","C","G","T")
    } else{
    PFM <- NULL
    }
    return(PFM)
}

### Extracts the PFM matrix from a jaspar style file
### filename: the file containing the motif
### return: a 4 row PFM matrix, with rows: A, C, G, T
.parseJASPARPFM <- function(filename){
    buffer <- readLines(filename)
    error <- FALSE
    for(i in seq_along(buffer)){
        if(!error){
            buffer[i] <- sub("[","",buffer[i],fixed=TRUE)
            buffer[i] <- sub("]","",buffer[i],fixed=TRUE)
            cells <- strsplit(buffer[i]," +",fixed=FALSE)[[1]]

            if(length(cells) > 1){
                values <- as.numeric(cells[2:(length(cells))])
                if(i == 1){
                    PFM <- values
                } else{
                    PFM <- rbind(PFM,values)
                }

            } else{
                error <- TRUE
            }
        }
    }
    if(error){
        PFM <- NULL
    } else {
        rownames(PFM) <- c("A","C","G","T")
    }
    return(PFM)
}




### Reads a list of sequences from a file and produces the PFM mtrix;
### fileame: the filename with the sequences
### return: a 4 row PFM matrix, with rows: A, C, G, T
.parsePFMFromFile <- function(fileame){
    sequences <- readDNAStringSet(fileame)
    if(!is.null(sequences) & length(sequences) > 0){
        PFMExtended <- consensusMatrix(sequences,as.prob = TRUE,
            baseOnly = TRUE)
        PFM <- PFMExtended[1:4,]
    } else{
        PFM <- NULL
    }
    return(PFM)
}

### generates a list of all PWM scores that are above a certain threshold
### PWMScore: a list a PWM scores vectors for each DNA segments
### PWMThreshold: the PWM score threshold
### minPWMScore: the min PWM score
### maxPWMScore: the max PWM score
### return: a list consisting of a vector with the indexes of the PWM
###scores higher than a threshold for each chromosome (or locus)
.getIndexOfPWMThresholded <- function(PWMScore, PWMThreshold){
    indexes <- vector("list",length(PWMScore))
    names(indexes) <- names(PWMScore)
    for(i in seq_along(PWMScore)){
        indexes[[i]] <- which(PWMScore[[i]] > PWMThreshold)
    }
    return(indexes)
}




.scoreDNAStringSet <- function(PWM,DNASequenceSet, strand="+",
    strandRule="max"){
    scoreSet <- vector("list", length(DNASequenceSet))
    IndexPositive <- vector("list",length(DNASequenceSet))
    IndexNegative <- vector("list",length(DNASequenceSet))
    names(scoreSet) <- names(DNASequenceSet)
    for(i in 1:length(DNASequenceSet)){
        sequenceLength <- (length(DNASequenceSet[[i]]) - ncol(PWM) + 1)
        scorePositive <- NULL
        scoreNegative <- NULL

        if(length(grep("+",strand))){
        ## Warnings from BioStrings Package
        scorePositive <- PWMscoreStartingAt(PWM,DNASequenceSet[[i]],
            starting.at = seq_len(sequenceLength))
        }

        if(length(grep("-",strand))){
            scoreNegative <- rev(PWMscoreStartingAt(PWM,
                reverseComplement(DNASequenceSet[[i]]),
                starting.at = seq_len(sequenceLength)))
        }

        if(is.null(scorePositive) & is.null(scoreNegative)){
            # no strand is specified so consider the positive strand
            scoreSet[[i]] <- PWMscoreStartingAt(PWM,DNASequenceSet[[i]],
                starting.at = seq_len(sequenceLength))
        } else if(!is.null(scorePositive) & !is.null(scoreNegative)){
        # both strands are specified
            if(strandRule == "sum"){
                scoreSet[[i]] <- scorePositive + scoreNegative
            } else if(strandRule =="mean"){
                scoreSet[[i]] <- (scorePositive + scoreNegative)/2
            } else{
                scoreSet[[i]] <- pmax(scorePositive,scoreNegative)
                IndexPositive[[i]] <- which(scorePositive ==
                    pmax(scorePositive,scoreNegative))
                IndexNegative[[i]] <- which(scoreNegative ==
                    pmax(scorePositive,scoreNegative))

            }
        } else {
            if(!is.null(scoreNegative)){
                scoreSet[[i]] <- scoreNegative
            } else{
                scoreSet[[i]] <- scorePositive
            }
        }
    }
    return(list(scoreSet,IndexPositive,IndexNegative))
}


# Validity check function for Genomic Profile Parameters
.is.genomicProfiles<- function(object){
    if(class(object) == "genomicProfiles"){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Validity check function for occupancy Profile Parameters
.is.parameterOptions <- function(object){
    if(class(object) == "parameterOptions"){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

.updateGenomicProfiles<-function(genomicProfiles,parameterOptions){
    slots<-slotNames(parameterOptions)
    localEnvir<-environment()
    lapply(slots,function(slots){

           if(!any(slot(genomicProfiles,slots)%in%slot(parameterOptions,slots))){
           slot(genomicProfiles,slots,check=TRUE)<-slot(parameterOptions,slots)
           assign("genomicProfiles",genomicProfiles,envir=localEnvir)
         }})
    return(genomicProfiles)
}

# Extracting and computing Base Pair frequency from BSGenome when
#contructing genomicProfileParamter
.computeBPFrequency <- function(object){
    if(class(object) == "BSgenome"){
    object <- getSeq(object)
    }
    if(class(object) == "DNAStringSet"){
    object <- object
    }
    if(length(object)>1){
    BPFrequency <- apply((alphabetFrequency(object,baseOnly=TRUE)[,1:4]),
    2,sum)
    BPFrequency <- BPFrequency/sum(as.numeric(BPFrequency))
    } else {
    BPFrequency <- (alphabetFrequency(object,baseOnly=TRUE)[,1:4])
    BPFrequency <- BPFrequency/sum(as.numeric(BPFrequency))
    }
    return(BPFrequency)
}
# Compute Thetha from MSE and Corr
#Generate chip Profile
.generateChIPProfileRcpp <- function(inputVector, mean, sd, smooth = NULL,
    norm = TRUE, peakSignificantThreshold=NULL){

    if(is.null(peakSignificantThreshold)){
        peakSignificantThreshold <- 0.0001
    }
    originalMax=max(inputVector);

    var = sd^2;
    shp = mean^2/var
    scl = var/mean
    l = length(inputVector)


    f = dgamma(0:length(inputVector), shape = shp, scale = scl)
    bigF = rev(cumsum(rev(f)))
    peakSignificantSize = min(which(bigF<peakSignificantThreshold))

    kernel = c(rev(bigF[seq(2,peakSignificantSize)]),
    bigF[seq_len(peakSignificantSize)])
    leftKernel=(floor(length(kernel)/2) - 1)
    rightKernel=length(kernel)-leftKernel-1

    rawVector = inputVector;
    rawVector[which(inputVector <= mean(inputVector))]=0
    rawVector = c(rep(0,times=leftKernel), rawVector, rep(0,times=rightKernel))

    peaks = roll_sum(rawVector, length(kernel), weights = kernel)

    if(!is.null(smooth)){
        leftKernel=(floor(smooth/2) - 1)
        rightKernel=smooth-leftKernel-1
        peaks = roll_mean(c(rep(0,times=leftKernel),peaks,
        rep(0,times=rightKernel)), smooth, weights = rep(1, smooth))
    }

    #print((kernel))

    if(norm && max(peaks)!=0){
        peaks=originalMax*peaks/max(peaks)
    }
    peaks <- c(peaks[2:length(peaks)],0)
    return(peaks)
}

#Generate chip Profile
.generateChIPProfile <- function(inputVector, mean, sd,
    smooth = NULL, norm = TRUE, quick = TRUE,peakSignificantThreshold=NULL){

    if(is.null(peakSignificantThreshold)){
        peakSignificantSize <- 1250
    }

    originalMax <- max(inputVector)

    var <- sd^2
    shp <- mean^2/var
    scl <- var/mean
    l <- length(inputVector)



    f <- dgamma(0:length(inputVector), shape = shp, scale = scl)
    gamma <- rev(cumsum(rev(f)))
    sp <- c(rev(gamma[seq(2,peakSignificantSize)]),
    gamma[seq_len(peakSignificantSize)])
    lp <- length(sp)
    hp <- (lp - 1)/2

    peak.centres <- which(inputVector > mean(inputVector))
    peaks <- vector("numeric", l)

    if(!quick){
        for(pc in peak.centres) {
            thisPeak <- vector("numeric", l)
            thisPeak[pc:l] <- gamma[1:(l-pc+1)]
            if(pc!=1){
                thisPeak[1:(pc-1)] <- gamma[pc:2]
            }
            peaks <- peaks + thisPeak * inputVector[pc]
        }
    } else {
        for(pc in peak.centres) {
            thisPeak <- sp
            left <- pc - hp
            right <- pc + hp
            cutLeft <- max(0, 0 - left)
            cutRight <- (-min(0, l - right))
            if(cutLeft != 0){
                thisPeak <- thisPeak[-(1:(cutLeft+1))]
            }
            if(cutRight != 0){
                thisPeak <- thisPeak[-((lp-cutRight+1):lp)]
            }
            left <- max(1, left)
            right <- left + length(thisPeak) - 1
            peaks[left:right] <- peaks[left:right] + thisPeak * inputVector[pc]
        }
    }

    if(!is.null(smooth)){
        if((smooth %% 2) == 0){
            smooth <- smooth - 1
        }
        mid <- round(smooth/2,0) + 1
        d <- smooth - mid
        for(i in mid:(length(peaks) - d)) {
            peaks[i] <- mean(peaks[(i-d):(i+d)])
        }
    }


    if(norm && max(peaks) != 0){
        peaks <- originalMax*peaks/max(peaks)
    }

    return(peaks)
}


## Search sites in occupancy and Chip if a lot of them
## wow that is some good english right there

searchSites <- function(Sites,lambdaPWM="all",
    BoundMolecules="all",Locus="all"){

    if(!is.null(names(Sites)) & all(names(Sites) %in%c("Optimal","Occupancy","ChIPProfiles","goodnessOfFit"))){

       bufferSites<-list("Optimal"=Sites[[1]],
                         "Occupancy"=searchSites(unname(Sites$Occupancy),lambdaPWM=lambdaPWM,
                                                 BoundMolecules=BoundMolecules,
                                                 Locus=Locus),
                         "ChIPProfiles"=searchSites(unname(Sites$ChIPProfile),lambdaPWM=lambdaPWM,
                                                  BoundMolecules=BoundMolecules,
                                                  Locus=Locus),
                          "goodnessOfFit"=searchSites(unname(Sites$goodnessOfFit),lambdaPWM=lambdaPWM,
                                                  BoundMolecules=BoundMolecules,
                                                  Locus=Locus))
       return(bufferSites)
    }else if(class(Sites)=="genomicProfiles" | class(Sites)=="list"){

        if(class(Sites)=="genomicProfiles"){
            Sites <- profiles(Sites)
        }






    #Setting Up for search

    lambda <- "lambda = "
    bound <- "boundMolecules = "

    buffer <- c()
    bufferSites <- Sites

    #Search for Scaling Factor
    if(all(lambdaPWM=="all")){
        bufferSites <- bufferSites
        buffer <- 0
    } else {
        localNames <- sapply(strsplit(sapply(strsplit(names(bufferSites)," &"),
        "[[",1),lambda),"[[",2)
        for( i in seq_along(lambdaPWM)){
            buffer <- c(buffer,grep(paste0("^",lambdaPWM[i],"$"),
            localNames))

        }

        if(length(buffer)>=1){
        bufferSites <- bufferSites[buffer]
        }
    }

    #Search for boundMolecules
    buffer <- c()
    if(all(BoundMolecules=="all")){
        bufferSites <-bufferSites
        buffer <- 0
    } else {
        localNames <- sapply(strsplit(names(bufferSites),bound),"[[",2)
        for( j in seq_along(BoundMolecules)){
            buffer <- c(buffer, grep(paste0("^",as.character(as.integer(BoundMolecules[j])),"$"),
            as.character(as.integer(localNames))))

        }
#browser()
        if(length(buffer)>=1){
            bufferSites <- bufferSites[buffer]
        }
    }
    #Search for Loci

    if(all(Locus == "all")){
        bufferSites <- bufferSites
    } else {
        if(class(bufferSites)=="list"){
        for( k in seq_along(bufferSites)){
            LocalLocus <- names(bufferSites[[k]])
            buffer <- which(LocalLocus==Locus)
            bufferSites[[k]] <- bufferSites[[k]][buffer]
        }
        } else {
            LocalLocus <- names(bufferSites)
            buffer <- which(Locus %in% LocalLocus)
            bufferSites <- bufferSites[buffer]
        }
    }
    return(bufferSites)
  } else{
    stop("Oops Something went wrong. We are not sure what you are parsing to Sites")
  }
}

### Extracting Non Accesible DNA
.AccessExtract<-function(subject,query,isCS=FALSE){
    setLocal<-vector("list",length(subject))

    for(i in seq_along(subject)){
        localIntersect<-setdiff(subject[i], query)
        if(isCS){
          setLocal[[i]]<-data.frame("chr"=as.character(seqnames(localIntersect)),
                                    "start"=start(localIntersect),
                                    "end"=end(localIntersect),
                                    "CS" = mcols(localIntersect)[,1])
        }else {
          setLocal[[i]]<-data.frame("chr"=as.character(seqnames(localIntersect)),
                                    "start"=start(localIntersect),
                                    "end"=end(localIntersect))
        }

    }
    names(setLocal)<-names(subject)
    return(setLocal)
}


#### building a skelton object for accuracy estimates

.cleanGoF <- function(genomicProfiles,ValidationSet,stepSize=10){
     #just in case it needs some reordering

     ValidationSet<-ValidationSet[match(names(genomicProfiles[[1]]),names(ValidationSet))]

     setUpGoF<-lapply(genomicProfiles,function(gp,chip,stepSize){

                      predicted <-lapply(gp,function(x){
                                         return(x$ChIP)})

                      idx<-lapply(gp,function(x){
                                  return(floor(seq(1,sum(width(x)),length.out=length(x))))})
                      validation<-mapply(function(x,idx){
                                  return(x[idx])},x=chip,idx=idx,SIMPLIFY=FALSE)
                      ## internalCheck
                      internalCheck <- all.equal(sapply(predicted,length),sapply(validation,length))
                      if(!internalCheck){
                         stop("predicted and valdiation set are not the same length")
                      }

                      paralleValues<-mapply(function(x,y){
                                     local<-list("prediction"=x,"validation"=y)},
                                     predicted,validation,SIMPLIFY=FALSE)
                      return(paralleValues)
       },ValidationSet,stepSize)

       return(setUpGoF)

}


.meanGoFScore <-function(GoF){
    tags<-names(GoF[[1]])
    idx<-seq_along(tags)

    meanGoF <- sapply(idx,function(idx,GoF){
                      res<-mean(sapply(GoF,"[[",idx),na.rm=TRUE)
                      return(res)
    },GoF=GoF)
    names(meanGoF)<-paste0(tags,"Mean")


    for(i in seq_along(GoF)){
       GoF[[i]]<-c(GoF[[i]],meanGoF)
    }
    return(GoF)
}



## Fscore computation for all regions
.peakExtractionFScore<-function(predicted,locusProfile,region=TRUE){
    # making sure that they are same length just for idx purposes
    if(length(predicted)!=length(locusProfile)){
        stop("precited Profile and locusProfile are not the same length")
    }
    ## threshold increments
    high<-max(max(predicted),max(locusProfile))
    low<-min(min(predicted),min(locusProfile))
    if(low==0){
        subs<-seq(low,high,length.out=21)
        subs<-sapply(subs,"^",2)[seq(2,21)]
    } else{
        subs<-seq(low,high,length.out=20)
        subs<-sapply(subs,"^",2)
    }

    subprof<-matrix(0,ncol=11,nrow=length(subs))

    #rownames(paste0("thresh=",subs))
    #browser()
    matpred<-matrix(0,ncol=length(subs), nrow=length(locusProfile))
    matloc<-matrix(0,ncol=length(subs), nrow=length(locusProfile))

    for(i in seq_along(subs)){

        localProfile<-rep(0,length(locusProfile))
        ## thresh extraction

        localProfile[locusProfile>=subs[i]]<-1
        if(all(localProfile==1)){
             localProfile[1]<-0
        }
        if(all(localProfile==0)){
             localProfile[1]<-1
        }
        matpred[,i]<-factor(predicted)
        matloc[,i]<-factor(localProfile)
    }

    ## Predictions and performance
    pred<-prediction(matpred,matloc)
    prec<-mean(sapply(performance(pred,"prec")@y.values,mean,na.rm=T))
    rec<-mean(sapply(performance(pred,"rec")@y.values,mean,na.rm=T))
    fscore<-mean(sapply(performance(pred,"f")@y.values,mean,na.rm=T))
    acc<-mean(sapply(performance(pred,"acc")@y.values,mean,na.rm=T))
    mcc<-mean(sapply(performance(pred,"mat")@y.values,mean,na.rm=T))

    auc<-mean(sapply(performance(pred,"auc")@y.values,mean,na.rm=T))


    return(c("precision"=prec,"recall"=rec,"f1"=fscore,"accuracy"=acc,"MCC"=mcc,"AUC"=auc))
}

##### geometricRatio

.geometricRatio<-function(predcited,locusProfile,step){
       diff<-abs(locusProfile-predcited)
       leftEdge<-diff[seq(1,length.out=length(diff)-1)]
       rightEdge<-diff[seq(2,length.out=length(diff)-1)]

       areaDiff<-rep(0,length(diff)-1)
       areaShared<-rep(0,length(diff)-1)
       area<-rep(0,length(diff)-1)

       for(i in seq_along(areaDiff)){
       ## Computing Diff and shared
          ## Diff
          areaDiff[i]<-step*((leftEdge[i]+rightEdge[i])/2)
          ## Shared
          if(locusProfile[i]>=predcited[i] & locusProfile[i+1]>=predcited[i+1]){
              areaShared[i]<-step*((predcited[i]+predcited[i+1])/2)
          } else {
              areaShared[i]<-step*((locusProfile[i]+locusProfile[i+1])/2)
          }
        }
        ## Considering zeroes
        if(sum(areaShared)!=0){
            area<-sum(areaDiff)/sum(areaShared)
        } else {

            area<-sum(areaDiff)
        }

    return(area)

}


### all metric extraction just to make things cleaner(should be using switch )

.allMetrics<-function(predicted,locusProfile,stepSize){
    # Checking for bullshit and doing some correlation stuff

    if((sd(predicted,na.rm=TRUE)==0) | (sd(locusProfile,na.rm=TRUE)==0)){
        corrP<-0
        corrS<-0
        corrK<-0
    } else {

        corrP<-cor(na.omit(predicted),na.omit(locusProfile),method="pearson")

        corrS<-cor(na.omit(predicted),na.omit(locusProfile),method="spearman")

        corrK<-cor(na.omit(predicted),na.omit(locusProfile),method="kendall")

    }

    ## ks distance
    ks<-ks.test(predicted,locusProfile)
    D<-unname(ks[[1]])

    ## Geometric
    geo<-.geometricRatio(predicted,locusProfile,stepSize)

    ## fscore
    AUC<-.peakExtractionFScore(predicted,locusProfile)

    metrics<-c("pearson"=corrP,
        "spearman"=corrS,
        "kendall"=corrK,
        "KsDist"=D,
        "geometric"=geo,
         AUC)
    return(metrics)
}



## rankBased method
.OptimalFormatConversion<-function(goodnessOfFit){
     ## for some reason arrays return this weird thing..
     ## things are reordered but not in the way i want them
     ## Im not entirely sure how they work
     ## For now we will use this method
     ## just easier to build
     fitVec<-seq_along(goodnessOfFit[[1]][[1]])
     lociVec<-seq_along(goodnessOfFit[[1]])
     paramVec<-seq_along(goodnessOfFit)

     convertedFormat<-lapply(fitVec,function(fitVec,lociVec,paramVec,goodnessOfFit){


                            param<-lapply(paramVec,function(paramVec,fitVec,lociVec,goodnessOfFit){
                                      loci<-sapply(goodnessOfFit[[paramVec]],"[[",fitVec)
                                      return(loci)

                            },fitVec,lociVec,goodnessOfFit)
                            mat<-do.call("rbind",param)
                            rownames(mat)<-names(goodnessOfFit)
                            return(mat)
     },lociVec,paramVec,goodnessOfFit)
     names(convertedFormat)<-names(goodnessOfFit[[1]][[1]])
     return(convertedFormat)
}


### Processing optimal matricies and parameters


.optimalExtraction <-function(goodnessOfFit,rank=FALSE){

    accu<-profiles(goodnessOfFit)
    bm<-boundMolecules(goodnessOfFit)
    lambda<-lambdaPWM(goodnessOfFit)

    ## converting format
    metrics<-.OptimalFormatConversion(accu)


    if(rank){
         metrics<-metrics[!grepl("Mean",names(metrics))]
         tag<-names(metrics)
         idx<-seq_along(tag)

         # extracting Rank based matrics
         metrics<-lapply(idx,function(idx,metrics,tag){
                              local<-metrics[[idx]]
                              localTag<-tag[idx]
                              internalIdx<-seq_len(ncol(local))
                              TopHits<-lapply(internalIdx,function(idx,met,tags){
                                              if(tags %in% c("geometric","MSE","KsDist")){
                                                 best<-which(met[,idx]==min(met[,idx]))
                                                 buffer<-as.vector(met[best,idx])
                                                 names(buffer)<-rownames(met)[best]


                                              }else{
                                                best<-which(met[,idx]==max(met[,idx]))
                                                buffer<-as.vector(met[best,idx])
                                                names(buffer)<-rownames(met)[best]
                                              }
                                              return(buffer)
                              },local,localTag)
                              #names(TopHits)<-colnames(local)
                              TopHits<-sort(table(names(unlist(TopHits))),decreasing=TRUE)
                              return(TopHits)},metrics,tag)
          ## extracting optimal parameter
          optimalParama<-mapply(function(met,tags){
                              unpacked<-.unpackMetrics(met)
                              if(tags %in% c("geometric","MSE","KsDist")){
                                  optimal<-unpacked[which(unpacked[,"score"]==min(unpacked[,"score"])),c("lambda","boundMolecules")]
                              } else{
                                  optimal<-unpacked[which(unpacked[,"score"]==max(unpacked[,"score"])),c("lambda","boundMolecules")]
                              }
                              return(optimal)

          },met=metrics,tags=tag,SIMPLIFY=FALSE)
          names(optimalParama)<-tag

          # building rank based matricies
          mattemp<-matrix(0,ncol=length(bm),nrow=length(lambda))
          rownames(mattemp)<-lambda
          colnames(mattemp)<-bm
          Loc<-lapply(metrics,function(met,mattemp){
                      unpacked<-.unpackMetrics(met)

                      idx<-seq_len(nrow(unpacked))


                      for(i in idx){
                        mattemp[as.character(unpacked[i,"lambda"]),as.character(unpacked[i,"boundMolecules"])]<-as.numeric(unpacked[i,"score"])
                      }

                      return(mattemp)},mattemp)
          names(Loc)<-tag

          FinalOptimalBuild<-list("OptimalParameters"=optimalParama,"OptimalMatrix"=Loc,"method"="rank")

    } else {
        metrics<-metrics[!grepl("Mean",names(metrics))]
        tag<-names(metrics)
        idx<-seq_along(tag)
        metrics<-lapply(metrics,function(met){
                        return(apply(met,1,mean))
        })

        ## Optimal Paramters
        optimalParama<-mapply(function(met,tag){
                              unpacked<-.unpackMetrics(met)
                              if(tag %in% c("geometric","MSE","KsDist")){
                                  optimal<-unpacked[which(unpacked[,"score"]==min(unpacked[,"score"])),c("lambda","boundMolecules")]
                              } else{
                                  optimal<-unpacked[which(unpacked[,"score"]==max(unpacked[,"score"])),c("lambda","boundMolecules")]
                              }
                              return(optimal) },met=metrics,tag=tag,SIMPLIFY=FALSE)
       ## building matricies
       mattemp<-matrix(0,ncol=length(bm),nrow=length(lambda))
       rownames(mattemp)<-lambda
       colnames(mattemp)<-bm
       Loc<-lapply(metrics,function(met,mattemp){
                   unpacked<-.unpackMetrics(met)

                   idx<-seq_len(nrow(unpacked))


                   for(i in idx){
                      mattemp[as.character(unpacked[i,"lambda"]),as.character(unpacked[i,"boundMolecules"])]<-as.numeric(unpacked[i,"score"])
                   }

                   return(mattemp)},mattemp)
       names(Loc)<-tag

       FinalOptimalBuild<-list("OptimalParameters"=optimalParama,"OptimalMatrix"=Loc,"method"="MeanScore")
    }

   return(FinalOptimalBuild)


}

.unpackMetrics<-function(metrics){

      lambda <- "lambda = "
      bound <- "boundMolecules = "
      lambdaLoc <- sapply(strsplit(sapply(strsplit(names(metrics)," &"),
      "[[",1),lambda),"[[",2)
      bmLoc <- sapply(strsplit(sapply(strsplit(names(metrics)," &"),
      "[[",2),bound),"[[",2)

      unpacked<-cbind("lambda"=as.numeric(lambdaLoc),"boundMolecules"=as.numeric(bmLoc),"score"=as.numeric(metrics))

      return(unpacked)
}




## cleaning function to drop unneccesary objects

.cleanUpAfterYourself<-function(...){
    args<-list(...)
    rm(args)
    gc()
}

### So for some reason unlist
.internalUnlist<-function(object){
    localList<-GRangesList()
    for(i in seq_along(object)){
      localList<-c(localList,object[[i]])
    }
    return(localList)
}



####

.what.is.predictedProfile<-function(predictedProfile){
    if(.is.genomicProfiles(predictedProfile)){
       if(.tags(predictedProfile)=="ChIPProfile"){
         predictedProfile <-profiles(predictedProfile)
    }
  }

    if(class(predictedProfile)=="list" &
              any(grepl("lambda",names(predictedProfile)))){
        stepSize<-unique(width(predictedProfile[[1]][[1]]))
        loci <- unique(unlist(GRangesList(lapply(predictedProfile[[1]], function(x){
                return(GRanges(seqnames=as.character(seqnames(x)),
                               ranges=IRanges(start(x)[1], end(x)[length(x)])))
        }))))

         loci<-rep(loci,length(predictedProfile))

         predictedProfile<-.internalUnlist(predictedProfile)
         predictedProfile<-lapply(predictedProfile, function(x){return(x$ChIP)})

         #warning("No Parameter Combination selected - We will just plot everything",immediate.=TRUE)
    } else if(class(predictedProfile)=="list" &
              any(!grepl("lambda",names(predictedProfile)))){
           stepSize<-unique(width(predictedProfile[[1]]))
           loci <- unique(unlist(GRangesList(lapply(predictedProfile, function(x){
                   return(GRanges(seqnames=as.character(seqnames(x)),
                                  ranges=IRanges(start(x)[1], end(x)[length(x)])))
           }))))
           predictedProfile <- lapply(predictedProfile, function(x){return(x$ChIP)})
    } else if(class(predictedProfile)=="GRanges"){
         stepSize<-unique(width(predictedProfile))
         loci<- GRanges(seqnames=as.character(predictedProfile),
                        ranges=IRanges(start(predictedProfile)[1],end(predictedProfile)[length(predictedProfile)]))
         predictedProfiles<-list(predictedProfile$ChIP)
    } else{
       stop("Oops Somthing went wrong. We are not sure what your are parsing to predictedProfile")
    }
    # building loci from predicted
    return(list("loci"=loci,"predictedProfile"=predictedProfile,"stepSize"=stepSize))
}


.what.is.occupancy<-function(occupancy,PWM=FALSE){
  if(.is.genomicProfiles(occupancy)){
     if(.tags(occupancy)=="Occupancy"){
       occupancy <-profiles(occupancy)
  }
}

    if(class(occupancy)=="list" &
       any(grepl("lambda",names(occupancy)))){
            occupancy<-.internalUnlist(occupancy)

           #warning("No Parameter Combination selected - We will just plot everything",immediate.=TRUE)
      } else if(class(occupancy)=="list" &
                any(!grepl("lambda",names(occupancy)))){


          occupancy-occupancy
      } else if(class(occupancy)=="GRanges"){

           occupancy<-list(occupancy)
      } else{
         stop("Oops Somthing went wrong. We are not sure what your are parsing to predictedProfile")
      }

    return(occupancy)
}



.what.is.ChIPScore <- function(ChIPScore){
    if(class(ChIPScore)=="ChIPScore"){

        ChIPScore<-scores(ChIPScore)


    }

      ## once un packed we can move twords checking the other elements
    if(class(ChIPScore)=="list"){
         ChIPScore<-ChIPScore
    } else if(class(ChIPScore)=="numeric"){
          ChIPScore<-list(ChIPScore)
    }else{
       stop("Oops Somthing went wrong. We are not sure what you are parsing to ChIPScore")
    }
    return(ChIPScore)
}


.what.is.goodnessOfFit<-function(goodnessOfFit){
  if(.is.genomicProfiles(goodnessOfFit)){
     if(.tags(goodnessOfFit)=="GoF"){
       goodnessOfFit <-profiles(goodnessOfFit)
  }
}

    if(class(goodnessOfFit)=="list" &
       any(grepl("lambda",names(goodnessOfFit)))){

                  goodnessOfFit<-unlist(goodnessOfFit,recursive=FALSE)

           #warning("No Parameter Combination selected - We will just plot everything",immediate.=TRUE)
      } else if(class(goodnessOfFit)=="list" &
                any(!grepl("lambda",names(goodnessOfFit)))){
          goodnessOfFit<-goodnessOfFit

      } else if(class(goodnessOfFit)=="numeric") {
           goodnessOfFit <-list(goodnessOfFit)
      } else{
         stop("Oops Somthing went wrong. We are not sure what your are parsing to predictedProfile")
      }
}

.what.is.chromatinState <- function(cs,lociLocal){
    
    if(class(cs) == "GRanges"){
        csList <- list()
        for(i in seq_along(lociLocal)){
            tmp <- cs[queryHits(findOverlaps(cs, lociLocal[i]))]
           
            start(tmp) <- pmax(start(tmp),start(lociLocal[i]))
            end(tmp) <- pmin(end(tmp),end(lociLocal[i]))
            if(ncol(mcols(tmp))==0){
                tmp$stateID <- "Inaccessible DNA"
            } else {
                # Making sure that it is correctly named 
                colnames(mcols(tmp))[1] <- "stateID"
            }
            csList[[i]] <- as.data.frame(tmp)
            
        }
    } else {
        stop("Oopss Somthing went wrong. Please provide GRanges for chromatinState")
    }
    return(csList)
}

.what.is.geneRef <- function(gr,lociLocal){
    if(class(gr) == "GRanges"){
        grList <- list()
        for(i in seq_along(lociLocal)){
           tmp <- gr[queryHits(findOverlaps(gr, lociLocal[i]))]
           start(tmp) <- pmax(start(tmp),start(lociLocal[i]))
           end(tmp) <- pmin(end(tmp),end(lociLocal[i]))
           tmp <- as.data.frame(tmp)
           grList[[i]] <- tmp[,c("seqnames","start","end","width","strand","type","gene_id")]
        }
    } else {
        stop("Oopss Somthing went wrong. Please provide GRanges for chromatinState")
    }
    return(grList)
}
