
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




.scoreDNAStringSet <- function(PWM, DNASequenceSet, strand="+",
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
# ZeroBackground setting and accessing
.ZeroBackground <- function(object){
    object@ZeroBackground
}

.ZeroBackgroundReplace <- function(object,value){
    object@ZeroBackground <- value
    return(object)
}
# AllSitesAboveThreshold setting
.AllSitesAboveThresholdReplace <- function(object,value){
    object@AllSitesAboveThreshold <- value
    return(object)
}
#NoAccess setting
.NoAccessReplace <- function(object,value){
    object@NoAccess <- value
    return(object)
}
#averageExpPWMScore setting
.averageExpPWMScoreReplace <- function(object,value){
    object@averageExpPWMScore <- value
    return(object)
}
#maxPWMScore setting
.maxPWMScoreReplace <- function(object,value){
    object@maxPWMScore <- value
    return(object)
}
#minPWMScore setting
.minPWMScoreReplace <- function(object,value){
    object@minPWMScore <- value
    return(object)
}

# Validity check function for Genomic Profile Parameters
.is.genomicProfileParameter <- function(object){
    if(class(object) == "genomicProfileParameters"){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Validity check function for occupancy Profile Parameters
.is.occupancyProfileParameters <- function(object){
    if(class(object) == "occupancyProfileParameters"){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Validity check for updated genomic Profile Paramter before
#computing computePWMScore
.is.genomeWideComputed <- function(object){
    if (length(maxPWMScore(object)) > 0 & length(minPWMScore(object)) > 0){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

# Validity check for updated genomic Profile Paramter before
#computing computeOccupancy
.is.PWMScoreComputed <- function(object){
    if(length(AllSitesAboveThreshold(object)) > 0){
        return(TRUE)
    } else {
        return(FALSE)
    }
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
    BPFrequency <- BPFrequency/sum(BPFrequency)
    } else {
    BPFrequency <- (alphabetFrequency(object,baseOnly=TRUE)[,1:4])
    BPFrequency <- BPFrequency/sum(BPFrequency)
    }
    return(BPFrequency)
}
# Compute Thetha from MSE and Corr
.computeTheta <- function(OPP,AccuracyEstimate){
    regionsThreshold <- thetaThreshold(OPP)
    bufferMSE <- unlist(lapply(AccuracyEstimate,function(x){x<-x[[1]][4]}))
    thresholdMSE <- 10^(min(log10(bufferMSE))+regionsThreshold*
        (max(log10(bufferMSE))-min(log10(bufferMSE))))
    bufferCorr <- unlist(lapply(AccuracyEstimate,function(x){x<-x[[1]][3]}))
    thresholdCorrelation <- (min(bufferCorr)+(1-regionsThreshold)*
        (max(bufferCorr)-min(bufferCorr)))

    bufferCorr[ bufferCorr < thresholdCorrelation ] <- thresholdCorrelation
    bufferMSE[ bufferMSE > thresholdMSE ] <- thresholdMSE
    meanTheta <- bufferCorr/bufferMSE

    for(i in 1:length(AccuracyEstimate)){
        for(j in 1:length(AccuracyEstimate[[i]])){
            AccuracyEstimate[[i]][[j]]["meanTheta"] <- meanTheta[[i]]
        }
    }
    return(AccuracyEstimate)
}

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

searchSites <- function(Sites,ScalingFactor="all",
    BoundMolecules="all",Locus="all"){

    #Validity check
    if(class(Sites)!="list" & class(Sites)!="genomicProfileParameters"){
        stop(paste0(deparse(substitute(Sites)),
        " must be a list or a genomicProfileParameters object."))
    }

    if(class(Sites)=="genomicProfileParameters"){
        Sites <- AllSitesAboveThreshold(Sites)
    }

    #Setting Up for search

    lambda <- "lambda = "
    bound <- "boundMolecules = "

    buffer <- c()
    bufferSites <- Sites

    #Search for Scaling Factor
    if(all(ScalingFactor=="all")){
        bufferSites <- bufferSites
        buffer <- 0
    } else {
        localNames <- sapply(strsplit(sapply(strsplit(names(bufferSites)," &"),
        "[[",1),lambda),"[[",2)
        for( i in seq_along(ScalingFactor)){
            buffer <- c(buffer,grep(paste0("^",ScalingFactor[i],"$"),
            localNames))

        }

        if(length(buffer)>1){
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
            buffer <- c(buffer, grep(paste0("^",BoundMolecules[j],"$"),
            localNames))

        }

        if(length(buffer)>1){
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

}
