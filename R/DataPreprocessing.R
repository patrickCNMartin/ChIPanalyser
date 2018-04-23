####################################
##### Pre-processing Functions #####
####################################

## ChIP-seq Pre-processing

processingChIPseq <- function(profile,loci=NULL,reduce=NULL,
    occupancyProfileParameters=NULL,peaks=NULL,Access=NULL,cores=1){

    #Validity Check
    if(!is.null(loci) & !is(loci,"GRanges")){
        stop(paste0(deparse(substitute(loci))," must be a GRanges Object"))
    }
    if(!is.null(occupancyProfileParameters) &
        class(occupancyProfileParameters)!="occupancyProfileParameters"){
        stop(paste0(deparse(substitute(occupancyProfileParameters)),
        " must be an occupancyProfileParameters object"))
    }

    # Loading ChIP-seq profile
    if(!is.null(profile) & is(profile,"character")){
        profile <- import(file)
        if(length(grep(x=as.character(seqnames(profile)), pattern="chr"))==0){
            profile <- data.frame("chr"=paste0("chr",
            as.character(seqnames(profile))),"start"=start(profile),
            "end"=end(profile),"score"=profile$score)
        } else {
            profile <- data.frame("chr"=as.character(seqnames(profile)),
            "start"=start(profile), "end"=end(profile),"score"=profile$score)
        }
    } else if(is(profile,"GRanges")){
        if(length(grep(x=as.character(seqnames(profile)), pattern="chr"))==0){
            profile <- data.frame("chr"=paste0("chr",
            as.character(seqnames(profile))),"start"=start(profile),
            "end"=end(profile),"score"=profile$score)
        } else {
            profile <- data.frame("chr"=as.character(seqnames(profile)),
            "start"=start(profile), "end"=end(profile),"score"=profile$score)
        }
    } else {
        profile<-profile
    }

    ## Loading Peaks
    if(!is.null(peaks) & is(peaks,"character")){
        peaks<-import(peaks)
        if(length(grep(x=as.character(seqnames(peaks)), pattern="chr"))==0){
            peaks<-GRanges(seqnames=paste0("chr",as.character(seqnames(peaks))),
                ranges=IRanges(start(peaks),end(peaks)))
        }
    } else if(is(peaks,"GRanges")){
        if(length(grep(x=as.character(seqnames(peaks)), pattern="chr"))==0){
            peaks<-GRanges(seqnames=paste0("chr",as.character(seqnames(peaks))),
                ranges=IRanges(start(peaks),end(peaks)))
        }
    }

    #Setting Sequence if loci is NULL
    if(!is.null(peaks) & is.null(loci)){
        localSeqnames <- unique(profile[,1])
        rangeMatrix <- matrix(0,ncol=2,nrow=length(localSeqnames))
        colnames(rangeMatrix)<-c("start","end")
        for( i in seq_along(localSeqnames)){
            rangeBuffer <- which(profile[,1]==localSeqnames[i])
            rangeMatrix[i,1] <- profile[rangeBuffer[1],2]
            rangeMatrix[i,2] <- profile[rangeBuffer[length(rangeBuffer)],3]
        }
        loci <- GRanges(seqnames=localSeqnames,
            ranges=IRanges(start=rangeMatrix[,"start"],end=rangeMatrix[,"end"]))
        names(loci)<-localSeqnames
    }

    #Occupancy Profile Parameter object
    if(is.null(occupancyProfileParameters)){
        occupancyProfileParameters <- occupancyProfileParameters()
        backgroundSignal(occupancyProfileParameters) <- mean(as.numeric(
            profile[,4]),na.rm=TRUE)
        maxSignal(occupancyProfileParameters) <- max(as.numeric(
            profile[,4]),na.rm=TRUE)
    }




    ChIPProfle<-.internalLociExtraction(profile=profile,
        setSequence=loci,reduce=reduce,
        maxSignal=maxSignal(occupancyProfileParameters),
        removeBackground=removeBackground(occupancyProfileParameters),
        chipSmooth=chipSmooth(occupancyProfileParameters),
        peaks=peaks,Access=Access,cores=cores)
    if(!is.null(reduce)){
        names(ChIPProfle[[2]]) <- names(ChIPProfle[[1]])
    }

    return(list(ChIPProfle,occupancyProfileParameters))

}
