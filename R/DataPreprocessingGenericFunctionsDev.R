############################################
##### Pre-processing Generic Functions #####
############################################

## checking file format of input data
.formatCheck<-function(profile,type){

    if(type=="ChIP"){
      if(!is.null(colnames(profile))){
          cols<-c(grep("chr",colnames(profile),ignore.case=TRUE),
                grep("start",colnames(profile),ignore.case=TRUE),
                grep("end", colnames(profile),ignore.case=TRUE),
                grep("score", colnames(profile), ignore.case=TRUE))


          profile<-GRanges(seqnames = profile[,cols[1]],
                         ranges = IRanges(profile[,cols[2]],profile[,cols[3]]),
                         score = profile[,cols[4]])
          seqlevelsStyle(profile) <- "UCSC"
      } else{
          if(!all(grepl("chr", profile[,1]))){
            profile[,1]<-paste("chr",profile[,1])
          }
          profile<-GRanges(seqnames=as.character(profile[,1]),
                         ranges= IRanges(profile[,2],profile[,3]),
                         score=profile[,4])
          seqlevelsStyle(profile) <- "UCSC"
    }


  }else{
    if(!is.null(colnames(profile))){
        cols<-c(grep("chr",colnames(profile),ignore.case=TRUE),
              grep("start",colnames(profile),ignore.case=TRUE),
              grep("end", colnames(profile),ignore.case=TRUE))
        if(!all(grepl("chr",cols[1]))){
         cols[1]<-paste("chr",cols[1])
        }

        profile<-GRanges(seqnames = profile[,cols[1]],
                       ranges = IRanges(profile[,cols[2]],profile[,cols[3]]))
    } else{
        if(!all(grepl("chr", profile[,1]))){
          profile[,1]<-paste("chr",profile[,1])
        }
        profile<-GRanges(seqnames=as.character(profile[,1]),
                       ranges= IRanges(profile[,2],profile[,3]))

    }

  }
  return(profile)
}


## profiles smoothing
.vectorSmooth <- function(x, smooth=NULL,norm=FALSE) {

			  ## checking if all values are 0
				## Otherwise it tends to crash
				## RcppRoll doesnt work well with just zeroes

				if(!all(x==0)){
					  # smoothing the vector
					  if(!is.null(smooth)){
							  result<-x
						  	if((smooth %% 2) == 0){
							  		smooth = smooth - 1
						  	}
							  mid = round(smooth/2,0) + 1
							  d = smooth - mid
							  localSmoothing <- roll_mean(x,smooth,align="center")
							  result[mid:(length(x)-d)]<-localSmoothing
					  }
						#normalising the vector
						# norm only works if it NOT all zeroes
						# 0/0 ain't possible fam
						if(norm) {
					      originalMax <- max(x)
				        result <- originalMax*result/max(result)
			      }
				} else{
					  result<-x
				}



    return(result)
}

## extracting peak data
.peakParametersExt<-function(peaks){
	  res<-c(mean(width(peaks)),sd(width(peaks)))
		return(res)
}


## getting sweet profile at lowest level

.scoreReplace <- function(idx,loci, profile){
    loci <-loci[idx]
    score <- rep(0, width(loci))
    overLaps<- profile[queryHits(findOverlaps(profile,loci))]
    overLaps<-overLaps[order(start(overLaps))]


    start <- start(overLaps)-start(loci)
    end <- end(overLaps)-start(loci)

    if(length(overLaps)>0){
    #cutting ends

        if(any(start <1)){
            start[which(start<1)] <-1
        }
        if(end[length(end)]>width(loci)){
            end[length(end)] <- width(loci)
        }


        ranges<-mapply(function(x,y){x:y},start,end,SIMPLIFY=F)
        ChIP<-overLaps$score
        idx<-seq_along(ranges)
        localEnvir<-environment()

        lapply(idx,function(idx){score[ranges[[idx]]]<-ChIP[idx]
                                    assign("score",score,envir=localEnvir)})

    }


    return(score)

}

# sigmoid weighting
.sigmoidWeight<-function(scores,midpoint=1,maxPoint=2, steepness=1){

    scores<-scores*((maxPoint)/(1+exp(-steepness*(scores-midpoint))))
    return(scores)
}


### Getting dem sweet sweet profiles at highest level

.internalChIPExtraction<-function(profile,loci,peaks=NULL,chromatinState=NULL,
                                  reduce=NULL,parameterOptions=NULL,noiseFilter,
                                  cores=cores){


    ## Extracting top values
    maxSignal <- maxSignal(parameterOptions)
    chipSmooth <- chipSmooth(parameterOptions)
    ## first round of reduction
    if(!is.null(peaks) & !is.null(chromatinState)){

        localRanges<-loci[unique(queryHits(findOverlaps(loci,peaks)))]
        localRanges<-localRanges[unique(queryHits(findOverlaps(localRanges,chromatinState)))]

    } else if(!is.null(peaks) & is.null(chromatinState)){

        localRanges<-loci[unique(queryHits(findOverlaps(loci,peaks)))]

    } else if(!is.null(chromatinState) & is.null(peaks)){

        localRanges<-loci[unique(queryHits(findOverlaps(loci,chromatinState)))]

    } else{
        localRanges<-loci
    }

    #split the remaining loci into list just for parallel
    buffer<-.splitRanges(localRanges,cores)
    cores<-buffer$cores
    localRanges<-buffer$rangeSet

    ExtractedProfiles <- unlist(parallel::mclapply(localRanges,FUN=.extractionChIP,
                                   profile=profile,noiseFilter=noiseFilter,
                                   mc.cores=cores),recursive=FALSE)


    if(!is.null(maxSignal) & noiseFilter!="sigmoid"){

        ExtractedProfiles <-parallel::mclapply(ExtractedProfiles,function(x,maxSignal){
                                  x/maxSignal}, maxSignal,mc.cores=cores)
    } else if(!is.null(maxSignal) & noiseFilter =="sigmoid"){
        maxSignal<-max(sapply(ExtractedProfiles,max))
        ExtractedProfiles <-parallel::mclapply(ExtractedProfiles,function(x,maxSignal){
                                  x/maxSignal}, maxSignal,mc.cores=cores)
    }

    ## Reducing if needed
    if(!is.null(reduce)){
        if(reduce > length(ExtractedProfiles)){
           reduce<-length(ExtractedProfiles)
        }

        topScore<-sapply(ExtractedProfiles,max)
        ExtractedProfiles <- ExtractedProfiles[head(order(topScore,decreasing=TRUE),reduce)]

        GenomicRange <-names(ExtractedProfiles)
        #Final setSequence Build
        chr<-sapply(strsplit(GenomicRange,":"),"[[",1)
        start<-sapply(strsplit(sapply(strsplit(GenomicRange,":"),"[[",2),"\\.."),"[[",1)
        end<-sapply(strsplit(sapply(strsplit(GenomicRange,":"),"[[",2),"\\.."),"[[",2)
        loci<-GRanges(seqnames=chr,ranges=IRanges(as.numeric(start), as.numeric(end)))
    }
    # mem drop
    gc()



    # Smoothing dem profiles
    # Smooooooooth operator
    # https://www.youtube.com/watch?v=4TYv2PhG89A

    if(!is.null(chipSmooth)){
       ExtractedProfiles <- parallel::mclapply(ExtractedProfiles,.vectorSmooth,
                                               chipSmooth,TRUE,mc.cores=cores)
    }

    # Replacing NA because vector smooth sucks
    areThereNAs<-which(sapply(ExtractedProfiles,function(x){sum(is.na(x))})>0)
    if(length(areThereNAs) >0){
      for( i in areThereNAs){
          ExtractedProfiles[[i]]<- rep(0, length(ExtractedProfiles[[i]]))
        }
    }



    if(is.null(reduce)){
        return(ExtractedProfiles)
    } else{
        return(list(ExtractedProfiles,loci))
    }


}
