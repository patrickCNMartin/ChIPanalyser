############################################
##### Pre-processing Generic Functions #####
############################################
.vectorSmooth <- function(x, smooth=NULL,norm=FALSE) {
	  originalMax=max(x)
    result=x
	  if(!is.null(smooth)){
        if((smooth %% 2) == 0){smooth = smooth - 1}
        mid = round(smooth/2,0) + 1
        d = smooth - mid
        localSmoothing <- roll_mean(x,smooth,align="center")
				result[mid:(length(x)-d)]<-localSmoothing
    }

    if(norm){result=originalMax*result/max(result)
    }
    return(result)
}

### Doesn't work
.peakParametersExt<-function(peaks){
	  res<-c(median(width(peaks)),sd(width(peaks)))
		return(res)
}


.extractOccupancyDataAtLoci <- function(profile, setSequence, maxSignal=NULL,
    removeBackground=0){
    #compute overlayed occupancy

    # create a list with the occupancy at each loci
    occupancyNorm = vector("list",length(setSequence))
    if(!is.null(setSequence) & length(names(setSequence))>0){
        names(occupancyNorm)=names(setSequence);
    } else{
        names(occupancyNorm)=paste0(seqnames(setSequence),":",
        start(ranges(setSequence)),"..",end(ranges(setSequence)),sep="")
    }

     for(gr in 1:length(setSequence)){
        chr=as.vector(seqnames(setSequence))[gr]
        range=ranges(setSequence)
        startPositon=as.vector(range@start)[gr]
        endPosition = as.vector(end(range))[gr]
        startPos=startPositon
        endPos=endPosition
        positions=startPos:endPos
        occupancyNorm[[gr]]=rep(0,(endPos-startPos+1))

		if(!is.null(profile) & length(profile)>0){
    profile[,2] = as.numeric(profile[,2])
    profile[,3] = as.numeric(profile[,3])
    profile[,4] = as.numeric(profile[,4])
    indexes = intersect(intersect(which(profile[,2]<endPos),
    		which(profile[,3]>=startPos)), which(profile[,1]==chr))
    if(length(indexes)>0){

			  #extract the occupancy
			  for(sector in indexes){
			  minPos = max(profile[sector,2]-startPos,1)
			  maxPos = min(profile[sector,3]-startPos,length(occupancyNorm[[gr]]))
			  occupancyNorm[[gr]][minPos:maxPos] = profile[sector,4]

			  }


			  #normalise the occupancy signal
			  if(is.null(maxSignal)){
					  occupancyNorm[[gr]]=occupancyNorm[[gr]]/max(profile[,4])
			  } else{
					  occupancyNorm[[gr]]=occupancyNorm[[gr]]/maxSignal
			  }

			  #remove the occupancy lower than a threshold
			  occupancyNorm[[gr]][which(occupancyNorm[[gr]] < removeBackground)]=0

			  # smooth the occupancy vector is required
			#  if(!is.null(chipSmooth)){
					#  occupancyNorm[[gr]]=.vectorSmooth(occupancyNorm[[gr]],
						#	  chipSmooth,TRUE)

			  #}
			  }
		}
	  }
	  return(occupancyNorm)
}

.internalLociExtraction<- function(profile, setSequence,reduce=NULL,
	   maxSignal=NULL, removeBackground=0, chipSmooth=NULL,
		 peaks=NULL,Access=NULL,cores=1){


	  ## split for parralle
    subProfile<-.splitDNARanges(setSequence,cores)

    ## parallel
		sub<-parallel::mclapply(subProfile,.internalExtraction,
		profile=profile,maxSignal=maxSignal,
		removeBackground=removeBackground,mc.cores=cores)

    sub<-unlist(sub, recursive=FALSE)

    ## Vector Smoothing if required
    if(!is.null(chipSmooth) & is.null(reduce)){
        for(subset in seq_along(sub)){

            sub[[subset]][which(sub[[subset]]< removeBackground)]<-0
            sub[[subset]]<-.vectorSmooth(sub[[subset]],chipSmooth,TRUE)
        }
				sub<-sub[match(names(setSequence),names(sub))]
    }
    #browser()
    ## If reduce is not NULL only take top hits as defined by reduce
    ## Top average signals for set sequence
    if(!is.null(reduce)){

			  if(is.null(peaks)){

					  if(!is.null(Access)){

                buffer <- which(sapply(lapply(setSequence,function(x,y){
								intersect(x,y)},Access),length)>0)
								sub<-sub[buffer]
								setsub<-setSequence[buffer]
								LocalMean<-sapply(sub,mean)
						    sub<-sub[head(order(LocalMean,decreasing=TRUE),reduce)]
						    setsub<-setSequence[head(order(LocalMean,decreasing=TRUE),
								reduce)]
					  } else {

							  LocalMean <- sapply(sub,mean)
							  sub<-sub[head(order(LocalMean,decreasing=TRUE),reduce)]
							  setsub<-setSequence[head(order(LocalMean,decreasing=TRUE),
								reduce)]
					  }
			  } else {

            if(is.null(Access)){

							  buffer <- sapply(lapply(setSequence,function(x,y){
								intersect(x,y)},peaks),length)
								if(reduce > length(buffer)){
									  reduce <-length(buffer)
								}

								sub<-sub[head(order(buffer,decreasing=TRUE),reduce)]
								setsub<-setSequence[head(order(buffer,decreasing=TRUE),reduce)]
						} else {


								bufferAcces<-which(sapply(lapply(setSequence,function(x,y){
								intersect(x,y)},Access),length)>0)
								setSequence<-setSequence[bufferAcces]
								sub<-sub[bufferAcces]
								bufferPeaks<-sapply(lapply(setSequence,function(x,y){
								intersect(x,y)},peaks),length)
								bufferEnrich<-sapply(sub,function(x){return(max(x))})

								if(reduce > length(bufferPeaks)){
									  reduce <-length(bufferPeaks)
								}

								sub<-sub[head(order(bufferEnrich,
								bufferPeaks,decreasing=TRUE),reduce)]
								setsub<-setSequence[head(order(bufferEnrich,
								bufferPeaks,decreasing=TRUE),reduce)]
						}
			  }

        if(!is.null(chipSmooth)){
            for(subset in seq_along(sub)){
            sub[[subset]][which(sub[[subset]]< removeBackground)]<-0
            sub[[subset]]<-.vectorSmooth(sub[[subset]],chipSmooth,TRUE)
            }
        }
				# Replacing NA by zero (vector smooth arifact)
        if(length(which(sapply(sub, function(x) sum(is.na(x)))>0)) >0){
            for( l in which(sapply(sub, function(x) sum(is.na(x))) >0)){
                sub[[l]]<- rep(0, length(sub[[l]]))
            }
        }


        return(list("ChIPProfle"=sub,"setSequence"=setsub))
    } else {
        # Replacing NA by zero (vector smooth arifact)
        if(length(which(sapply(sub, function(x) sum(is.na(x)))>0)) >0){
        for( l in which(sapply(sub, function(x) sum(is.na(x))) >0)){
            sub[[l]]<- rep(0, length(sub[[l]]))
            }
        }
    return(sub)
    }
}
