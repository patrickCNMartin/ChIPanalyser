############################################
##### Pre-processing Generic Functions #####
############################################
.vectorSmooth <- function(x, smooth=NULL,norm=FALSE) {
    # validity checks for this as it seems to crash a built



			  ## checking if all values are 0
				## Otherwise it tends to crash
				## RcppRoll doesnt work well with just zeroes
				if(all(x!=0)){
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

### Doesn't work
.peakParametersExt<-function(peaks){
	  res<-c(median(width(peaks)),sd(width(peaks)))
		return(res)
}

######
.sigmoidWeight<-function(scores,midpoint=1,maxPoint=2, steepness=1){

    scores<-scores*((maxPoint)/(1+exp(-steepness*(scores-midpoint))))
    return(scores)
}


### Extractiong loci from chip data
.extractOccupancyDataAtLoci <- function(profile, setSequence, maxSignal=NULL,backgroundMethod="zero",peaks=NULL){


	  # create a list with the occupancy at each loci
	  occupancyNorm = vector("list",length(setSequence))
		if(!is.null(setSequence) & length(names(setSequence))>0){
			  name<-names(setSequence)
		    names(occupancyNorm)<-name
	  } else{
				name<-paste0(seqnames(setSequence),":",
				start(ranges(setSequence)),"..",end(ranges(setSequence)),sep="")
		    names(occupancyNorm)<-name
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
			  indexes = intersect(intersect(which(profile[,2]<endPos), which(profile[,3]>=startPos)), which(profile[,1]==chr))
			  if(length(indexes)>0){

			  #extract the occupancy
			  for(sector in indexes){
			  minPos = max(profile[sector,2]-startPos,1)
			  maxPos = min(profile[sector,3]-startPos,length(occupancyNorm[[gr]]))
			  occupancyNorm[[gr]][minPos:maxPos] = profile[sector,4]

			  }


			  #normalise the occupancy signal
			  		if(is.null(maxSignal)){
							  maxSignal<-max(profile[,4])
					  		occupancyNorm[[gr]]=occupancyNorm[[gr]]/maxSignal
			  		} else{
					  		occupancyNorm[[gr]]=occupancyNorm[[gr]]/maxSignal
			  		}
			  }
		}
	  }
		# Noise filtering
		## Method selection

		if(grepl("zero",backgroundMethod,ignore.case=TRUE)){Background<-0}
	  if(grepl("mean",backgroundMethod,ignore.case=TRUE)){
		Background<-mean(sapply(occupancyNorm,mean),na.rm=T)

		}
		if(grepl("median",backgroundMethod,ignore.case=TRUE)){
			Background<-median(sapply(occupancyNorm,median),na.rm=T)

		}
		if(grepl("sigmoid",backgroundMethod,ignore.case=TRUE)){
			#Background<-mean(sapply(occupancyNorm,mean),na.rm=T)
      Background<-0
			#midpoint<-mean(sapply(occupancyNorm,median),na.rm=T)
			midpoint<-mean(sapply(occupancyNorm,quantile,0.95),na.rm=T)
			occupancyNorm<-lapply(occupancyNorm,.sigmoidWeight,midpoint=midpoint)
			NewMax<-max(sapply(occupancyNorm,max))
			occupancyNorm<-lapply(occupancyNorm,"/",NewMax)
		}


		# Filtering background

		occupancyNorm<-lapply(occupancyNorm, function(x){x[x < Background]<-0;return(x)})
    names(occupancyNorm)<-name

		#return occupancy norm
	  return(occupancyNorm)
}

.internalLociExtraction<- function(profile, setSequence,reduce=NULL,occupancyProfileParameters,noiseFilter,peaks=NULL,Access=NULL,cores=1){

    ## Parameters extraction
		maxSignal<-maxSignal(occupancyProfileParameters)
		backgroundMethod<-noiseFilter
		chipSmooth<-chipSmooth(occupancyProfileParameters)
	  ## split for parralle
    subProfile<-.splitDNARanges(setSequence,cores)
    # checking for names
		if(!is.null(setSequence) & is.null(names(setSequence))){
			  names(setSequence)<-paste0(seqnames(setSequence),":",
				start(ranges(setSequence)),"..",end(ranges(setSequence)),sep="")
		}
    ## parallel
		sub<-parallel::mclapply(subProfile,FUN=.internalExtraction,profile=profile,maxSignal=maxSignal,backgroundMethod=backgroundMethod,mc.cores=cores)
    #sub<-.extractOccupancyDataAtLoci(profile, setSequence, maxSignal=maxSignal, removeBackground=removeBackground)
    #browser()
    sub<-unlist(sub, recursive=F)

    ## Vector Smoothing if required
    if(!is.null(chipSmooth) & is.null(reduce)){
        for(subs in seq_along(sub)){
            if(!is.null(sub[[subs]])){
            #sub[[subset]][which(sub[[subset]]< removeBackground)]<-removeBackground
            sub[[subs]]<-.vectorSmooth(sub[[subs]],chipSmooth,TRUE)
          } else{
             next
          }
        }
				seqMatch <- match(names(setSequence),names(sub))
				sub<-sub[seqMatch[which(!is.na(seqMatch))]]
    }

    ## If reduce is not NULL only take top hits as defined by reduce
    ## Top average signals for set sequence
    if(!is.null(reduce)){

			  if(is.null(peaks)){

					  if(!is.null(Access)){
                # something chaned ig Granes in version 3.5...
								# maybe if you had proper unit test this wouldnt happen you piece of shit
								splitSeq<-split(setSequence,names(setSequence))
                buffer <- which(sapply(lapply(splitSeq,function(x,y){intersect(x,y)},Access),length)>0)
								sub<-sub[buffer]
								setsub<-setSequence[buffer]
								LocalMean<-sapply(sub,mean)
						    sub<-sub[head(order(LocalMean,decreasing=T),reduce)]
						    setsub<-setSequence[head(order(LocalMean,decreasing=T),reduce)]
					  } else {

							  LocalMean <- sapply(sub,mean)
							  sub<-sub[head(order(LocalMean,decreasing=T),reduce)]
							  setsub<-setSequence[head(order(LocalMean,decreasing=T),reduce)]
					  }
			  } else {

            if(is.null(Access)){
                splitSeq<-split(setSequence,names(setSequence))
							  buffer <- sapply(lapply(splitSeq,function(x,y){intersect(x,y)},peaks),length)
								if(reduce > length(buffer)){
									  reduce <-length(buffer)
								}

								sub<-sub[head(order(buffer,decreasing=T),reduce)]
								setsub<-setSequence[head(order(buffer,decreasing=T),reduce)]
						} else {

                splitSeq<-split(setSequence,names(setSequence))
							  #buffer <- sapply(lapply(setSequence,function(x,y,z){intersect(intersect(x,y),z)},peaks,Access),length)
								bufferAcces<-which(sapply(lapply(splitSeq,function(x,y){intersect(x,y)},Access),length)>0)
								setSequence<-setSequence[bufferAcces]
								splitSeq<-split(setSequence,names(setSequence))
								sub<-sub[bufferAcces]

								bufferPeaks<-sapply(lapply(splitSeq,function(x,y){intersect(x,y)},peaks),length)
								bufferEnrich<-sapply(sub,max)

								if(reduce > length(bufferPeaks)){
									  reduce <-length(bufferPeaks)
								}

								sub<-sub[head(order(bufferEnrich,bufferPeaks,decreasing=T),reduce)]
								setsub<-setSequence[head(order(bufferEnrich,bufferPeaks,decreasing=T),reduce)]
						}
			  }

        if(!is.null(chipSmooth)){
            for(subs in seq_along(sub)){
                if(!is.null(sub[subs])){
                   sub[[subs]]<-.vectorSmooth(sub[[subs]],chipSmooth,TRUE)
                } else{
                    next
                }
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
