####################################
##### Pre-processing Functions #####
####################################

## ChIP-seq Pre-processing


processingChIP <- function(profile,loci=NULL,reduce=NULL,
     peaks=NULL,chromatinState=NULL,parameterOptions=NULL,cores=1){

     #### Validity checking
     ## will do after , first need to know how to do this
     ## also make some god damn test units


     ## Loading ChIP profile
     ## check if RDA files as well

     if(is(profile,"character")){
       if(!grepl(".Rda",profile)){
        fileFormat<-unlist(strsplit(profile,"\\."))
        fileFormat<-fileFormat[length(fileFormat)]
        # Just in case
        # File format are sometimes inconsistent
        # Might need to add more of these,only one i know of atm
        if(fileFormat %in% c("bdg")){
            fileFormat<-"bedGraph"
        }

        profile<-import(profile,format=fileFormat)
        seqlevelsStyle(profile) <- "UCSC"
      } else if(grepl(".Rda",profile)){
        profile<-get(load(profile))
        seqlevelsStyle(profile) <- "UCSC"
      }
     }else if(is(profile, "data.frame")|is(profile,"matrix")){

        profile<-.formatCheck(profile,type="ChIP")


     } else if( is(profile,"GRanges")) {

        profile<-profile
        seqlevelsStyle(profile) <- "UCSC"

     } else {
        stop("Profile: unsuported format")
     }

     # Extracting raw signal metrics
     if(is.null(parameterOptions)){
         parameterOptions<-parameterOptions()
         maxSignal(parameterOptions)<-max(profile$score)
         backgroundSignal(parameterOptions)<-mean(profile$score)

     }else{
         maxSignal(parameterOptions)<-max(profile$score)
         backgroundSignal(parameterOptions)<-mean(profile$score)
     }

     noiseFilter<-noiseFilter(parameterOptions)
    ## Loading peaks
    if(!is.null(peaks)){
        if(is(peaks,"GRanges")){
            peaks<-peaks
            seqlevelsStyle(peaks) <- "UCSC"
        }else if(is(peaks,"character")){
          if(!grepl(".Rda",peaks)){
              fileFormat<-unlist(strsplit(peaks,"\\."))
              fileFormat<-fileFormat[length(fileFormat)]
       # Just in case
       # File format are sometimes inconsistent
       # Might need to add more of these,only one i know of atm
              if(fileFormat %in% c("bdg")){
                  fileFormat<-"bedGraph"
                }

                peaks<-import(peaks,format=fileFormat)
                seqlevelsStyle(peaks) <- "UCSC"
          }else if(grepl(".Rda",peaks)){
              peaks<-get(load(peaks))
              seqlevelsStyle(peaks) <- "UCSC"
        }
      } else if(is(peaks, "data.frame")|is(peaks,"matrix")){
          peaks<-.formatCheck(peaks,type="peaks")
          seqlevelsStyle(peaks) <- "UCSC"
      }else{

          stop("Peaks: unsuported format")

      }
        peakParam<-.peakParametersExt(peaks)
        chipMean(parameterOptions)<-peakParam[1]
        chipSd(parameterOptions)<-peakParam[2]
    }

    ## Loading access
    if(!is.null(chromatinState)){
      if(is(chromatinState,"GRanges")){
        chromatinState<-chromatinState
        seqlevelsStyle(chromatinState) <- "UCSC"
      } else if(is(chromatinState,"character" )){
          if(!grepl(".Rda",chromatinState)){
             fileFormat<-unlist(strsplit(chromatinState,"\\."))
             fileFormat<-fileFormat[length(fileFormat)]
       # Just in case
       # File format are sometimes inconsistent
       # Might need to add more of these,only one i know of atm
             if(fileFormat %in% c("bdg")){
                 fileFormat<-"bedGraph"
             }

             chromatinState<-import(chromatinState,format=fileFormat)
             seqlevelsStyle(peaks) <- "UCSC"
       }else if(grepl(".Rda",chromatinState)){
             chromatinState<-get(load(chromatinState))
             seqlevelsStyle(chromatinState) <- "UCSC"
       }else if(is(chromatinState, "data.frame")| is(chromatinState,"matrix")){
         chromatinState<-.formatCheck(chromatinState,type="peaks")
         seqlevelsStyle(chromatinState) <- "UCSC"
        }
      } else {
          stop("chromatinState: unsuported format")

      }
    }
    ## Creating set loci if none are provided
    ## will need to make change to OPP as well
    if(is.null(loci)){
        if(is.null(parameterOptions)){
            parameterOptions<-parameterOptions()
            lociWidth<-lociWidth(parameterOptions)
        } else {
            lociWidth<-lociWidth(parameterOptions)
        }

        localRanges<-split(profile,seqnames(profile))
        # building range per chromosome
        # tiling is easier
        localRanges<-lapply(localRanges,function(x){
            x<-GRanges(seqnames=unique(as.character(seqnames(x))),
                       ranges=IRanges(min(start(x)),max(end(x))))
            return(x)
        })

        loci<-GRanges()
        for(i in seq_along(localRanges)){
            loci<-suppressWarnings(c(loci, unlist(tile(localRanges[[i]],width=lociWidth))))
        }
        .cleanUpAfterYourself(localRanges)
    }  else if(is(loci,"GRanges")){
       loci <-loci
     } else if(is(loci,"character") & !grepl(".Rda", loci )){
        loci<-import(loci)
        seqlevelsStyle(loci) <- "UCSC"

    } else if(is(loci,"character") & grepl(".Rda",loci)){
        loci <-get(load(loci))
        seqlevelsStyle(loci) <- "UCSC"

    } else{
      stop("Loci: unsuported format.")
    }

    # cores change if required
    if(length(loci)<cores){
      cores<-length(loci)
      warning("Number of cores requested higher than number of loci provided - some cores will be dropped")
    }



    ## Alright lets get this extraction going
    ChIPProfile <- .internalChIPExtraction(profile=profile,
                                          loci=loci,
                                          peaks=peaks,
                                          chromatinState=chromatinState,
                                          reduce=reduce,
                                          parameterOptions=parameterOptions,
                                          noiseFilter=noiseFilter,
                                          cores=cores)



    if(!is.null(reduce)){
        ChIPProfile<-.ChIPScore(scores=ChIPProfile[[1]],loci=ChIPProfile[[2]],
                               maxSignal=max(sapply(ChIPProfile[[1]],max)),
                               backgroundSignal=mean(sapply(ChIPProfile[[1]],mean)),
                               lociWidth=lociWidth(parameterOptions),
                               paramTag="ChIPExt")

    } else{
      ChIPProfile<-.ChIPScore(scores=ChIPProfile,loci=loci,
                             maxSignal=max(sapply(ChIPProfile,max)),
                             backgroundSignal=mean(sapply(ChIPProfile[[1]],mean)),
                             lociWidth=lociWidth(parameterOptions),
                             paramTag="ChIPExt")
    }
    ### crate new Class for custom show method
    ##

    return(ChIPProfile)

}
