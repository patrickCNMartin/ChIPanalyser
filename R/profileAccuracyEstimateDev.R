#######################################################################
# Functions profile Accuracy
#######################################################################


profileAccuracyEstimate <- function(genomicProfiles,ChIPScore,
    parameterOptions=NULL,method="all",cores=1){

      #Validity checking
      if(!.is.genomicProfiles(genomicProfiles)){
      stop(paste0(deparse(substitute(genomicProfiles)),
      " is not a Genomic Profiles Object"))
      }
      if(!.is.parameterOptions(parameterOptions) &
      !is.null(parameterOptions)){
      stop(paste0(deparse(substitute(parameterOptions)),
      " is not a parameterOptions Object."))
      }
      if(class(ChIPScore)!="ChIPScore" &
      !is.null(ChIPScore)){
      stop(paste0(deparse(substitute(ChIPScore)),
      " is not a ChIPScore Object."))
      }

      if(!is.null(parameterOptions)){
          genomicProfiles<-.updateGenomicProfiles(genomicProfiles,parameterOptions)
      }

      #If Sites are not accesible or have no overlapp with chipprofile
      dropLoci<-drop(genomicProfiles)
      if(dropLoci!="No loci dropped"){
        widthDisplay<-round(options()$width*0.5)
        cat("No Accessible DNA in: ",paste(rep(" ",
           times=(widthDisplay-nchar("StepSize: ")-nchar(dropLoci[1]))),collapse=''),
           dropLoci,"\n","\n")
      }

     ## Extract and reorder
     predictionSet<-profiles(genomicProfiles)

     ValidationSet<-scores(ChIPScore)

     stepSize <- stepSize(genomicProfiles)

     GoF<-.cleanGoF(predictionSet,ValidationSet,stepSize)
     ## paralle : choosing which one to split over

     if(length(predictionSet) > length(ValidationSet)){
        GoF<-parallel::mclapply(GoF,.GoFPred,method,step=stepSize,mc.cores=cores)
     }else{
        GoF <-lapply(GoF,.GoFLoci,method=method,step=stepSize, cores=cores)
     }

     ### computing mean over all regions
     GoF <- lapply(GoF,.meanGoFScore)

     .profiles(genomicProfiles)<-GoF
     .tags(genomicProfiles)<-"GoF"

     return(genomicProfiles)
}
