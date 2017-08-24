profileAccuracyEstimate <- function(LocusProfile,predictedProfile,
    occupancyProfileParameters = NULL){
    #Validity checking
    if(!.is.occupancyProfileParameters(occupancyProfileParameters) &
    !is.null(occupancyProfileParameters)){
    stop(paste0(deparse(substitute(occupancyProfileParameters)),
    "is not an occupancyProfileParamaters Object."))
    }

    if(is.null(occupancyProfileParameters)){
    occupancyProfileParameters <- occupancyProfileParameters()
    }

    #If Sites are not accesible or have no overlapp with chipprofile
    NoAccess<-names(LocusProfile)[(names(LocusProfile) %in%
        names(predictedProfile[[1]])==FALSE)]
    if(length(NoAccess)>0){
    cat("No Profile for:",NoAccess,
    "  --  Do Not Contain Accessible Sites","\n", sep=" ")
    }
    #Reduce ChipProfile to only contain sequences present in predictedProfile
    LocusProfile<-LocusProfile[names(LocusProfile) %in%
        names(predictedProfile[[1]])]

    stepSize<-stepSize(occupancyProfileParameters)
    ProfileSub<-c(TRUE,rep(FALSE,stepSize-1))


    #Calculating Corr, MSE, meanCorr, meanMSE and meanTheta
    AccuracyEstimate<-vector("list", length(predictedProfile))

    for(i in 1:length(AccuracyEstimate)){
        AccuracyEstimate[[i]]<-vector("list", length(predictedProfile[[i]]))
        for(j in 1:length(AccuracyEstimate[[i]])){
            AccuracyEstimate[[i]][[j]]<-as.numeric(as.matrix(
                mcols(predictedProfile[[i]][[j]])))

            correlation<-cor(LocusProfile[[j]][ProfileSub],
                AccuracyEstimate[[i]][[j]])

            MSE<-sum((LocusProfile[[j]][ProfileSub]-
                AccuracyEstimate[[i]][[j]])^2)/
                length(LocusProfile[[j]][ProfileSub])

            AccuracyEstimate[[i]][[j]]<-c("Corr"=correlation,"MSE"=MSE)
        }

        meanCorr<-mean(unlist(lapply(AccuracyEstimate[[i]],"[[",1)),
            na.rm=TRUE)
        meanMSE<-mean(unlist(lapply(AccuracyEstimate[[i]],"[[",2)),
            na.rm=TRUE)*1000
        meanTheta<-meanCorr/meanMSE

        AccuracyEstimate[[i]]<-lapply(AccuracyEstimate[[i]],function(x){
            x<-c(x,"meanCorr"=meanCorr,"meanMSE"=meanMSE,
            "meanTheta"=meanTheta)})
        names(AccuracyEstimate[[i]])<-names(predictedProfile[[i]])
    }

    AccuracyEstimate<-.computeTheta(occupancyProfileParameters,
        AccuracyEstimate)
    names(AccuracyEstimate)<-names(predictedProfile)
    return(AccuracyEstimate)
}
