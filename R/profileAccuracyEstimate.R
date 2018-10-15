profileAccuracyEstimate <- function(LocusProfile,predictedProfile,
    occupancyProfileParameters = NULL,method="all"){

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
    ## Note this loops over all regions selected
    ## Loop over paramter combination
    for(i in 1:length(AccuracyEstimate)){
        AccuracyEstimate[[i]]<-vector("list", length(predictedProfile[[i]]))


        ## Loop over selected Loci

        for(j in 1:length(AccuracyEstimate[[i]])){
            ## Extracting predicted scores
            AccuracyEstimate[[i]][[j]]<-as.numeric(as.matrix(
                mcols(predictedProfile[[i]][[j]])))

            ## MSE betwenn the two curves
            MSE<-sum((LocusProfile[[j]][ProfileSub]-
                AccuracyEstimate[[i]][[j]])^2)/
                length(LocusProfile[[j]][ProfileSub])
            ## Pearson correlation Method
            if(grepl("pearson",method,ignore.case=T)){
                if((sd(LocusProfile[[j]][ProfileSub],na.rm=TRUE)==0) |
                    (sd(AccuracyEstimate[[i]][[j]],na.rm=TRUE)==0)){
                    correlation <-0
                } else {
                    correlation<-cor(na.omit(LocusProfile[[j]][ProfileSub]),
                    na.omit(AccuracyEstimate[[i]][[j]]),method="pearson")
                }
                  AccuracyEstimate[[i]][[j]]<-c("pearson"=correlation,"MSE"=MSE)
            } else if(grepl("spearman",method,ignore.case=T)){
                if((sd(LocusProfile[[j]][ProfileSub],na.rm=TRUE)==0) |
                    (sd(AccuracyEstimate[[i]][[j]],na.rm=TRUE)==0)){
                    correlation <-0
                } else {
                    correlation<-cor(na.omit(LocusProfile[[j]][ProfileSub]),
                    na.omit(AccuracyEstimate[[i]][[j]]),method="spearman")
                }
                AccuracyEstimate[[i]][[j]]<-c("spearman"=correlation,"MSE"=MSE)
            }else if(grepl("kendall",method,ignore.case=T)){
                if((sd(LocusProfile[[j]][ProfileSub],na.rm=TRUE)==0) |
                    (sd(AccuracyEstimate[[i]][[j]],na.rm=TRUE)==0)){
                    correlation <-0
                } else {
                    correlation<-cor(na.omit(LocusProfile[[j]][ProfileSub]),
                    na.omit(AccuracyEstimate[[i]][[j]]),method="kendall")
              }
              AccuracyEstimate[[i]][[j]]<-c("kendall"=correlation,"MSE"=MSE)
            }else if(grepl("ks",method,ignore.case=T)){
                ks<-ks.test(LocusProfile[[j]][ProfileSub],AccuracyEstimate[[i]][[j]])
                D<-ks[[1]]
                pval<-ks[[2]]
                AccuracyEstimate[[i]][[j]]<-c("MSE"=MSE,"KsDist"=ks[[1]],"KsPval"=ks[[2]])
            }else if(grepl("fscore",method,ignore.case=T)){
                AccuracyEstimate[[i]][[j]]<-.peakExtractionFScore(AccuracyEstimate[[i]][[j]],LocusProfile[[j]][ProfileSub])

            }else if(grepl("geometric",method,ignore.case=T)){
                AccuracyEstimate[[i]][[j]]<-c("MSE"=MSE,"geometric"=.geometricRatio(AccuracyEstimate[[i]][[j]],LocusProfile[[j]][ProfileSub],step=stepSize))
            }else if(grepl("all",method,ignore.case=T)){
                allMetrics<-.allMetrics(AccuracyEstimate[[i]][[j]],LocusProfile[[j]][ProfileSub],stepSize=stepSize)
                AccuracyEstimate[[i]][[j]]<-list(c("MSE"=MSE,allMetrics[[1]]),allMetrics[[2]])
            }

        }


        if(any(method %in% c("pearson","spearman","kendall"))){
        meanCorr<-mean(sapply(AccuracyEstimate[[i]],"[[",1), na.rm=TRUE)
        meanMSE<-mean(sapply(AccuracyEstimate[[i]],"[[",2), na.rm=TRUE)
        AccuracyEstimate[[i]]<-lapply(AccuracyEstimate[[i]],function(x){
                x<-c(x,"meanCorr"=meanCorr,"MSEMean"=meanMSE)})
        }else if(any(method=="ks")){
            meanDstat<-mean(sapply(AccuracyEstimate[[i]],"[[",2),na.rm=TRUE)
            meanMSE<-mean(sapply(AccuracyEstimate[[i]],"[[",1),na.rm=TRUE)
            AccuracyEstimate[[i]]<-lapply(AccuracyEstimate[[i]],function(x){
                x<-c(x,"meanDstat"=meanDstat,"MSEMean"=meanMSE)})
        }else if(any(method=="fscore")){
            precisionMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["precision"]}),na.rm=T)
            recallMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["recall"]}),na.rm=T)
            FscoreMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["f1"]}),na.rm=T)
            MCCMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["MCC"]}),na.rm=T)
            AccuracyMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["accuracy"]}),na.rm=T)
            AucMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["AUC"]}),na.rm=T)

            AccuracyEstimate[[i]]<-lapply(AccuracyEstimate[[i]],function(x){

                res<-c("precisionMean"=precisionMean,
                         "recallMean"=recallMean,
                         "FscoreMean"=FscoreMean,
                         "AccuracyMean"=AccuracyMean,
                         "MCCMean"=MCCMean,
                         "AUCMean"=AucMean,x[[1]])
                return(list(res,x[[2]]))
            })
        }else if(any(method=="geometric")){
            meanGeoRatio<-mean(sapply(AccuracyEstimate[[i]],"[[",2), na.rm=TRUE)
            meanMSE<-mean(sapply(AccuracyEstimate[[i]],"[[",1),na.rm=TRUE)
            AccuracyEstimate[[i]]<-lapply(AccuracyEstimate[[i]],function(x){
                    x<-c(x,"meanGeoRatio"=meanGeoRatio,"MSEMean"=meanMSE)})
        }else if(any(method=="all")){
          precisionMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["precision"]}),na.rm=T)
          recallMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["recall"]}),na.rm=T)
          FscoreMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["f1"]}),na.rm=T)
          MCCMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["MCC"]}),na.rm=T)
          AccuracyMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["accuracy"]}),na.rm=T)
          AucMean<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["AUC"]}),na.rm=T)
          meanPearson<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["pearson"]}),na.rm=T)
          meanSpearman<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["spearman"]}),na.rm=T)
          meanKendall<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["kendall"]}),na.rm=T)
          meanksdist<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["KsDist"]}),na.rm=T)
          meanGeoRatio<-mean(sapply(AccuracyEstimate[[i]],function(x){x[[1]]["geometric"]}),na.rm=T)
          meanMSE<-mean(sapply(AccuracyEstimate[[i]],"[[",1),na.rm=TRUE)
          AccuracyEstimate[[i]]<-lapply(AccuracyEstimate[[i]],function(x){

              res<-c(x[[1]][1:7],
                       "pearsonMean"=meanPearson,
                       "spearmanMean"=meanSpearman,
                       "kendallMean"=meanKendall,
                       "MSEMean"=meanMSE,
                       "ksMean"=meanksdist,
                       "geometricMean"=meanGeoRatio,
                       "precisionMean"=precisionMean,
                       "recallMean"=recallMean,
                       "FscoreMean"=FscoreMean,
                       "AccuracyMean"=AccuracyMean,
                       "MCCMean"=MCCMean,
                       "AUCMean"=AucMean,x[[1]][8:length(x[[1]])])
              return(list(res,x[[2]]))
          })

        }


        names(AccuracyEstimate[[i]])<-names(predictedProfile[[i]])
    }

    #AccuracyEstimate<-.computeTheta(occupancyProfileParameters,
        #AccuracyEstimate)
    names(AccuracyEstimate)<-names(predictedProfile)
    return(AccuracyEstimate)
}
