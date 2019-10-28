##############################################
#############      S4 Methods    #############
##############################################

########### Show Methods ##########
setMethod("show",
    signature = "genomicProfiles",
    definition = function(object){
        widthDisplay<-round(options()$width*0.75)
        if(object@tags == "empty"){
          cat(rep("_", widthDisplay),"\n","\n")
          title <- as.character("genomicProfiles object with no internal data")
          titleWidth<-nchar(title, type="byte")
          cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
          cat(rep("_", widthDisplay),"\n")
          slots <- slotNames(object)[1:4]

          line1<-paste("Suggested Slots : ",paste(rep(" ",
            times=(widthDisplay-nchar("Suggested Slots : ")-nchar(slots[1]))),collapse=''),
            slots[1],"\n")
          line2<-paste("Suggested Slots : ",paste(rep(" ",
            times=(widthDisplay-nchar("Suggested Slots : ")-nchar(slots[2]))),collapse=''),
            slots[2],"\n")
          line3<-paste("Suggested Slots : ",paste(rep(" ",
            times=(widthDisplay-nchar("Suggested Slots : ")-nchar(slots[3]))),collapse=''),
            slots[3],"\n")
          line4<-paste("Suggested Slots : ",paste(rep(" ",
              times=(widthDisplay-nchar("Suggested Slots : ")-nchar(slots[4]))),collapse=''),
              slots[4],"\n")
          cat("",line1,line2,line3,line4)
        }else  if(object@tags == "PWM_stage"){
          cat(rep("_", widthDisplay),"\n","\n")
          title <- as.character("genomicProfiles object: Position Weight Matrix")
          titleWidth<-nchar(title, type="byte")
          cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
          cat(rep("_", widthDisplay),"\n","\n")
          bpf <- object@BPFrequency
          names(bpf)<-c("A" ,"C" ,"G" ,"T")

          cat("Position Weight Matrix (PWM): \n","\n")
          print(object@PWM)
          cat(rep("_", widthDisplay),"\n","\n")
          cat("Base pairs frequency for PWM weighting : \n","\n")
          print(bpf)
          cat(rep("_", widthDisplay),"\n","\n")
          cat("PWM built from (PFM - Format ",object@PFMFormat, ") : \n","\n")
          print(object@PFM)
          cat("\n",rep("_", widthDisplay),"\n")

        }else  if(grepl("genomeWide",object@tags)){
          cat(rep("_", widthDisplay),"\n","\n")
          if(grepl("CS",object@tags)){
              title <- as.character("genomicProfiles object: Accessible Chromatin Score")
          } else{
              title <- as.character("genomicProfiles object: Genome Wide Score")
          }
          titleWidth<-nchar(title, type="byte")
          cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
          cat(rep("_", widthDisplay),"\n","\n")

          cat("","Position Weight Matrix (PWM): \n","\n")
          print(object@PWM)
          cat(rep("_", widthDisplay),"\n","\n")
          cat("",paste0("maxPWMScore: ",paste(rep(" ",
            times=(widthDisplay-nchar("maxPWMScore: ")-nchar(object@maxPWMScore))),collapse=''),
            object@maxPWMScore,"\n"))

          cat("",paste0("minPWMScore: ",paste(rep(" ",
            times=(widthDisplay-nchar("minPWMScore: ")-nchar(object@minPWMScore))),collapse=''),
            object@minPWMScore,"\n"))

          cat("",paste0("averageExpPWMScore: ",paste(rep(" ",
            times=(widthDisplay-nchar("averageExpPWMScore: ")-nchar(paste0(object@averageExpPWMScore[1],"(lambda = ",object@lambdaPWM[1],")")))),collapse=''),
            object@averageExpPWMScore," (lambda = ",object@lambdaPWM,")\n"))

          cat("",paste0("DNASequenceLength: ",paste(rep(" ",
             times=(widthDisplay-nchar("DNASequenceLength: ")-nchar(object@DNASequenceLength))),collapse=''),
             object@DNASequenceLength,"\n"))
         cat(rep("_", widthDisplay),"\n","\n")

       }else if(grepl("PWMScore",object@tags)){
         cat(rep("_", widthDisplay),"\n","\n")
         title <- as.character("genomicProfiles object: PWM scores above Threshold")
         titleWidth<-nchar(title, type="byte")
         cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
         cat(rep("_", widthDisplay),"\n","\n")
         profile<-object@profiles

         if(length(profile)<2){
             cat("","Scores above Threshold in locus: ","\n")
             print(profile)
         }else{
             cat("","Scores above Threshold in locus: ","\n")
             print(profile[[1]])
             cat("\n")
             cat("",length(profile)-1," remaining loci","\n")
             cat("",paste0("Locus: ",paste(rep(" ",
               times=(widthDisplay-nchar("Locus: ")-nchar(paste0(names(profile)[1],"(length = ",length(profile[[1]]),")")))),collapse=''),
               names(profile)[2:length(profile)]," (length = ",sapply(profile, length)[2:length(profile)],")\n"))

         }
         cat(rep("_", widthDisplay),"\n","\n")
         title <- as.character("Associated parameters")
         titleWidth<-nchar(title, type="byte")
         cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
         cat(rep("_", widthDisplay),"\n","\n")
         cat("",paste0("minPWMScore: ",paste(rep(" ",
           times=(widthDisplay-nchar("minPWMScore: ")-nchar(object@minPWMScore))),collapse=''),
           object@minPWMScore,"\n"))
         cat("",paste0("maxPWMScore: ",paste(rep(" ",
           times=(widthDisplay-nchar("minPWMScore: ")-nchar(object@maxPWMScore))),collapse=''),
           object@maxPWMScore,"\n"))
         cat("",paste0("lambdaPWM: ",paste(rep(" ",
           times=(widthDisplay-nchar("lambdaPWM: ")-nchar(object@lambdaPWM[1]))),collapse=''),
           object@lambdaPWM,"\n"))
         cat("",paste0("PWMThreshold: ",paste(rep(" ",
             times=(widthDisplay-nchar("PWMThreshold: ")-nchar(object@PWMThreshold))),collapse=''),
             object@PWMThreshold,"\n"))
      }else if(grepl("Occupancy",object@tags)){
        cat(rep("_", widthDisplay),"\n","\n")
        title <- as.character("genomicProfiles object: Occupancy at sites above Threshold")
        titleWidth<-nchar(title, type="byte")
        cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
        cat(rep("_", widthDisplay),"\n","\n")
        profile<-object@profiles

        if(length(profile)<2){
            cat("","Parameter combination: ",names(profile),"\n","\n")
            print(profile)
        }else{
            cat("","Parameter combination: ",names(profile)[1],"\n","\n")
            print(profile[[1]])
            cat("\n")
            cat("",length(profile)-1," remaining combinations","\n","\n")
            cat("",paste0("Parameter combination: ",paste(rep(" ",
              times=(widthDisplay-nchar(paste0(names(profile)[1],"(loci = ",length(profile[[1]]),")")))),collapse=''),
              names(profile)[2:length(profile)]," (loci = ",sapply(profile, length)[2:length(profile)],")\n"))

        }
        cat(rep("_", widthDisplay),"\n","\n")
        title <- as.character("Associated parameters")
        titleWidth<-nchar(title, type="byte")
        cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
        cat(rep("_", widthDisplay),"\n","\n")

        cat("",paste0("lambdaPWM: ",paste(rep(" ",
          times=(widthDisplay-nchar("lambdaPWM: ")-nchar(object@lambdaPWM[1]))),collapse=''),
          object@lambdaPWM,"\n"))
        cat("",paste0("boundMolecules: ",paste(rep(" ",
          times=(widthDisplay-nchar("boundMolecules: ")-nchar(object@boundMolecules[1]))),collapse=''),
          object@boundMolecules,"\n"))

     }else if(grepl("ChIPProfile",object@tags)){
       cat(rep("_", widthDisplay),"\n","\n")
       title <- as.character("genomicProfiles object: ChIP like Profiles")
       titleWidth<-nchar(title, type="byte")
       cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
       cat(rep("_", widthDisplay),"\n","\n")
       profile<-object@profiles

       if(length(profile)<2){
           cat("","Parameter combination: ",names(profile),"\n","\n")
           print(profile)
       }else{
           cat("","Parameter combination: ",names(profile)[1],"\n","\n")
           print(profile[[1]])
           cat("\n")
           cat("",length(profile)-1," remaining combinations","\n","\n")
           cat("",paste0("Parameter combination: ",paste(rep(" ",
             times=(widthDisplay-nchar(paste0(names(profile)[1],"(loci = ",length(profile[[1]]),")")))),collapse=''),
             names(profile)[2:length(profile)]," (loci = ",sapply(profile, length)[2:length(profile)],")\n"))

       }
       cat(rep("_", widthDisplay),"\n","\n")
       title <- as.character("Associated parameter Options")
       titleWidth<-nchar(title, type="byte")
       cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
       cat(rep("_", widthDisplay),"\n","\n")

       cat("",paste0("lambdaPWM: ",paste(rep(" ",
         times=(widthDisplay-nchar("lambdaPWM: ")-nchar(object@lambdaPWM[1]))),collapse=''),
         object@lambdaPWM,"\n"))
       cat("",paste0("boundMolecules: ",paste(rep(" ",
         times=(widthDisplay-nchar("boundMolecules: ")-nchar(object@boundMolecules[1]))),collapse=''),
         object@boundMolecules,"\n"))
        cat("",paste0("stepSize: ",paste(rep(" ",
             times=(widthDisplay-nchar("stepSize: ")-nchar(object@stepSize))),collapse=''),
             object@stepSize,"\n"))
        cat("",paste0("backgroundSignal: ",paste(rep(" ",
           times=(widthDisplay-nchar("backgroundSignal: ")-nchar(object@backgroundSignal))),collapse=''),
           object@backgroundSignal,"\n"))
        cat("",paste0("maxSignal: ",paste(rep(" ",
           times=(widthDisplay-nchar("maxSignal: ")-nchar(object@maxSignal))),collapse=''),
           object@maxSignal,"\n"))
        cat("",paste0("removeBackground: ",paste(rep(" ",
           times=(widthDisplay-nchar("removeBackground: ")-nchar(object@removeBackground))),collapse=''),
           object@removeBackground,"\n"))
        cat("",paste0("chipMean: ",paste(rep(" ",
            times=(widthDisplay-nchar("chipMean: ")-nchar(object@chipMean))),collapse=''),
            object@chipMean,"\n"))
        cat("",paste0("chipSd: ",paste(rep(" ",
             times=(widthDisplay-nchar("chipSd: ")-nchar(object@chipSd))),collapse=''),
             object@chipSd,"\n"))
        cat("",paste0("chipSmooth: ",paste(rep(" ",
            times=(widthDisplay-nchar("chipSmooth: ")-nchar(object@chipSmooth))),collapse=''),
            object@chipSmooth,"\n"))

    } else if(grepl("GoF",object@tags)){
      cat(rep("_", widthDisplay),"\n","\n")
      title <- as.character("genomicProfiles object: Goodness of Fit ")
      titleWidth<-nchar(title, type="byte")
      cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
      cat(rep("_", widthDisplay),"\n","\n")
      profile<-object@profiles

      metrics<- names(profile[[1]][[1]])[!grepl("Mean",names(profile[[1]][[1]]))]
      cat("","Selected Goodness of Fit metrics: ","\n")
      for(i in seq_along(metrics)){
         cat("",metrics[i],"\n")
      }
      cat(rep("_", widthDisplay),"\n","\n")
      cat("",paste(length(profile),"parameter combinations for",length(profile[[1]]),"loci"),"\n")
      cat(rep("_", widthDisplay),"\n","\n")
      if(length(profile)<2){
          cat("","Parameter combination: ",names(profile),"\n","\n")
          print(profile[[1]][[1]][grep("Mean",names(profile[[1]][[1]]))])
          cat("\n",paste(names(profile[[1]]),"\n"))
          cat("\n")
            cat(rep("_", widthDisplay),"\n","\n")
      }else{
          cat("","Parameter combination: ",names(profile)[1],"\n","\n")
          print(profile[[1]][[1]][grep("Mean",names(profile[[1]][[1]]))])
          cat("\n",paste(names(profile[[1]]),"\n"))
          cat("\n")
          cat("",length(profile)-1," remaining combinations","\n","\n")
          cat("",paste0("Parameter combination: ",paste(rep(" ",
            times=(widthDisplay-nchar(paste0(names(profile)[1],"(loci = ",length(profile[[1]]),")")))),collapse=''),
            names(profile)[2:length(profile)]," (loci = ",sapply(profile, length)[2:length(profile)],")\n"))
            cat(rep("_", widthDisplay),"\n","\n")
      }

    }
    }
)


setMethod("show",
    signature = "parameterOptions",
    definition = function(object){
        widthDisplay<-round(options()$width*0.75)
        line1<-class(object)
        line2<-paste("Ploidy: ",paste(rep(" ",
          times=(widthDisplay-nchar("Ploidy: ")-nchar(object@ploidy))),collapse=''),
          object@ploidy,"\n")

        line3<-paste("boundMolecules: ",paste(rep(" ",
          times=(widthDisplay-nchar("boundMolecules: ")-nchar(object@boundMolecules[1]))),collapse=''),
          object@boundMolecules,"\n")

        line4<-paste("maxSignal: ",paste(rep(" ",
          times=(widthDisplay-nchar("maxSignal: ")-nchar(object@maxSignal))),collapse=''),
          object@maxSignal,"\n")

        line5<-paste("backgroundSignal: ",paste(rep(" ",
           times=(widthDisplay-nchar("backgroundSignal: ")-nchar(object@backgroundSignal))),collapse=''),
           object@backgroundSignal,"\n")

        line6<-paste("chipMean: ",paste(rep(" ",
           times=(widthDisplay-nchar("chipMean: ")-nchar(object@chipMean))),collapse=''),
           object@chipMean,"\n")

        line7<-paste("chipSd: ",paste(rep(" ",
            times=(widthDisplay-nchar("chipSd: ")-nchar(object@chipSd))),collapse=''),
            object@chipSd,"\n")

        line8<-paste("chipSmooth:",paste(rep(" ",
           times=(widthDisplay-nchar("chipSmooth:")-nchar(object@chipSmooth))),collapse=''),
           object@chipSmooth,"\n")

        line9<-paste("stepSize: ",paste(rep(" ",
           times=(widthDisplay-nchar("StepSize: ")-nchar(object@stepSize))),collapse=''),
           object@stepSize,"\n")
        line10<-paste("lociWidth: ",paste(rep(" ",
           times=(widthDisplay-nchar("lociWidth: ")-nchar(object@lociWidth))),collapse=''),
           object@lociWidth,"\n")
        line11<-paste("removeBackground: ",paste(rep(" ",
           times=(widthDisplay-nchar("removeBackground: ")-nchar(object@removeBackground))),collapse=''),
          object@removeBackground,"\n")
        line12<-paste("noiseFilter: ",paste(rep(" ",
           times=(widthDisplay-nchar("noiseFilter: ")-nchar(object@noiseFilter))),collapse=''),
           object@noiseFilter,"\n")
        line13<-paste("naturalLog: ",paste(rep(" ",
           times=(widthDisplay-nchar("naturalLog: ")-nchar(object@naturalLog))),collapse=''),
           object@naturalLog,"\n")
        line14<-paste("noOfSites: ",paste(rep(" ",
           times=(widthDisplay-nchar("noOfSites: ")-nchar(object@noOfSites))),collapse=''),
           object@noOfSites,"\n")
        line15<-paste("PWMThreshold: ",paste(rep(" ",
          times=(widthDisplay-nchar("PWMThreshold: ")-nchar(object@PWMThreshold))),collapse=''),
          object@PWMThreshold,"\n")
        line16<-paste("strandRule: ",paste(rep(" ",
          times=(widthDisplay-nchar("strandRule: ")-nchar(object@strandRule))),collapse=''),
          object@strandRule,"\n")
        line17<-paste("whichstrand: ",paste(rep(" ",
           times=(widthDisplay-nchar("whichstrand: ")-nchar(object@whichstrand))),collapse=''),
           object@whichstrand,"\n")
        line18<-paste("PWMpseudocount: ",paste(rep(" ",
           times=(widthDisplay-nchar("PWMpseudocount: ")-nchar(object@PWMpseudocount))),collapse=''),
           object@PWMpseudocount,"\n")
       line19<-paste("lambdaPWM: ",paste(rep(" ",
          times=(widthDisplay-nchar("lambdaPWM: ")-nchar(object@lambdaPWM[1]))),collapse=''),
         object@lambdaPWM,"\n")

      ### conditional prinint
      if(object@paramTag=="empty"){
        cat(rep("_", widthDisplay),"\n")
        title <- line1
        titleWidth<-nchar(title, type="byte")

        cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
        cat(rep("_", widthDisplay),"\n")

        cat("",paste("processingChIP options","\n","\n"),line6,line7,
          line8, line10,line12,"\n")
          cat(rep("_", widthDisplay),"\n")
        cat("",paste("processingChIP options Updated","\n","\n"),line4,line5,"\n")
        cat(rep("_", widthDisplay),"\n")
        cat("",paste("Position Weight Matrix Options","\n","\n"),line13,line14,line18,"\n")
        cat(rep("_", widthDisplay),"\n")
        cat("",paste("Genome Wide Score options","\n","\n"),line16,line17,line19,"\n")
        cat(rep("_", widthDisplay),"\n")
        cat("",paste("PWM Scores above Threshold options","\n","\n"),line16,line17,"\n")
        cat(rep("_", widthDisplay),"\n")
        cat("",paste("Occupancy options","\n","\n"),line2,line19,line3,line4,line5,"\n")
        cat(rep("_", widthDisplay),"\n")
        cat("",paste("ChIP Profile options","\n","\n"),line6,line7,line8,line9,"\n")
        cat(rep("_", widthDisplay),"\n")

      }




    }
)



setMethod("show",
    signature = "ChIPScore",
    definition = function(object){

    widthDisplay<-round(options()$width*0.75)


        scores<-object@scores
        cat(rep("_", widthDisplay),"\n")
        title <- as.character(paste("ChIP score from", length(scores),"regions"))
        titleWidth<-nchar(title, type="byte")

      cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
        cat(rep("_", widthDisplay),"\n")

        if(length(scores)>6){
          topScore <- c(paste(round(scores[[1]][1:6],digits=3),sep='',collapse=' '),
                        paste(round(scores[[2]][1:6],digits=3),sep='',collapse=' '),
                        paste(round(scores[[3]][1:6],digits=3),sep='',collapse=' '))
          len<-length(scores)
          bottomScore <- c(paste(round(scores[[len-2]][1:6],digits=3),sep='',collapse=' '),
                        paste(round(scores[[len-1]][1:6],digits=3),sep='',collapse=' '),
                        paste(round(scores[[len]][1:6],digits=3),sep='',collapse=' '))

           cat("",paste("[",1:3,"] ",head(names(scores),3)," : ",topScore,"...\n"),
               paste("     . . .","\n"),
               paste("     . . .","\n"),
               paste("     . . .","\n"),
               paste("[",(len-2):len,"] ",tail(names(scores),3)," : ",bottomScore,"...\n")
               )
        } else{
           cat(paste("[",seq_along(scores),"] ",names(scores)),"\n")
        }
        cat(rep("_", widthDisplay),"\n","\n")

        cat(rep("_", widthDisplay),"\n","\n")
        title <- paste("Top ", length(object@loci)," regions")
        titleWidth<-nchar(title, type="char")
        cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
        cat(rep("_", widthDisplay),"\n")
        print(object@loci)
        cat(rep("_", widthDisplay),"\n","\n")

        cat(rep("_", widthDisplay),"\n","\n")
        title <- "Associated Options"
        titleWidth<-nchar(title, type="byte")
        cat(rep(" ",ceiling((widthDisplay)/5)),title,"\n")
        cat(rep("_", widthDisplay),"\n")
        line<-paste("","backgroundSignal: ",paste(rep(" ",
          times=(widthDisplay-nchar("backgroundSignal: ")-nchar(object@backgroundSignal))),collapse=''),
          object@backgroundSignal,"\n")
        line2<-paste("maxSignal: ",paste(rep(" ",
          times=(widthDisplay-nchar("maxSignal: ")-nchar(object@maxSignal))),collapse=''),
          object@maxSignal,"\n")
        line3<-paste("chipMean: ",paste(rep(" ",
          times=(widthDisplay-nchar("chipMean: ")-nchar(object@chipMean))),collapse=''),
          object@chipMean,"\n")
        line4<-paste("chipSd: ",paste(rep(" ",
          times=(widthDisplay-nchar("chipSd: ")-nchar(object@chipSd))),collapse=''),
          object@chipSd,"\n")
        line5<-paste("chipSmooth: ",paste(rep(" ",
            times=(widthDisplay-nchar("chipSmooth: ")-nchar(object@chipSmooth))),collapse=''),
            object@chipSmooth,"\n")
        line6<-paste("lociWidth: ",paste(rep(" ",
          times=(widthDisplay-nchar("lociWidth: ")-nchar(object@lociWidth))),collapse=''),
          object@lociWidth,"\n")
        cat(line,line2,line3,line4,line5,line6)
        cat(rep("_", widthDisplay),"\n")

      }

)
