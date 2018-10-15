##############################################
#############      S4 Methods    #############
##############################################

########### Show Methods ##########
setMethod("show",
    signature = "genomicProfileParameters",
    definition = function(object){
        message("Object Class:", class(object), "\n", sep = " ")
        message("\n","PWM:", "\n");print(object@PWM)
        message("\n","PFM:", "\n");print(object@PFM)
        message("\n","PFMFormat: ", object@PFMFormat,"\n")
        message("\n","PWM Scores at Sites higher than Threshold: \n")
        print(object@AllSitesAboveThreshold)
        message("\n","No Accessible DNA at Loci:","\n")
        if(length(object@NoAccess)>20){
            rem<-length(object@NoAccess) -20
            message(head(object@NoAccess,20),"\n","\n",rem," Loci Omitted ")
        }else {
            message(object@NoAccess)
        }
        message("\n","Genomic Profile Parameters:","\n")
        cat("Lambda:", object@ScalingFactorPWM,"\n", sep="\t")
        cat("BP Frequency:", object@BPFrequency,"\n",sep="\t")
        message("Pseudocount: ", object@PWMpseudocount,"\n",
            "Natural log: ", object@naturalLog,"\n",
            "Number Of Sites: ",object@noOfSites,"\n",
            "maxPWMScore: ",object@maxPWMScore,"\n",
            "minPWMScore: ",object@minPWMScore, "\n",
            "PWMThreshold: ",object@PWMThreshold)
        cat("Average Exponential PWM Score: ", object@averageExpPWMScore,
            "\n",sep="\t")
        message("DNA Sequence Length: ",object@DNASequenceLength,"\n",
            "Strand Rule: ",object@strandRule,"\n",
            "Strand: ",object@whichstrand,"\n")
    }
)

setMethod("show",
    signature = "occupancyProfileParameters",
    definition = function(object){
        message("Object Class:", class(object), "\n", sep = " ")
        message( "Ploidy: ", object@ploidy )
        cat("boundMolecules: ", object@boundMolecules,"\n",sep="\t")
        message("backgroundSignal: ", object@backgroundSignal,"\n",
            "maxSignal: ", object@maxSignal,"\n",
            "chipMean: ",object@chipMean,"\n",
            "chipSd: ",object@chipSd,"\n",
            "chipSmooth: ", object@chipSmooth,"\n",
            "Step Size: ",object@stepSize,"\n")
    }
)
