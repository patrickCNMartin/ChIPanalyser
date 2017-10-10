plotOccupancyProfile<-function(predictedProfile,setSequence,
    profileAccuracy = NULL, chipProfile = NULL, occupancy = NULL,PWM=FALSE,
    DNAAccessibility = NULL, occupancyProfileParameters = NULL,geneRef = NULL){

    # restricting setSequence to sequences present in profile
    ##( case of non accessible DNA)


    if(!is.null(DNAAccessibility)){
        localAccessibility <- intersect(setSequence, DNAAccessibility)
    }
    names(localAccessibility) <- names(setSequence)


    #Extracting chromosome Information.
    chromosome <- as.character(seqnames(setSequence))
    startPoint <- start(setSequence)
    ranges <- ranges(setSequence)
    endPoint <- end(setSequence)
    widthSeq <- width(setSequence)

    stepSize <- stepSize(occupancyProfileParameters)

    # Extracting gene info (intron, exon, UTR's and enhancers) from geneRef
    if(!is.null(geneRef)) {
        geneList <- vector("list", length(setSequence))
        enhancerList <- vector("list", length(setSequence))
        if(length(which(names(geneRef)=="exon"))>0){
            idsExons <- intersect(intersect(which(start(
                geneRef[["exon"]])<endPoint),
                which(end(geneRef[["exon"]])>=startPoint)),
                which(as.vector(seqnames(geneRef[["exon"]])) ==
                chromosome))
            ids5UTR <- intersect(intersect(which(start(
                geneRef[["5UTR"]])<endPoint),
                which(end(geneRef[["5UTR"]])>=startPoint)),
                which(as.vector(seqnames(geneRef[["5UTR"]])) ==
                chromosome))
            ids3UTR <- intersect(intersect(which(start(
                geneRef[["3UTR"]])<endPoint),
                which(end(geneRef[["3UTR"]])>=startPoint)),
                which(as.vector(seqnames(geneRef[["3UTR"]])) ==
                chromosome))
            geneList <- list("location"=c(geneRef[["exon"]][idsExons],
                geneRef[["5UTR"]][ids5UTR],geneRef[["3UTR"]][ids3UTR]),
                "name"=c(names(geneRef$exon[idsExons]),
                names(geneRef[["5UTR"]][ids5UTR]),
                names(geneRef[["3UTR"]][ids3UTR])))
        }
        if(length(which(names(geneRef)=="enhancer"))>0){
            ids <- intersect(intersect(which(start(
                geneRef[["enchancer"]])<endPoint),
                which(end(geneRef[["enchancer"]])>=startPoint)),
                which(as.vector(seqnames(geneRef[["enhancer"]])) ==
                chromosome))
            enhancerList <- list("location"=geneRef$enhancer[ids],
                "name"=names(geneRef$enhancer[ids]))
        }
    }
    #Extracting postions and predictedProfile
    stepIndex <- seq(1, widthSeq,by=stepSize)
    localPosition <- 1:widthSeq
    localPosition <- c(0,1:length(localPosition[stepIndex]),
        length(localPosition[stepIndex])+1)
    localPredcitedPropfile <- predictedProfile$ChIP
    localPredcitedPropfile <- c(0,localPredcitedPropfile,0)


    xlabels<-paste("DNA Position",chromosome,
        paste(startPoint,":",endPoint,sep=""),sep=" ")
    xaxislabels<-seq(from=startPoint,to=endPoint, by=1000)
    xaxis<-seq(from=1,to=widthSeq,by=1000)/stepSize




    if(length(chipProfile)>1){
        localchipProfile <- chipProfile[stepIndex]
        localchipProfile <- c(0,localchipProfile,0)
        yaxislabels <- signif(seq(min(localchipProfile),
        max(localchipProfile), length.out=6),3)

        if(!is.null(geneRef)){
            ylimMin<- (-0.5)
            } else {
            ylimMin<- 0
            }

        ylimMax <- 1

        plot(localPosition,localchipProfile,type="n",axes=FALSE,
            xlab="",ylab="",ylim=c(ylimMin,ylimMax))
        axis(side=BELOW<-1,at=xaxis,labels=xaxislabels)
        #axis(side=LEFT<-2, at=yaxislabels, labels=yaxislabels,las=1,
        #cex.axis=0.8)
        box()
        title(xlab=xlabels)

        ### DNA Access
        if(length(localAccessibility)>1){
            localAccesStart <- start(localAccessibility)
            localAccesStart <- c(localAccesStart[2:length(localAccesStart)],
                endPoint)
            localAccesEnd <- end(localAccessibility)

            for(n in 1:length(localAccesStart)){
            rect(round((localAccesEnd[n]-startPoint)/stepSize),
                min(localPredcitedPropfile),
                round((localAccesStart[n]-startPoint)/stepSize),
                ylimMax*0.8,col="#F0E442", border = NA)
            }
        }

        polygon(localPosition,localchipProfile,col="#999999",border=NA)
        lines(localPosition,localPredcitedPropfile,col="#D55E00",
            lwd=2,type="l")
        } else {
        yaxislabels <- signif(seq(min(localPredcitedPropfile),
            max(localPredcitedPropfile), length.out=6),3)
        ylimMax <- max(localPredcitedPropfile)+
            (max(localPredcitedPropfile)/5)
        if(!is.null(geneRef)){
            ylimMin <- (-0.5)
        } else {
            ylimMin<- 0
        }

        plot(localPosition,localPredcitedPropfile,type="n",axes=FALSE,xlab="",
            ylab="",ylim=c(ylimMin,ylimMax))
        axis(side=BELOW<-1,at=xaxis,labels=xaxislabels)
        axis(side=LEFT<-2, at=yaxislabels, labels=yaxislabels, las =1,
            cex.axis=0.8)
        box()

        title(xlab=xlabels)
        lines(localPosition,localPredcitedPropfile,xlab=xlabels,type="l",
            lwd=2,col="#D55E00")
        polygon(localPosition,localPredcitedPropfile,col="#D55E00",border=NA)
        }
    ## Occupnacy plotting
    if(!is.null(occupancy)){
        if(PWM=="FALSE"){
        occupancyBuffer <- occupancy[which(
            occupancy$Occupancy > (ylimMax*0.9)*0.05)]
            for(k in 1:length(occupancyBuffer)){
            lines(x=round((start(occupancyBuffer)-startPoint)/stepSize),
            y=(occupancyBuffer$Occupancy),type="h",col="#56B4E9",lwd=2)
            }
        } else {
        PWMScaling <- (occupancy$PWMScore + min(occupancy$PWMScore))/
            max(occupancy$PWMScore)
        PWMBuffer <- occupancy[which(occupancy$PWMScore > PWMScaling)]
            for(k in 1:length(PWMBuffer)){
            lines(x=round((start(PWMBuffer)-startPoint)/stepSize),
            y=(PWMBuffer$PWMScore),type="h",col="#56B4E9",lwd=2)
            }
        }
    }

    ##profileAccuracy
    if(!is.null(profileAccuracy)){
        corr <- round(profileAccuracy["Corr"],3)
        MSE <- round(profileAccuracy["MSE"],3)
        text(max(localPosition)*0.9,ylimMax,
            paste0("Correlation: ",corr),col="black")
        text(max(localPosition)*0.9,ylimMax*0.95,
            paste0("MSE: ",MSE), col="black")
    }
    ### geneRef plotting
    if(!is.null(geneRef)){
        ylimMin <- (-0.5)
        for(geneName in unique(geneList$name)){
            lastExonXend <- (-1)
            lastExonYmed <- (-1)
            firstExonXStart <- (-1)
            for(k in which(geneList$name==geneName)){
                exonXStart <- round((start(geneList$location[k])-
                startPoint)/stepSize)
                exonXEnd <- round((end(geneList$location[k])-
                startPoint)/stepSize)
                exonYStart<- ylimMin*0.25
                if(as.vector(strand(geneList$location[k])) == "+"){
                    exonYStart <- ylimMin*0.25
                }
                exonYEnd<-exonYStart+(ylimMin*0.15)
                rect(exonXStart, exonYStart, exonXEnd, exonYEnd,
                col = "#009E73", border = "#009E73",lwd = 0)
                if(lastExonXend >= 0){
                    lines(c(lastExonXend,exonXStart),c(lastExonYmed,
                    ((exonYStart + exonYEnd)/2)), col = "#009E73",
                    lty=1, lwd=0.8)
                } else{
                    firstExonXStart <- exonXStart
                }
                lastExonXend <- exonXEnd
                lastExonYmed <- (exonYStart + exonYEnd)/2;
            }

        textPos<- ylimMin*0.15
        text(((firstExonXStart+lastExonXend)/2), (textPos), geneName,
        col="black")
        }
        ## strand plotting
    lines(c(min(localPosition),max(localPosition)),c((ylimMin*0.5),
    (ylimMin*0.5)), col = "black",lty=1)
    text(min(localPosition), (ylimMin*0.25), "+", col="black")
    text(min(localPosition), (ylimMin*0.75), "-", col="black")
    ##Emhancer plotting
    if(length(enhancerList)>0){
        for(enhancerName in unique(enhancerList$name)){
        lastEnhancerXend <- (-1)
        lastEnhancerYmed <- (-1)
        firstEnhancerXStart <- (-1)
        for(k in which(enhancerList$name==enhancerName)){
            enhancerXStart <- round((start(
            enhancerList$location[k])-startPoint)/stepSize)
            enhancerXEnd <- round((end(
            enhancerList$location[k])-startPoint)/stepSize)
            enhancerYStart<- (ylimMin*0.25)
            enhancerYEnd<-enhancerYStart+(ylimMin*0.1)
            rect(enhancerXStart, enhancerYStart, enhancerXEnd,
                enhancerYEnd, col = "#CC79A7", border = "#CC79A7",
                lwd = 0)
            if(lastEnhancerXend >= 0){
                lines(c(lastEnhancerXend,enhancerXStart),
                c(lastEnhancerYmed,((enhancerYStart+enhancerYEnd)/2)),
                col = "#CC79A7",lty=1, lwd=0.8)
            } else{
                firstEnhancerXStart <- enhancerXStart
            }
            lastEnhancerXend <- enhancerXEnd
            lastEnhancerYmed <- (enhancerYStart + enhancerYEnd)/2;
        }
    # gene name plotting
        textPos <- ylimMin*0.15
        text(((firstEnhancerXStart+lastEnhancerXend)/2), (textPos),
        enhancerName, col="black")

        }
        if(!is.null(geneList)){
            lines(c(min(localPosition),max(localPosition)),
            c((ylimMin*0.5),(ylimMin*0.5)),
            col = "black",lty=1, lwd=0.5)
            text(min(localPosition), (ylimMin*0.25), "+",
            col="black")
            text(min(localPosition), (ylimMin*0.75), "-",
            col="black")
        }
    }
    }
}
