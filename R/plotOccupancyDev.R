###### Ploting occupancy Internal ######
.argsSwap<-function(input, args){

    paramBuffer<-rbind(rep(0,length(args)),args)
    unNamed<-c()
    for(i in seq_along(input)){
        if(names(input)[i] %in% names(args)){
            paramBuffer[,match(names(input)[i],names(args))]<-c(1,input[i])
        }else{
            unNamed<-c(unNamed,i)
        }

    }

    if(length(unNamed)>0){
    paramBuffer[,which(paramBuffer[1,]==0)[length(unNamed)]]<-input[unNamed]
    }
    buffer<-as.vector(paramBuffer[2,])
    names(buffer)<-names(args)
    return(buffer)
}

.prepCS<-function(cs,loc){
    ## cutting cs
    tmp <- data.frame(cs)
    tmp2 <- data.frame(loc)
    if(tmp[1,"start"]< tmp2[1,"start"]){
        tmp[1,"start"]<-tmp2[1,"start"]
    }
    if(tmp[nrow(tmp),"end"]> tmp2[nrow(tmp2),"end"]){
        tmp[nrow(tmp),"end"]<-tmp2[nrow(tmp2),"end"]
    }

    ## replace cs state by numeric
    tmp[,"stateID"]<-gsub("CS","",tmp[,"stateID"])
    return(tmp)
}

.parameterShuffle<-function(args,geneRefColour,csColours,chrinfo,set){
    ## Setting Plot order for name retrieval
    #argsOrder<-c("chromatinState","chipProfile","predictedProfile","occupancy","text",names(GeneRefColour))
    argsOrder<-c("predictedProfile",
                 "chipProfile",
                 "chromatinState",
                 "occupancy",
                 "text",
                 names(geneRefColour),
                 names(csColours))
    ## Setting Defaults for color argument
    colour<-c("#D55E00","#999999","#F0E442","#56B4E9","black",geneRefColour,csColours)
    names(colour)<-argsOrder

    ## Setting default for density
    density<-c(NA,NA,rep(NA,length(geneRefColour)),rep(NA,length(csColours)))
    names(density)<-c("chipProfile","chromatinState",names(geneRefColour),names(csColours))

    ## Setting default borders
    border<-c(NA,NA,rep(NA,length(geneRefColour)),rep(NA,length(csColours)))
    names(border)<-c("chipProfile","chromatinState",names(geneRefColour),names(csColours))

    ## Setting Default for line type
    lineType<-c(1,1,1,1,1,rep(1,length(geneRefColour)),rep(1,length(csColours)))
    names(lineType)<-argsOrder

    ## Setting Default for lwd
    lineWidth<-c(3,0.8,1,2,0.5,rep(1,length(geneRefColour)),rep(1,length(csColours)))
    names(lineWidth)<-argsOrder

    ## setting up font size for plot
    fontsSizePlot<-c(1,1,0.5,0.5,0.5)
    names(fontsSizePlot)<-c("xlab","ylab","strand","geneRef","cs")

    ##setting up font size for Axis
    fontsSizeAxis<-0.8

    ## Setting up labels
    main<-"Occupancy Profile"
    xlab<-paste("DNA Position",chrinfo[[1]][set],paste(chrinfo[[2]][set],":",chrinfo[[4]][set],sep=""),sep=" ")
    ylab<-"Occupancy"


    ## reshaping localPredcitedPropfile to remove 0 on the endPos



    lableOri<-c(1,1)
    names(lableOri)<-c("xaxis","yaxis")
    ## plot limit set up

    xlim<-c(chrinfo[[2]][set], chrinfo[[4]][set])
    xaxislabels<-round(seq(from=xlim[1],to=xlim[2], length.out=10))


    ## Axis set up
    if(any(geneRefColour=="None") & any(csColours =="None")){
        ylim<-c(0,1)
    }else if(!any(geneRefColour=="None") & !any(csColours =="None")){
        ylim<-c(-0.8,1)
        track <- c(-0.2,-0.4)
    }else {
        ylim<-c(-0.4,1)
        track <- c(-0.2)
    }

    ## Swaping if arguments supplied
    if(any(names(args)=="col")){
        buffer<-args$col
        if(is.null(names(buffer))){
            colour[seq_along(buffer)]<-buffer
        }else{
            colour<-.argsSwap(buffer,colour)
        }
    }
    if(any(names(args)=="density")){
        buffer<-args$density
        if(is.null(names(buffer))){
            density[seq_along(buffer)]<-buffer
        }else{
            density<-.argsSwap(buffer,density)
        }
    }
    if(any(names(args)=="border")){
        buffer<-args$border
        if(is.null(names(buffer))){
            border[seq_along(buffer)]<-buffer
        }else{
            border<-.argsSwap(buffer,border)
        }
    }
    if(any(names(args)=="lty")){
        buffer<-args$lty
        if(is.null(names(buffer))){
            lineType[seq_along(buffer)]<-buffer
        }else{
            lineType<-.argsSwap(buffer,lineType)
        }
    }
    if(any(names(args)=="lwd")){
        buffer<-args$lwd
        if(is.null(names(buffer))){
            lineWidth[seq_along(buffer)]<-buffer
        }else{
            lineWidth<-.argsSwap(buffer,lineWidth)
        }
    }
    if(any(names(args)=="cex")){
        buffer<-args$cex
        if(is.null(names(buffer))){
            fontsSizePlot[seq_along(buffer)]<-buffer
        }else{
        fontsSizePlot<-.argsSwap(buffer,fontsSizePlot)
        }
    }
    if(any(names(args)=="cex.axis")){
        fontsSizeAxis<-  args$cex.axis
    }
    if(any(names(args)=="las")){
        buffer<-args$las
        if(is.null(names(buffer))){
            lableOri[seq_along(buffer)]<-buffer
        }else{
        lableOri<-.argsSwap(buffer,lableOri)
        }
    }
    if(any(names(args)=="main")){
        main<-args$main
    }
    if(any(names(args)=="xlab")){
        xlab<-args$xlab
    }
    if(any(names(args)=="ylab")){
        ylab<-args$ylab
    }
    if(any(names(args)=="xlim")){

        buffer<-args$xlim
        xlim[seq_along(buffer)]<-buffer
        xaxislabels<-round(seq(from=xlim[1],to=xlim[2], length.out=10))

    }
    if(any(names(args)=="ylim")){
        buffer<-args$ylim
        ylim[2]<-max(buffer)

    }




    ##Name szwpping and multicolor handeling
    ### Genetic element coding regions
    ### NOTE: Does this need to be done for CS plotting
    names(colour)[grepl("intron",names(colour)) | grepl("exon",names(colour)) | grepl("gene",names(colour),ignore.case=TRUE)]<-"gene"
    colour[names(colour)=="gene"]<-colour[names(colour)=="gene"][1]

    names(lineType)[grepl("intron",names(lineType)) | grepl("exon",names(lineType)) | grepl("gene",names(lineType),ignore.case=TRUE)]<-"gene"
    lineType[names(lineType)=="gene"]<-lineType[names(lineType)=="gene"][1]

    names(lineWidth)[grepl("intron",names(lineWidth)) | grepl("exon",names(lineWidth)) | grepl("gene",names(lineWidth),ignore.case=TRUE)]<-"gene"
    lineWidth[names(lineWidth)=="gene"]<-lineWidth[names(lineWidth)=="gene"][1]

    names(density)[grepl("intron",names(density)) | grepl("exon",names(density)) | grepl("gene",names(density),ignore.case=TRUE)]<-"gene"
    density[names(density)=="gene"]<-density[names(density)=="gene"][1]

    names(border)[grepl("intron",names(border)) | grepl("exon",names(border)) | grepl("gene",names(border),ignore.case=TRUE)]<-"gene"
    border[names(border)=="gene"]<-border[names(border)=="gene"][1]





    ### genetic element UTR
    names(colour)[grepl("UTR",names(colour))]<-"UTR"
    colour[names(colour)=="UTR"]<-colour[names(colour)=="UTR"][1]

    names(lineType)[grepl("UTR",names(lineType))]<-"UTR"
    lineType[names(lineType)=="UTR"]<-lineType[names(lineType)=="UTR"][1]

    names(lineWidth)[grepl("UTR",names(lineWidth))]<-"UTR"
    lineWidth[names(lineWidth)=="UTR"]<-lineWidth[names(lineWidth)=="UTR"][1]

    names(density)[grepl("UTR",names(density))]<-"UTR"
    density[names(density)=="UTR"]<-density[names(density)=="UTR"][1]

    names(border)[grepl("UTR",names(border))]<-"UTR"
    border[names(border)=="UTR"]<-border[names(border)=="UTR"][1]



    return(list("colour"=colour,"density"=density,"border"=border,"lineType"=lineType,"lineWidth"=lineWidth,"fontsSizePlot"=fontsSizePlot,"fontsSizeAxis"=fontsSizeAxis,"lableOri"=lableOri,"main"=main,"xlab"=xlab,"ylab"=ylab,"xlim"=xlim,"ylim"=ylim,"xaxislabels"=xaxislabels,"track" =track))
}






###### Plotiing occupancy Profile #######


plotOccupancyProfile<-function(predictedProfile,
                               ChIPScore = NULL,
                               chromatinState = NULL,
                               asCS = FALSE,
                               occupancy = NULL,
                               goodnessOfFit = NULL,
                               PWM=FALSE,
                               geneRef = NULL,
                               axis=TRUE,...){

     ### handling paramter checks and building up object for plotting
     buffer<- .what.is.predictedProfile(predictedProfile)
     lociLocal<-buffer$loci
     predictedProfileLocal<-buffer$predictedProfile
     stepSize<-buffer$stepSize
     .cleanUpAfterYourself(buffer)

     ## Extracting Chromosome information
     chromosome <- as.character(seqnames(lociLocal))
     startPoint <- start(lociLocal)
     ranges <- ranges(lociLocal)
     endPoint <- end(lociLocal)
     widthSeq <- width(lociLocal)
     chrinfo<-list(chromosome,startPoint,ranges,endPoint,widthSeq)



     if(!is.null(occupancy)){
        occupancyLocal<-.what.is.occupancy(occupancy)
     }


     if(!is.null(ChIPScore)){
         ChIPScoreLocal<-.what.is.ChIPScore(ChIPScore)

         ChIPScoreLocal<-rep(ChIPScoreLocal,length(predictedProfileLocal)/length(ChIPScoreLocal))
     }

     if(!is.null(goodnessOfFit)){
        goodnessOfFitLocal<-.what.is.goodnessOfFit(goodnessOfFit)
     }

    if(!is.null(chromatinState) & asCS == FALSE){
        localchromatineState <- .AccessExtract(lociLocal, chromatinState,asCS)
        csColours<-"None"
        names(csColours)<-"default"
    } else if(!is.null(chromatinState) & asCS == TRUE){
        localchromatineState <- .AccessExtract(lociLocal, chromatinState,asCS)
        csPal <- colorRampPalette(c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
        csColours <- csPal(length(unique(localchromatineState$CS)))
        names(csColours) <- unique(localchromatineState$CS)

    } else{
      localchromatineState <- lociLocal
      names(localchromatineState)<-names(lociLocal)
      csColours<-"None"
      names(csColours)<-"default"
    }



    ## Graphical Parameters set up!
    args<-list(...)


    ## Parsing gene ref if present

    if(!is.null(geneRef) & is(geneRef,"GRanges")){
        ##Seting Graphical bounderies

        # Selecting genes in regions of Interest
        genes <- geneRef[queryHits(findOverlaps(geneRef, lociLocal))]

        # Fitting to window
        start(genes) <- pmax(start(genes), start(lociLocal))
        end(genes) <- pmin(end(genes), end(lociLocal))
        geneRefColour<-heat.colors(length(unique(genes$type)))
        names(geneRefColour)<-unique(genes$type)
        elementBuffer<-split(genes,genes$type)

    } else {
        geneRefColour<-"None"
        names(geneRefColour)<-"default"
    }


    for(i in seq_along(predictedProfileLocal)){
    ### Extracting data to be plotted

    stepIndex <- seq(1, width(lociLocal[i]),by=stepSize)
    localPosition<-seq(start(lociLocal[i]), end(lociLocal[i]),by=stepSize)
    localPosition<-c(localPosition[1]-1,localPosition,localPosition[length(localPosition)]+1)
    localPredcitedPropfile <- predictedProfileLocal[[i]]
    localPredcitedPropfile <- c(0,localPredcitedPropfile,0)

    ## Setting Defaults

    param<-.parameterShuffle(args,geneRefColour,csColours,chrinfo,i)

    ## setting up empty plot with and without Gene rangeBuffer
    if(!is.null(goodnessOfFit)){
        par(xpd=T)
        par(mar=c(6,2,4,10))
    }
    par(family="mono")
    plot(0,type="n", axes=FALSE,xlab="",ylab="",ylim=param$ylim,xlim=c(param$xlim[1],param$xlim[2]))
    title(xlab=param$xlab, cex.lab=param$fontsSizePlot["xlab"])
    title(ylab=param$ylab, cex.lab=param$fontsSizePlot["ylab"],line=0.8)
    title(param$main)
    if(axis){axis(side=BELOW<-1,at=param$xaxislabels,labels=param$xaxislabels,cex.axis=param$fontsSizeAxis)}




    ## Accesibility plotting
    if(!is.null(chromatinState) & asCS == FALSE){
        localAccessibility<-localchromatineState[[i]]

        for(localAccess in seq_len(nrow(localAccessibility))){
            rect(localAccessibility[localAccess,"start"],0,localAccessibility[localAccess,"end"],param$ylim[2],
            col=param$colour[["chromatinState"]],density=param$density[["chromatinState"]], border=param$borders[["chromatinState"]],lwd=param$lineWidth[["chromatinState"]],lty=param$lineType[["chromatinState"]])
        }
    } else if(!is.null(chromatinState) & asCS == TRUE) {
      localAccessibility<-localchromatineState[[i]]

      for(localAccess in seq_len(nrow(localAccessibility))){
          rect(localAccessibility[localAccess,"start"],param$track[1],localAccessibility[localAccess,"end"],param$ylim[2],
          col=param$colour[["chromatinState"]],density=param$density[["chromatinState"]], border=param$borders[["chromatinState"]],lwd=param$lineWidth[["chromatinState"]],lty=param$lineType[["chromatinState"]])
          #text(start(elementBuffer[[elem]])[pos],(trackPosY + 0.03),elementBuffer[[elem]]$ID,pos=4, col=param$colour["text"],cex=param$fontsSizePlot["geneRef"])
      }
    }
    ## ChIP Profile Plotting
    if(!is.null(ChIPScore)){

        chipProfile<-ChIPScoreLocal[[i]]

        localchipProfile <- c(0,chipProfile[stepIndex],0)
        polygon(localPosition,localchipProfile,col=param$colour[["chipProfile"]],density=param$density[["chipProfile"]],border=param$border[["chipProfile"]],lty=param$lineType[["chipProfile"]])
    }
    ## Predicted Profile plotting

    lines(localPosition,localPredcitedPropfile,type="l",col=param$colour[["predictedProfile"]],lwd=param$lineWidth[["predictedProfile"]],lty=param$lineType[["predictedProfile"]])

    ## Occupancy or PWM Scores
    if(!is.null(occupancy)){
        occupancy<-occupancyLocal[[i]]

        if(PWM){
            PWMScaling <- occupancy[head(order(occupancy$PWMScore,decreasing=T), round(0.9*length(occupancy$PWMScore)))]
            ReScale<-((PWMScaling$PWMScore+abs(min(PWMScaling$PWMScore)))/(max(PWMScaling$PWMScore)+abs(min(PWMScaling$PWMScore))))*(param$ylim[2]*0.5)

            lines(x=start(PWMScaling),y=ReScale,type="h",lty=param$lineType[["occupancy"]],
                col=param$colour[["occupancy"]],lwd=param$lineWidth[["occupancy"]])

        }else{

            OccupScaling <- occupancy[head(order(occupancy$Occupancy,decreasing=T), round(0.9*length(occupancy$Occupancy)))]

            ReScale<-(OccupScaling$Occupancy/max(OccupScaling$Occupancy))*(param$ylim[2]*0.5)
            lines(x=start(OccupScaling),y=ReScale,type="h",lty=param$lineType[["occupancy"]],
                col=param$colour[["occupancy"]],lwd=param$lineWidth[["occupancy"]])

        }
    }
    ## Adding accuracy estimate
    if(!is.null(goodnessOfFit)){
      goodnessOfFit<-goodnessOfFitLocal[[i]]

      leg<-paste(names(goodnessOfFit),"=",signif(goodnessOfFit,4))
      legend(x=(param$xlim[2])+0.06*((param$xlim[2])-(param$xlim[1])),y=max(param$ylim)+0.2,
      legend=leg,cex=0.68)

    }

    ## Plotting geneRef
    if(!is.null(geneRef)){
        trackPosY <- ifelse(length(param$track)==1,param$track[1],param$track[2])
        lines(c(min(localPosition),max(localPosition)),c(trackPosY,trackPosY), col = param$colour[["text"]],lty=param$lineType[["text"]], lwd=param$lineWidth[["text"]])
        text(min(localPosition), (trackPosY + 0.1), "+", col=param$colour[["text"]],cex=param$fontsSizePlot["strand"])
        text(min(localPosition), (trackPosY - 0.1), "-", col=param$colour[["text"]],cex=param$fontsSizePlot["strand"])
        ## Genetic Element plotting

        for(elem in seq_along(elementBuffer)){
            ## initialising count for "other elem" elements

            if(names(elementBuffer)[elem]=="intron"){

                pos<-which(as.logical(strand(elementBuffer[[elem]])=="+" | strand(elementBuffer[[elem]])=="*" & strand(elementBuffer[[elem]])!="-"))
                neg<-which(as.logical(strand(elementBuffer[[elem]])=="-" | strand(elementBuffer[[elem]])=="*" & strand(elementBuffer[[elem]])!="+"))
                if(length(pos)>0){
                    segments(x0=start(elementBuffer[[elem]])[pos],x1=end(elementBuffer[[elem]])[pos],y0=(trackPosY + 0.1),lwd=param$lineWidth["gene"],lty=param$lineType["gene"],col=param$colour["gene"])
                }
                if(length(neg)>0){
                    segments(x0=start(elementBuffer[[elem]])[neg],x1=end(elementBuffer[[elem]])[neg],y0=(trackPosY - 0.1),lwd=param$lineWidth["gene"],lty=param$lineType["gene"],col=param$colour["gene"])
                }

            } else if(names(elementBuffer)[elem]=="gene" | names(elementBuffer)[elem]=="exon"){

                pos<-which(as.logical(strand(elementBuffer[[elem]])=="+" | strand(elementBuffer[[elem]])=="*" & strand(elementBuffer[[elem]])!="-"))
                neg<-which(as.logical(strand(elementBuffer[[elem]])=="-" | strand(elementBuffer[[elem]])=="*" & strand(elementBuffer[[elem]])!="+"))
                if(length(pos)>0){
                    rect(start(elementBuffer[[elem]])[pos],(trackPosY + 0.15),end(elementBuffer[[elem]])[pos],(trackPosY + 0.05),col=param$colour["gene"],density=param$density["gene"],lty=param$lineType["gene"],lwd=param$lineWidth["gene"],border=param$border["gene"])
                    text(start(elementBuffer[[elem]])[pos],(trackPosY + 0.03),elementBuffer[[elem]]$ID,pos=4, col=param$colour["text"],cex=param$fontsSizePlot["geneRef"])
                }
                if(length(neg)>0){
                    rect(start(elementBuffer[[elem]])[neg],(trackPosY - 0.05),end(elementBuffer[[elem]])[neg],(trackPosY - 0.15),col=param$colour["gene"],density=param$density["gene"],lty=param$lineType["gene"],lwd=param$lineWidth["gene"],border=param$border["gene"])
                    text(start(elementBuffer[[elem]])[neg],(trackPosY - 0.18),elementBuffer[[elem]]$ID,pos=4, col=param$colour["text"],cex=param$fontsSizePlot["geneRef"])
                }

            } else if(length(grep(x=names(elementBuffer)[elem],pattern="UTR",ignore.case=TRUE))>0){

                pos<-which(as.logical(strand(elementBuffer[[elem]])=="+" | strand(elementBuffer[[elem]])=="*" & strand(elementBuffer[[elem]])!="-"))
                neg<-which(as.logical(strand(elementBuffer[[elem]])=="-" | strand(elementBuffer[[elem]])=="*" & strand(elementBuffer[[elem]])!="+"))

                if(length(pos)>0){
                    rect(start(elementBuffer[[elem]])[pos],(trackPosY + 0.15),end(elementBuffer[[elem]])[pos],(trackPosY + 0.05),col=param$colour["UTR"],density=param$density["UTR"],lty=param$lineType["UTR"],lwd=param$lineWidth["UTR"],border=param$border["UTR"])
                    text(start(elementBuffer[[elem]])[pos],(trackPosY + 0.03),elementBuffer[[elem]]$ID,pos=4, col=param$colour["text"],cex=param$fontsSizePlot["geneRef"])
                }
                if(length(neg)>0){
                    rect(start(elementBuffer[[elem]])[neg],(trackPosY - 0.05),end(elementBuffer[[elem]])[neg],(trackPosY - 0.15),col=param$colour["UTR"],density=param$density["UTR"],lty=param$lineType["UTR"],lwd=param$lineWidth["UTR"],border=param$border["UTR"])
                    text(start(elementBuffer[[elem]])[neg],(trackPosY - 0.18),elementBuffer[[elem]]$ID,pos=4, col=param$colour["text"],cex=param$fontsSizePlot["geneRef"])
                }

            } else {
                # Shifting through multiple colors if different other elem elements
                  nameBuffer<-names(elementBuffer)[elem]
                  pos<-which(as.logical(strand(elementBuffer[[elem]])=="+" | strand(elementBuffer[[elem]])=="*" & strand(elementBuffer[[elem]])!="-"))
                  neg<-which(as.logical(strand(elementBuffer[[elem]])=="-" | strand(elementBuffer[[elem]])=="*" & strand(elementBuffer[[elem]])!="+"))

                  if(length(pos)>0){
                      rect(start(elementBuffer[[elem]])[pos],(trackPosY + 0.15),end(elementBuffer[[elem]])[pos],(trackPosY + 0.05),col=param$colour[nameBuffer],density=param$density[nameBuffer],lty=param$lineType[nameBuffer],lwd=param$lineWidth[nameBuffer],border=param$border[nameBuffer])
                      text(start(elementBuffer[[elem]])[pos],(trackPosY + 0.03),elementBuffer[[elem]]$ID,pos=4, col=param$colour["text"],cex=param$fontsSizePlot["geneRef"])
                  }
                  if(length(neg)>0){
                      rect(start(elementBuffer[[elem]])[neg],(trackPosY - 0.05),end(elementBuffer[[elem]])[neg],(trackPosY - 0.15),col=param$colour[nameBuffer],density=param$density[nameBuffer],lty=param$lineType[nameBuffer],lwd=param$lineWidth[nameBuffer],border=param$border[nameBuffer])
                      text(start(elementBuffer[[elem]])[neg],(trackPosY - 0.18),elementBuffer[[elem]]$ID,pos=4, col=param$colour["text"],cex=param$fontsSizePlot["geneRef"])
                  }
            }
        }
    }
  }
}
