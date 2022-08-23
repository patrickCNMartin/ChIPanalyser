# Cleaner plotting function for occupancy profiles. 

## Plotting occupency main function 

plotOccupancyProfile <- function(predictedProfile,
    ChIPScore = NULL,
    chromatinState = NULL,
    occupancy = NULL,
    goodnessOfFit = NULL,
    PWM = FALSE,
    geneRef = NULL,
    addLegend = TRUE,
    ...){
    ## First we check that all objects contain what they need to contain 
    ### handling paramter checks and building up object for plotting
    buffer<- .what.is.predictedProfile(predictedProfile)
    lociLocal<-buffer$loci
    predictedProfileLocal<-buffer$predictedProfile
    stepSize<-buffer$stepSize
    .cleanUpAfterYourself(buffer)

    ## Extracting Chromosome information
    chromosome <- as.character(seqnames(lociLocal))
    startPoint <- start(lociLocal)
    #ranges <- ranges(lociLocal)
    endPoint <- end(lociLocal)
    

    ## Exrtracting Occupancy if present 
    if(!is.null(occupancy)){
        occupancyLocal<-.what.is.occupancy(occupancy)
        occup <- TRUE
    } else {
        occup <- FALSE
    }
    
    ## Extracting ChIPscore if present 
    if(!is.null(ChIPScore)){
        ChIPScoreLocal <- .what.is.ChIPScore(ChIPScore)
        ChIPScoreLocal <- rep(ChIPScoreLocal,
            length(predictedProfileLocal)/length(ChIPScoreLocal))
        chip <- TRUE
    } else {
        chip <- FALSE
    }
    ## Extrating GoF score if present 
    if(!is.null(goodnessOfFit)){
        goodnessOfFitLocal<-.what.is.goodnessOfFit(goodnessOfFit)
        gof <- TRUE
    } else {
        gof <- FALSE
    }

    ## Extracting chromatin states if present 
    if(!is.null(chromatinState)){
        chromaState <- .what.is.chromatinState(chromatinState,lociLocal)
        cs <- TRUE
    } else {
        cs <- FALSE
    }

    ## Extracting geneRef if present 
    if(!is.null(geneRef)){
        genes <- .what.is.geneRef(geneRef,lociLocal)
        gr <- TRUE
    } else {
        gr <- FALSE
    }

    # Get elipsis arguments 
    graphical <- list(...)
    ### Looping over profiles 
    for(i in seq_along(predictedProfileLocal)){
        # Get local graphical paramters 
        param <- .dispatchGraphical(graphical,
            xlim = c(startPoint[i],endPoint[i]),
            chr = chromosome[i],
            cs = cs,
            geneRef = gr,
            occupancy = occup,
            gof = gof,
            legend = addLegend)
        # Extracting data to be plotted 
        stepIndex <- seq(1, width(lociLocal[i]),by=stepSize)
        localPosition<-seq(start(lociLocal[i]), end(lociLocal[i]),by=stepSize)
        localPosition<-c(localPosition[1]-1,localPosition,localPosition[length(localPosition)]+1)
        localPredcitedPropfile <- predictedProfileLocal[[i]]
        localPredcitedPropfile <- c(0,localPredcitedPropfile,0)

        ## Plot empty window 
        .drawBackground(param)

         ## Adding ChIPscore 
        if(chip){
            chipProfile<-ChIPScoreLocal[[i]]
            localchipProfile <- c(0,chipProfile[stepIndex],0)
            .drawChIP(localPosition,localchipProfile,param)
        }
        ## Adding Occupancy 
        if(occup){
            .drawOccup(occupancyLocal[[i]],param,PWM)
        }
        ## Adding Predicted 
        .drawPrediction(localPosition,localPredcitedPropfile,param)

         ## Adding chromatin states 
        if(cs){
            chroma <- .prepCS(chromaState[[i]],param)
            .drawCS(chroma)
        }

        ## Adding gene reference 
        if(gr){
            genes <- .prepGR(genes[[i]],param)
            .drawGeneRef(genes)
        }

        if(addLegend){
           
           leg <- .prepLegend(chip = chip,
                occup = occup,
                cs = ifelse(cs,list(chroma),"empty"),
                gof = ifelse(gof,goodnessOfFitLocal[i],"empty"),
                param)
            .drawLegend(leg)
        }
    }

    
}




## Creating param table/list for new plots and segment plots 
# we just want to make it easier and cleaner to add colors and what not
.dispatchGraphical <- function(graph,
    xlim = NULL,
    chr = NULL,
    cs = FALSE,
    geneRef = FALSE,
    occupancy = FALSE,
    gof = FALSE,
    legend = FALSE){
    ## First set up parameters for empty plot.
    param <- list()

    ## Expanding if legened needs to be added
    if(legend & cs & gof){
        param$xpd <- TRUE 
        param$mar <- c(6,8,4,15)
    } else if(legend & !cs & gof){
        param$xpd <- TRUE 
        param$mar <- c(6,8,4,10)
    } else if(legend & !cs & !gof){
        param$xpd <- TRUE 
        param$mar <- c(6,2,4,10)
    }
    ## Expanding ylim to have space for geneRef and/or CS
    if(cs & geneRef){
        param$ylim <- c(-0.8,1)
    }else if((cs & !geneRef) | (!cs & geneRef)){
        param$ylim <- c(-0.3,1)
    }else{
        param$ylim <-c(0,1)
    }
    # Adding xlims to param 
    param$xlim <- xlim
    # Dispatching font size 
    param$cex <- ifelse(any(names(graph) == "cex"),graph$cex,1)
    param$cex.lab <- ifelse(any(names(graph) == "cex.lab"),graph$cex.lab,1)
    # Dispatch main title font 
    param$cex.main <- ifelse(any(names(graph) == "cex.main"),graph$cex.main,1)
    # Dispatch CS density 
    param$densityCS <- ifelse(any(names(graph) == "densityCS"),graph$densityCS,50)
    # Dispatch geneRef Density 
    param$densityGR <- ifelse(any(names(graph) == "densityGR"),graph$densityGR,50)
    # Dispatch pred color 
    param$colPred <- ifelse(any(names(graph) == "colPred"),graph$colPred,"#E69F00")
    # Dispatch ChIPcolor 
    param$colChIP <- ifelse(any(names(graph) == "colPred"),graph$colChIP,"#999999")
    # Dispatch occup color 
    param$colOccup <- ifelse(any(names(graph) == "colOccup"),graph$colChIP,"#56B4E9")
    # Dispatch CS color
    param$colCS <- ifelse(any(names(graph) == "colPred"),
        colorRampPalette(graph$colCS),
        colorRampPalette(c("#F0E442","#999999","#E69F00", "#56B4E9", "#009E73",
            "#0072B2", "#D55E00", "#CC79A7")))
    # Dispatch Gene Ref color
    param$colGR <- ifelse(any(names(graph) == "colGR"),
        colorRampPalette(graph$colGR),
        colorRampPalette(brewer.pal(8, "Blues")))
    # Dispatch xlab names 
    param$xlab <- ifelse(any(names(graph) == "xlab"),graph$xlab,
        paste("Occupancy at Position",chr,paste(xlim[1],":",xlim[2],sep=""),sep=" "))
    # Dispatch ylab 
    param$ylab <- ifelse(any(names(graph) == "ylab"),graph$ylab," ")
    # Dispatch axis ticks 
    param$axis <- ifelse(any(names(graph) == "n_axis_ticks"),
        round(seq(from=xlim[1],to=xlim[2], length.out=n_axis_ticks)),
        round(seq(from=xlim[1],to=xlim[2], length.out=10)))
    
    return(param)

}



.drawBackground <- function(param){
    ## setting up empty plotting window only using the relevant parameters.
    ## Checking if the legend needs to be plotted 
    
    par(xpd = param$xpd)
    par(mar = param$mar)
    
    plot(0, type = "n",
        axes = FALSE,
        xlab = "", ylab = "",
        ylim = param$ylim,
        xlim = param$xlim)
    axis(side=BELOW<-1,at=param$axes,labels=param$axes,cex.axis=param$cex)
    title(xlab = param$xlab , cex.lab = param$cex.lab)
    title(ylab = param$ylab , cex.lab = param$cex.lab)
    title(main = param$main , cex = param$cex.main)
}

.prepCS <- function(CS, param){
    CS$y0 <- -0.3
    CS$y1 <- -0.1
    CS$density <- param$densityCS
    col_local <- param$colCS(length(unique(CS$stateID)))
    CS$col <- col_local[match(CS$stateID,unique(CS$stateID))]
    return(CS)
}

.drawCS <- function(CS){
    # looping over blocks 
    for(block in seq_len(nrow(CS))){
        rect(CS$start[block],
            CS$y0[block],
            CS$end[block],
            CS$y1[block],
            density = CS$density[block],
            col = CS$col[block])
    }
}

.drawChIP <- function(x,y,param){
    polygon(x,
        y,
        col = param$colChIP,
        border = NA,
        density = param$densityChIP)
}

.drawPrediction <- function(x,y,param){
    lines(x,
        y,
        type = "l",
        col = param$colPred)
}

.drawOccup <- function(occupancy,param,PWM){
    if(PWM){
        PWMScaling <- occupancy[head(order(occupancy$PWMScore,decreasing=T), 
            round(0.9*length(occupancy$PWMScore)))]
        ReScale<-((PWMScaling$PWMScore+abs(min(PWMScaling$PWMScore)))/
            (max(PWMScaling$PWMScore)+abs(min(PWMScaling$PWMScore))))*(param$ylim[2]*0.5)
        lines(x=start(PWMScaling),y=ReScale,type="h",col = param$colOccup,lwd =2)
    }else{
        OccupScaling <- occupancy[head(order(occupancy$Occupancy,decreasing=T), 
            round(0.9*length(occupancy$Occupancy)))]
        ReScale<-(OccupScaling$Occupancy/max(OccupScaling$Occupancy))*
            (param$ylim[2]*0.5)
        lines(x=start(OccupScaling),y=ReScale,type="h",col = param$colOccup,lwd =2)

    }
}

.prepGR <- function(gr, param){
   elem <- c("chr",param$xlim[1],param$xlim[2],NA,NA,"line")
   gr <- rbind(elem,gr)
   colnames(gr) <- c("chr","x0","x1","width","strand","element")
   elemY0 <- param$ylim[1]
   elemY0[is.na(gr$strand)] <- param$ylim[1] + 0.2
   elemY0[gr$strand == "+"] <- param$ylim[1] + 0.3
   elemY0[gr$strand == "-"] <- param$ylim[1] + 0.1

   elemY1 <- param$ylim[1]
   elemY1[is.na(gr$strand)] <- param$ylim[1] + 0.2
   elemY1[gr$strand == "+"] <- param$ylim[1] + 0.35
   elemY1[gr$strand == "-"] <- param$ylim[1] + 0.15

   gr$y0 <- elemY0
   gr$y1 <- elemY1

   ## Setting colours 
   ids <- match(gr$element,unique(gr$element))
   gr$col <- param$colGR(length(unique(gr$element)))[ids]
   gr$cex <-c(param$cex.lab,rep(gr$cex,nrow(gr)-1))
   return(gr)
}

.drawGeneRef  <- function(gr){
    for(seg in seq_len(nrow(gr))){
        if(gr$element[seg] == "line"){
            segments(gr$x0[seg],
                gr$x1[seg],
                gr$y0[seg],
                col = gr$col[seg])
            text(gr$x0[seg],
                 gr$y0[seg] + 0.15,
                "+",
                col = gr$col[seg],
                cex = gr$cex[seg])
            text(gr$x1[seg],
                gr$y1[seg] - 0.15,
                "-",
                col = gr$col[seg],
                cex = gr$cex[seg])
        }else if(gr$element[seg] == "intron"){
             segments(gr$x0[seg],
                gr$x1[seg],
                gr$y0[seg],
                col = gr$col[seg])
        } else {
            rect(gr$x0[seg],
                gr$y0[seg],
                gr$x1[seg],
                gr$y1[seg],
                col = gr$col[seg])
            text(gr$x0[seg],
                gr$y0[seg],
                gr$element[seg],
                pos = 4,
                cex = gr$cex[seg])
        }
    }
}

.prepLegend <- function(chip,
    occup,
    cs,
    gof,
    param){
    leg <- list()
    leg$pred <- list("x" = param$xlim[2],
        "y" = 0.5,
        "legend" = "Predicted Profile",
        "col" = param$colPred,
        "fill" = NA,
        "border" = NA,
        "cex" = param$cex)
    if(chip){
        leg$chip <- list("x" = param$xlim[2],
        "y" = 0.75,
        "legend" = "ChIP Profile",
        "col" = NA,
        "fill" = param$colChIP,
        "border" = NA,
        "cex" = param$cex)
    }
    if(occup){
        leg$occup <- list("x" = param$xlim[2],
        "y" = 0.25,
        "legend" = "Binding Site",
        "col" = param$colOccup,
        "fill" = NA,
        "border" = NA,
        "cex" = param$cex)
    }
    if(any(gof != "empty")){
        gof <- gof[[1]]
        tmp <- gof[names(gof) %in% c("MSE","AUC","pearson")]
        leg$gof <- list("x" = param$xlim[1],
        "y" = 0,
        "legend" = paste(paste(names(tmp),"=",signif(tmp,4)),collapse = "\n"),
        "col" = "black",
        "fill" = NA,
        "border" = NA,
        "cex" = param$cex)
    }
    if(any(cs != "empty")){
        cs <- cs[[1]]
        leg$cs <- list("x" = param$xlim[2] +((param$xlim[2] - param$xlim[1]) * 0.3),
        "y" = 0.0,
        "legend" = unique(cs$stateID),
        "col" = NA,
        "fill" = unique(cs$col),
        "border" = NA,
        "cex" = param$cex)
    }
    
    return(leg)

}
.drawLegend <- function(leg){
    for(l in seq_along(leg)){
        if(names(leg)[l] == "gof"){
            
            text(x = leg[[l]]$x,
            y = leg[[l]]$y,
            labels = leg[[l]]$legend,
            adj = c(1,0),
            col = leg[[l]]$col,
            cex = leg[[l]]$cex)
        }else {
            legend(x = leg[[l]]$x,
            y = leg[[l]]$y,
            legend = leg[[l]]$legend,
            bty = "n",
            col = leg[[l]]$col,
            fill = leg[[l]]$fill,
            border = leg[[l]]$border,
            lwd = 2,
            xjust = 0,
            yjust = 0.5,
            cex = leg[[l]]$cex)

        }
       
    }
}