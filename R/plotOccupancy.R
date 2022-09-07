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
        states <- unique(unlist(sapply(chromaState,function(x){return(x$stateID)})))
        cs <- TRUE
    } else {
        cs <- FALSE
    }

    ## Extracting geneRef if present 
    if(!is.null(geneRef)){
        genes <- .what.is.geneRef(geneRef,lociLocal)
        types <- unique(unlist(sapply(genes,function(x){return(x$type)})))
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
            chroma <- .prepCS(chromaState[[i]],states,param)
            .drawCS(chroma)
        }

        ## Adding gene reference 
        if(gr){
            genelist <- .prepGR(genes[[i]],types,param)
            .drawGeneRef(genelist)
        }

        if(addLegend){
           
           leg <- .prepLegend(chip = chip,
                occup = occup,
                cs = ifelse(cs,list(chroma),"empty"),
                gof = ifelse(gof,goodnessOfFitLocal[i],"empty"),
                gr = ifelse(gr,list(genelist),"empty"),
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
        param$mar <- c(6,8,4,20)
    } else if(legend & !cs & gof){
        param$xpd <- TRUE 
        param$mar <- c(6,8,4,10)
    } else if(legend & !cs & !gof){
        param$xpd <- TRUE 
        param$mar <- c(6,2,4,10)
    } else if(legend & !cs & geneRef){
        param$xpd <- TRUE 
        param$mar <- c(6,8,4,10)
    }
    ## Expanding ylim to have space for geneRef and/or CS
    if(cs & geneRef){
        param$ylim <- c(-1.2,1)
    }else if((cs & !geneRef)){
        param$ylim <- c(-0.3,1)
    }else if(!cs & geneRef){
        param$ylim <- c(-0.9,1)
    }else {
        param$ylim <-c(0,1)
    }
    # Adding xlims to param 
    param$xlim <- xlim
    # Dispatching font size 
    param$cex <- ifelse(any(names(graph) == "cex"),graph$cex,1)
    param$cex.lab <- ifelse(any(names(graph) == "cex.lab"),graph$cex.lab,1)
    # Dispatch main title font 
    param$cex.main <- ifelse(any(names(graph) == "cex.main"),graph$cex.main,1)
    # Dispatch prediction line width 
    param$lwd <- ifelse(any(names(graph) == "lwd"),graph$lwd,1)
    # Dispatch CS density 
    param$densityCS <- ifelse(any(names(graph) == "densityCS"),graph$densityCS,50)
    # Dispatch geneRef Density 
    param$densityGR <- ifelse(any(names(graph) == "densityGR"),graph$densityGR,50)
    # Dispatch pred color 
    param$colPred <- ifelse(any(names(graph) == "colPred"),graph$colPred,"#E69F00")
    # Dispatch ChIPcolor 
    param$colChIP <- ifelse(any(names(graph) == "colChIP"),graph$colChIP,"#999999")
    # Dispatch occup color 
    param$colOccup <- ifelse(any(names(graph) == "colOccup"),graph$colChIP,"#56B4E9")
    # Dispatch CS color
    param$colCS <- ifelse(any(names(graph) == "colCS"),
        colorRampPalette(graph$colCS),
        colorRampPalette(c("#F0E442","#999999","#E69F00", "#56B4E9", "#009E73",
            "#0072B2", "#D55E00", "#CC79A7")))
    # Dispatch Gene Ref color
    param$colGR <- ifelse(any(names(graph) == "colGR"),
        colorRampPalette(graph$colGR),
        colorRampPalette(brewer.pal(8, "Accent")))
    # Dispatch xlab names 
    param$xlab <- ifelse(any(names(graph) == "xlab"),graph$xlab,
        paste("Occupancy at Position",chr,paste(xlim[1],":",xlim[2],sep=""),sep=" "))
    # Dispatch ylab 
    param$ylab <- ifelse(any(names(graph) == "ylab"),graph$ylab," ")
    # Dispatch Main 
    param$main <- ifelse(any(names(graph) == "main"),graph$main," ")
    # Dispatch axis ticks 
    param$n_axis_ticks <- ifelse(any(names(graph) == "n_axis_ticks"),
        list(round(seq(from=xlim[1],to=xlim[2], length.out=graph$n_axis_ticks))),
        list(round(seq(from=xlim[1],to=xlim[2], length.out=10))))
    
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
    axis(side=BELOW<-1,at=param$n_axis_ticks[[1]],labels=param$n_axis_ticks[[1]],cex.axis=param$cex)
    title(xlab = param$xlab , cex.lab = param$cex.lab)
    title(ylab = param$ylab , cex.lab = param$cex.lab)
    title(main = param$main , cex = param$cex.main)
}

.prepCS <- function(CS,states, param){
    CS$y0 <- -0.3
    CS$y1 <- -0.1
    CS$density <- param$densityCS
    col_local <- param$colCS(length(states))
    CS$col <- col_local[match(CS$stateID,states)]
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
        lwd = param$lwd,
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

.prepGR <- function(gr, types,param){
    elem<- data.frame("chr" = "chr",
            "x0" = param$xlim[1],
            "x1" = param$xlim[2],
            "width" = param$xlim[2] - param$xlim[1],
            "strand" = "*",
            "type" = "line",
            "gene_id" = "None")
    if(!is.null(nrow(gr))){
        colnames(gr) <- c("chr","x0","x1","width","strand","type","gene_id")
        gr <- rbind(elem,gr)
        
    }else {
        gr <- elem
    }
   
   
   elemY0 <- param$ylim[1]
   elemY0[gr$strand == "*"] <- param$ylim[1] + 0.35
   elemY0[gr$strand == "+"] <- param$ylim[1] + 0.4
   elemY0[gr$strand == "-"] <- param$ylim[1] + 0.15

   elemY1 <- param$ylim[1]
   elemY1[gr$strand == "*"] <- param$ylim[1] + 0.35
   elemY1[gr$strand == "+"] <- param$ylim[1] + 0.55
   elemY1[gr$strand == "-"] <- param$ylim[1] + 0.3

   gr$y0 <- elemY0
   gr$y1 <- elemY1

   ## Setting colours 
   col_local <- param$colGR(length(types))
   gr$col <- col_local[match(gr$type,types)]
   
   gr$cex <-c(param$cex.lab,(rep(gr$cex,nrow(gr)-1)*0.5))
   return(gr)
}

.getGeneStarts  <- function(gr){
    genes <- gr[gr$type == "transcript",]
    if(!is.null(nrow(genes))){
       
        genes <- split(genes, genes$gene_id)
        genes <- lapply(genes,function(x)return(min(x$x0)))
        return(genes)
    }
}

.drawGeneRef  <- function(gr){
    # getting start pos of GR 
    genes <- .getGeneStarts(gr)
    for(seg in seq_len(nrow(gr))){
        
        if(gr$type[seg] == "line"){
            
            lines(c(gr$x0[seg],gr$x1[seg]),
                c(gr$y0[seg],gr$y0[seg]),
                col = "black",
                lwd=1)
            text(gr$x1[seg] + gr$width[seg] * 0.015,
                 gr$y0[seg] + 0.2,
                "+",
                col = "black",
                cex = gr$cex[seg])
            text(gr$x1[seg] + gr$width[seg] * 0.015,
                gr$y1[seg] - 0.2,
                "-",
                col = "black",
                cex = gr$cex[seg])
        }else if(gr$type[seg] == "exon"){
             lines(c(gr$x0[seg],gr$x1[seg]),
                c(gr$y0[seg],gr$y0[seg]),
                col = gr$col[seg])
        } else if(gr$type[seg] == "transcript"){
            rect(gr$x0[seg],
                gr$y0[seg],
                gr$x1[seg],
                gr$y1[seg],
                col = gr$col[seg])
            loc <- genes[[gr$gene_id[seg]]]
            if(loc == gr$x0[seg] & gr$strand[seg] == "-"){
                text(gr$x0[seg],
                gr$y0[seg] - 0.08,
                gr$gene_id[seg],
                col = "black",
                cex = gr$cex[seg] * 0.5)
            } else if(loc == gr$x0[seg] & gr$strand[seg] == "+"){
                text(gr$x0[seg],
                gr$y1[seg] + 0.08,
                gr$gene_id[seg],
                col = "black",
                cex = gr$cex[seg] * 0.5)
            }
            
        } else {
            rect(gr$x0[seg],
                gr$y0[seg],
                gr$x1[seg],
                gr$y1[seg],
                col = gr$col[seg])
            
        }
    }
}

.reorder <- function(cs,dat = "cs"){
    if(dat == "cs"){
        states <- unique(cs$stateID)
        cols <- rep(NA,length(states))
        for(i in seq_along(states)){
            cols[i] <- unique(cs$col[cs$stateID == states[i]])
        }
        return(list("stateID" = states, "col" = cols))
    }else{
        types <- unique(cs$type)
        cols <- rep(NA,length(types))
        for(i in seq_along(types)){
            cols[i] <- unique(cs$col[cs$type == types[i]])
        }
        return(list("type" = types, "col" = cols))
    }
    
}

.prepLegend <- function(chip,
    occup,
    cs,
    gof,
    gr,
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
        leg$gof <- list("x" = param$xlim[1] -((param$xlim[2] - param$xlim[1]) * 0.02),
        "y" = 1,
        "legend" = paste(paste(names(tmp),"=",signif(tmp,4)),collapse = "\n"),
        "col" = "black",
        "fill" = NA,
        "border" = NA,
        "cex" = param$cex)
    }
    if(any(cs != "empty")){
        cs <- cs[[1]]
        cs <- .reorder(cs,dat = "cs")
        leg$cs <- list("x" = param$xlim[2] +((param$xlim[2] - param$xlim[1]) * 0.2),
        "y" = 0.5,
        "legend" = cs$stateID,
        "col" = NA,
        "fill" = cs$col,
        "border" = NA,
        "cex" = param$cex)
    }
    if(any(gr != "empty")){
        gr <- gr[[1]]
        gr <- gr[gr$type != "line",]
        gr <- .reorder(gr,dat="gr")
        leg$gr <- list("x" = param$xlim[1] -((param$xlim[2] - param$xlim[1]) * 0.2),
        "y" = -0.3,
        "legend" = gr$type,
        "col" = NA,
        xjust = 1,
        "fill" = gr$col,
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