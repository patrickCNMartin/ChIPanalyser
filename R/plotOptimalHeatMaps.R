plotOptimalHeatMaps<-function(optimalParam,parameter="all",Contour=TRUE){
    #paramter extraction
    if(parameter=="correlation"){
        parameter<-"meanCorr"
    }
    if(parameter=="MSE"){
        parameter <- "meanMSE"
    }
    if(parameter=="theta"){
        parameter <- "meanTheta"
    }
    #Extracting Data from optimal Paramter
    if(length(names(optimalParam[[2]])) > 1 & parameter == "all"){
        corrMatrix <- optimalParam[[2]][[1]]
        MSEMatrix <- optimalParam[[2]][[2]]
        thetaMatrix <- optimalParam[[2]][[3]]
        lambdas <- rownames(corrMatrix)
        boundMolecules <- colnames(corrMatrix)
        parameter <- 'all'
    }
    if(length(names(optimalParam[[2]])) > 1 & parameter != "all"){
        optimalMatrix <- optimalParam[[2]][[parameter]]
        parameter <- parameter
        lambdas <- rownames(optimalMatrix)
        boundMolecules <- colnames(optimalMatrix)
    }

    if(length(names(optimalParam[[2]])) == 1 &&
    names(optimalParam[[2]]) == parameter){
        optimalMatrix <- optimalParam[[2]]
        parameter <- optimalParam[[3]]
        lambdas <- rownames(optimalMatrix)
        boundMolecules <- colnames(optimalMatrix)
    }
    if(length(names(optimalParam[[2]])) == 1 &&
    names(optimalParam[[2]]) != parameter){
        stop(paste0(deparse(substitute(optimalParam)),
        ": Parameter Mismatch. "))
    }
    #Bulding colour palette
    heatColors <- rev(heat.colors(100))

    xlabs <- boundMolecules
    ylabs <- lambdas
    contLev <- round(min(c(length(boundMolecules), length(lambdas))))
    selectedColor <- "gray20"
    optimumColor <- "darkolivegreen3"

    corrColor <- "#00979b"
    corrMax <- "#9b0400"

    mseColor <- "#f2ae72"
    mseMin <- "#72b6f2"

    thetaColor <- "#0072B2"
    thetaMax <- "#b24000"

    colfunc <- colorRampPalette(c("white", corrColor))
    correlationColors <- colfunc(100)
    colfunc <- colorRampPalette(c(mseColor,"white"))
    MSEColors <- colfunc(100)
    colfunc <- colorRampPalette(c("white",thetaColor))
    meanThetaColor <- colfunc(100)

    ## Plotting Correlation Heat maps
    if(parameter=="meanCorr"){
        par(cex=0.6)
        par(mar=c(6,6.5,6, 1)+0.1)
        graphics::image(1:length(xlabs),1:length(ylabs),t(optimalMatrix),
            axes = FALSE,
            xlab=" ", ylab=" ",
            col=correlationColors, cex.lab=3.5)
        title(main="Correlation",cex.main=4.3)
        title( ylab="Scaling Factor", line=4.3, cex.lab=2.5)
        title(xlab="Number of Bound Molecules",line=4, cex.lab=2.5)
        for(textXId in 1:length(xlabs)){
            for(textYId in 1:length(ylabs)){
                text(textXId,textYId,signif(optimalMatrix[textYId,textXId],3),
                cex=1.2)
            }
        }
        if(Contour){
            contour(1:length(xlabs),1:length(ylabs), t(optimalMatrix),
            nlevels = contLev, add = TRUE,drawlabels=FALSE,
            col = "black", lwd=1, labcex = 1.0)
        }
        maxCorr <- which(optimalMatrix==max(optimalMatrix), arr.ind=TRUE)
        rect(maxCorr[,2]-0.5,maxCorr[,1]-0.5,maxCorr[,2]+
            0.5,maxCorr[,1] + 0.5, border=corrMax, lwd=2)
        axis(BELOW<-1, at=1:length(xlabs), labels=xlabs, cex.axis=2)
        axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,
            las= HORIZONTAL<-1,cex.axis=2)

    }


    ## Plotting MSE heat Maps
    if(parameter=="meanMSE"){
        par(cex=0.6)
        par(mar=c(6,6.5,6, 1)+0.1)
        graphics::image(1:length(xlabs),1:length(ylabs),
            t(log10(optimalMatrix)),
            axes = FALSE, xlab=" ",
            ylab=" ",col=MSEColors,cex.lab=3.5)
        title(main="Mean Squared Error",cex.main=4.3)
        title( ylab="Scaling Factor", line=4.3, cex.lab=2.5)
        title(xlab="Number of Bound Molecules",line=4, cex.lab=2.5)
        for(textXId in 1:length(xlabs)){
            for(textYId in 1:length(ylabs)){
                text(textXId,textYId,signif(optimalMatrix[textYId,textXId],3),
                cex=1.2)
            }
        }
        if(Contour){
            contour(1:length(xlabs),1:length(ylabs), t(log10(
            optimalMatrix*1000)),nlevels=contLev, drawlabels=FALSE,
            add = TRUE, col = "black", lwd=1, labcex = 1.0)
        }
        minMSE <- which(optimalMatrix==min(optimalMatrix), arr.ind=TRUE)
        rect(minMSE[,2]-0.5,minMSE[,1]-0.5,minMSE[,2]+0.5,minMSE[,1]+0.5,
            border=mseMin, lwd=2)
        axis(BELOW<-1, at=1:length(xlabs), labels=xlabs, cex.axis=2)
        axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,
        las= HORIZONTAL<-1,cex.axis=2)

    }

    ## Plotting Theta heat maps

    if(parameter=="meanTheta"){
        par(cex=0.6)
        par(mar=c(4, 4, 4, 1)+0.1)
        graphics::image(1:length(xlabs),1:length(ylabs),
            t(optimalMatrix), axes = FALSE,
            xlab="Number of bound molecules",
            ylab=expression(lambda),col=meanThetaColor)
        title(main="Optimal Parameters - Theta (Corr/MSE)", cex.main = 1.6)
        if(Contour){
            contour(1:length(xlabs),1:length(ylabs),t(optimalMatrix),
                nlevels=contLev/2,labels=NULL, add = TRUE,
                drawlabels=FALSE,
                col = "#2F4F4F", lwd=1, labcex = 1.0)
        }
        maxTheta <- which(optimalMatrix==max(optimalMatrix), arr.ind=TRUE)
        rect(maxTheta[,2]-0.5,maxTheta[,1]-0.5,maxTheta[,2]+
            0.5,maxTheta[,1]+0.5, border=thetaMax, lwd=2)
        axis(BELOW<-1, at=1:length(xlabs), labels=xlabs, cex.axis=0.7)
        axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,
            las= HORIZONTAL<-1,cex.axis=0.7)

    }


    ## Plotting all heat maps with Theta transform

    if(parameter=="all"){
        par(mfrow=c(3,1))
        par(cex=0.6)
        par(mar=c(4, 4, 4, 1)+0.1)
        ##Correlation Map
        graphics::image(1:length(xlabs),1:length(ylabs),t(corrMatrix),
            axes = FALSE,
            xlab="Number of bound molecules",
            ylab=expression(lambda),col=correlationColors)
        title(main="Correlation",cex.main=1.6)
        for(textXId in 1:length(xlabs)){
            for(textYId in 1:length(ylabs)){
                text(textXId,textYId,signif(corrMatrix[textYId,textXId],3),
                cex=0.7);
            }
        }
        if(Contour){
            contour(1:length(xlabs),1:length(ylabs), t(corrMatrix),
                nlevels = contLev, add = TRUE,
                drawlabels=FALSE,
                col = "#2F4F4F", lwd=1, labcex = 1.0)
        }
        maxCorr<-which(corrMatrix==max(corrMatrix), arr.ind=TRUE)
        rect(maxCorr[,2]-0.5,maxCorr[,1]-0.5,maxCorr[,2]+
            0.5,maxCorr[,1]+0.5, border=corrMax, lwd=2)
        axis(BELOW<-1, at=1:length(xlabs), labels=xlabs, cex.axis=0.7)
        axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,
            las= HORIZONTAL<-1,cex.axis=0.7)

        ## MSE Map
        graphics::image(1:length(xlabs),1:length(ylabs),
            t(log10(MSEMatrix)),
            axes = FALSE,
            xlab="Number of bound molecules",
            ylab=expression(lambda),col=MSEColors)
        title(main="Mean Squared Error",cex.main=1.6)
        for(textXId in 1:length(xlabs)){
            for(textYId in 1:length(ylabs)){
                text(textXId,textYId,signif(MSEMatrix[textYId,textXId],3),
                cex=0.7)
            }
        }
        if(Contour){
            contour(1:length(xlabs),1:length(ylabs), t(log10(MSEMatrix*1000)),
                nlevels=contLev,labels=NULL,
                drawlabels=FALSE,
                add = TRUE, col = "#2F4F4F", lwd=1, labcex = 1.0)
        }

        minMSE<-which(MSEMatrix==min(MSEMatrix), arr.ind=TRUE)
        rect(minMSE[,2]-0.5,minMSE[,1]-0.5,minMSE[,2]+
            0.5,minMSE[,1]+0.5, border=mseMin, lwd=2)
        axis(BELOW<-1, at=1:length(xlabs), labels=xlabs, cex.axis=0.7)
        axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,
            las= HORIZONTAL<-1,cex.axis=0.7)


        ## Theta Map
        graphics::image(1:length(xlabs),1:length(ylabs),t(thetaMatrix),
            axes = FALSE,
            xlab="Number of bound molecules",
            ylab=expression(lambda),col=meanThetaColor)
        title(main="Optimal Parameters - Theta (Corr/MSE)", cex.main = 1.6)


        if(Contour){
            contour(1:length(xlabs),1:length(ylabs),t(thetaMatrix),
            nlevels= contLev/2 ,drawlabels=FALSE,
            labels=NULL, add = TRUE, col = "#2F4F4F", lwd=1, labcex = 1.0)
        }
        maxTheta<-which(thetaMatrix==max(thetaMatrix), arr.ind=TRUE)
        rect(maxTheta[,2]-0.5,maxTheta[,1]-0.5,maxTheta[,2]+
            0.5,maxTheta[,1]+0.5, border=thetaMax, lwd=2)
        axis(BELOW<-1, at=1:length(xlabs), labels=xlabs, cex.axis=0.7)
        axis(LEFT <-2, at=1:length(ylabs), labels=ylabs, las= HORIZONTAL<-1,
            cex.axis=0.7)

    }

}
