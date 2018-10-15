plotOptimalHeatMaps<-function(optimalParam,contour=TRUE,col=NULL,main=NULL,layout=TRUE){

    # Parsing matricies if everthing is provided
    if(all(names(optimalParam)%in%c("Optimal Parameters","Optimal Matrix","method"))){
        optimalParam<-optimalParam[[2]]
    }
    ## replacing NA or NAN by 0 because it messes things up for the plotting
    ## It be annoying

    nans<-lapply(optimalParam, function(x){return(which(is.na(x), arr.ind=T))})


    ## screw this
    for(i in seq_along(optimalParam)){
        optimalParam[[i]][nans[[i]]]<-min(optimalParam[[i]],na.rm=TRUE)
    }


    #Setting up some paramters
    if(class(optimalParam)!="list"){
        optimalParam[]
        if(is.null(main)){
            mainTitle<-"Optimal Paramters"
        } else{
            mainTitle<-main
        }

        ylabs<-as.numeric(rownames(optimalParam))
        xlabs<-as.numeric(colnames(optimalParam))

        ifelse(!is.null(col),cols<-col,cols<-"#00979b")

        optimalParam<-list(optimalParam)
    } else {
        if(is.null(main) & length(main)>0){
            mainTitle<-names(optimalParam)
        } else{
            mainTitle<-names(optimalParam)
            mainTitle[seq_along(main)]<-main
        }
        ylabs<-as.numeric(rownames(optimalParam[[1]]))
        xlabs<-as.numeric(colnames(optimalParam[[1]]))
        ifelse(is.null(col),cols<-rainbow(length(optimalParam)),cols<-col)
        if(length(cols)!=length(optimalParam)){
            cols<-rep(cols,ceiling(length(optimalParam)/length(cols)))
        }
    }

    ## plotting
    if(layout){
        layout(matrix(1:2,ncol=2,byrow=T), width = c(6,1),height = c(1,1))
    }

    for(i in seq_along(optimalParam)){
        #if(grepl("MSE",mainTitle[i])|grepl("geometric",mainTitle[i])|grepl("ks", mainTitle[i])){
            #colfunc<-colorRampPalette(c(cols[i],"white"))
        #}else{
            colfunc<-colorRampPalette(c("white",cols[i]))
      #  }

        Colors <- colfunc(20)
        legend_image<-as.raster(matrix(rev(Colors),ncol=1))
        par(mar=c(5.5,5.5,4.5, 0.5)+0.1)

        graphics::image(1:length(xlabs),1:length(ylabs),t(optimalParam[[i]]),
            axes = FALSE, xlab=" ", ylab=" ",col=Colors)

        title(main=mainTitle[i],cex.main=1.8)
        title( ylab="Scaling Factor", line=3.5, cex.lab=1.2)
        title(xlab="Number of Bound Molecules",line=4.5, cex.lab=1.2)
        axis(1,at=seq_along(xlabs),labels=F)
        text(seq_along(xlabs),y =(-0.2), srt = 45, adj = 1,labels = xlabs, xpd = TRUE,cex=1.15)
        axis(LEFT <-2, at=1:length(ylabs), labels=ylabs,las= HORIZONTAL<-1,cex.axis=1.15)
        if(contour){
          contour(1:length(xlabs),1:length(ylabs), t(log10(
          optimalParam[[i]]*1000)),nlevels=7, drawlabels=FALSE,
          add = TRUE, col = "black", lwd=1, labcex = 1.0)
        }

        # raster scacle
        par(mar=c(3,0.5,3.2,0.5))
        plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
        text(x=1.6, y =seq(0,1,l=5) , labels = round(seq(min(optimalParam[[i]]),max(optimalParam[[i]]),l=5),2),cex=0.95)
        rasterImage(legend_image, 0, 0, 1,1)
    }
}
