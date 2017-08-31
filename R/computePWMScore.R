#######################################################################
# Functions
#######################################################################

computePWMScore <- function(DNASequenceSet,genomicProfileParameters,
    setSequence = NULL,DNAAccessibility = NULL,verbose = TRUE){

    # Validity checking
    if(!.is.genomicProfileParameter(genomicProfileParameters)){
    stop(paste0(deparse(substitute(genomicProfileParameters)),
    " is not a Genomic Profile Paramter object.
    Please set a Genomic Profile Parameters Object."))
    }
    if(!.is.genomeWideComputed(genomicProfileParameters)){
    stop(paste0("Genome Wide PWM Score has not yet been computed in ",
    deparse(substitute(genomicProfileParameters))))
    }
    if(class(DNASequenceSet) != "BSgenome" &
    class(DNASequenceSet)!="DNAStringSet"){
    stop(paste0(deparse(substitute(DNASequenceSet)),
    " is not a BSgenome Object or a DNAStringSet"))
    }
    if(class(DNASequenceSet) == "BSgenome"){
    DNASequenceSet<-getSeq(DNASequenceSet)
    }
    if(class(DNAAccessibility) != "GRanges" & !is.null(DNAAccessibility)){
    stop(paste0(deparse(substitute(DNAAccessibility)),
    " must be a GRanges Object."))
    }
    if(class(setSequence) != "GRanges" & !is.null(setSequence)){
    stop(paste0(deparse(substitute(setSequence)),
    " must be a GRanges Object."))
    }

    #Extracting genomic Profile Parameters
    PWMThreshold <- PWMThreshold(genomicProfileParameters)
    minPWMScore <- minPWMScore(genomicProfileParameters)
    maxPWMScore <- maxPWMScore(genomicProfileParameters)
    PWM <- PositionWeightMatrix(genomicProfileParameters)
    strand <- whichstrand(genomicProfileParameters)
    strandRule <- strandRule(genomicProfileParameters)

    #Calculating Threshhold Value
    PWMThresholdLocal <- minPWMScore + PWMThreshold*(maxPWMScore-minPWMScore)



    # Creating a Granges Object from full sequence if setsequence is NULL
    if (is.null(setSequence)){
        setSequence <- GRanges(seqnames = names(DNASequenceSet),
            ranges = IRanges::IRanges(start = 1,
                end = sapply(DNASequenceSet, function(x) length(x))),
            strand = strand
        )
        names(setSequence)<-seqnames(setSequence)
    }

    #Scoring Sequence with PWM
    AccessibleSequence <- vector("list", length(setSequence))
    IntersectSequence <- vector("list", length(setSequence))
    NoAccess<-c()

    for(i in seq_along(setSequence)){
        # Refining SequenceSet if DNA accessibility Data is available
        if(!is.null(DNAAccessibility)){
            if(verbose==TRUE & i == 1){
            message("Processing DNA Acccssibility \n",
            "Extracting Sites Above threshold")
            }

            DNAAccessibility <- DNAAccessibility[as.character(
                seqnames(DNAAccessibility)) %in%
                as.character(seqnames(setSequence))]

            IntersectSequence[[i]] <- intersect(setSequence[i],
                DNAAccessibility)
            AccessibleSequence[[i]] <- DNAStringSet(DNASequenceSet[[
                which(names(DNASequenceSet)==(as.character(seqnames(
                    setSequence[i]))))]],
                start = start(IntersectSequence[[i]]),
                end = end(IntersectSequence[[i]])+ncol(PWM)-1)

    # Generating vector with Loci in setSequence that do not display
    # any accesible DNA (with DNA Access)
        if(length(AccessibleSequence[[i]])==0){
        NoAccess<-c(NoAccess," ",as.character(names(setSequence[i])))
        next
        } else {
        NoAccess<-c(NoAccess,"-")
        }

    # Score Loci of interest with PWM (with DNA Accessibility)
        AccessibleSequence[[i]] <- .scoreDNAStringSet(PWM,
            AccessibleSequence[[i]],
            strand = strand,strandRule = strandRule)

    # Extract regions within setSequence
    # that have a PWM score higher than Threshold (with DNA Accessibility)
        indexPWMThresholded <- .getIndexOfPWMThresholded(
            AccessibleSequence[[i]][[1]],PWMThresholdLocal)

    # Building GRanges with sites above Threshold (with DNA Accessibility)
        strandLocal <- vector("list",length(AccessibleSequence[[i]][[1]]))
        AllSites <- GRangesList()
    #Extracting Strand Information for sites above Threshold
        for(j in seq_along(AccessibleSequence[[i]][[1]])){
            if(strand == "+-" | strand == "-+"){
                strandLocal[[j]]<-rep("*",
                    length(AccessibleSequence[[i]][[1]][[j]]))
                strandLocal[[j]][AccessibleSequence[[i]][[2]][[j]]] <- "+"
                strandLocal[[j]][AccessibleSequence[[i]][[3]][[j]]] <- "-"
            }
            if(strand == "+"){
                strandLocal[[j]] <- rep("+",
                    length(AccessibleSequence[[i]][[1]][[j]]))
            }
            if(strand == "-"){
                strandLocal[[j]] <- rep("-",
                    length(AccessibleSequence[[i]][[1]][[j]]))
            }

            GRLocal <- GRanges(seqnames = S4Vectors::Rle(unique(as.character(
                seqnames(IntersectSequence[[i]]))),
                length(indexPWMThresholded[[j]])),
                ranges = IRanges::IRanges(start = (
                    start(IntersectSequence[[i]][j]) +
                    indexPWMThresholded[[j]]-1),
                    end = (
                    start(IntersectSequence[[i]][j]) +
                    indexPWMThresholded[[j]]-1+ncol(PWM)-1)),
                strand = strandLocal[[j]][indexPWMThresholded[[j]]],
                PWMScore=
                AccessibleSequence[[i]][[1]][[j]][indexPWMThresholded[[j]]])

            AllSites <- c(AllSites,GRangesList(GRLocal))
        }
        AccessibleSequence[[i]] <- AllSites
    } else {
        if(verbose == TRUE & i == 1){
        message( "Extracting Sites Above threshold")
        }
    # setSequence if no DNA Accesibility
        AccessibleSequence[[i]] <- DNAStringSet(DNASequenceSet[[
            which(names(DNASequenceSet)==(as.character(seqnames(
                setSequence[i]))))]],
            start = start(setSequence[i]),
            end = end(setSequence[i])+ncol(PWM)-1)

    # Generating vector with Loci in setSequence that do not display
    # any accesible DNA (with DNA Access)
        NoAccess<-c("-")

    # Scoring Loci of Interest with PWM (without DNA Accesibility)
        AccessibleSequence[[i]] <- .scoreDNAStringSet(PWM,
            AccessibleSequence[[i]], strand=strand,strandRule=strandRule)

    #Extracting sites above threshold withing Loci of interest
    #(Without DNA accessibility)
        indexPWMThresholded <- .getIndexOfPWMThresholded(
            AccessibleSequence[[i]][[1]],PWMThresholdLocal)
    #Building GRanges of sites above threshold (without DNA accessibility)
        strandLocal <- vector("list",length(AccessibleSequence[[i]][[1]]))
        AllSites <- GRangesList()
        for(j in seq_along(AccessibleSequence[[i]][[1]])){
            #Extracting strand Information for sites above threshold
            if(strand == "+-" | strand == "-+"){
                strandLocal[[j]] <- rep("*",
                    length(AccessibleSequence[[i]][[1]][[j]]))
                    strandLocal[[j]][AccessibleSequence[[i]][[2]][[j]]] <- "+"
                    strandLocal[[j]][AccessibleSequence[[i]][[3]][[j]]] <- "-"
            }
            if(strand == "+"){
                strandLocal[[j]] <- rep("+",
                    length(AccessibleSequence[[i]][[1]][[j]]))
            }
            if(strand == "-"){
                strandLocal[[j]] <- rep("-",
                    length(AccessibleSequence[[i]][[1]][[j]]))
            }

            GRLocal <- GRanges(seqnames = S4Vectors::Rle(unique(as.character(
                seqnames(setSequence[i]))),
                    length(indexPWMThresholded[[j]])),
                ranges = IRanges::IRanges(start=(
                    start(setSequence[i][j]) +
                    indexPWMThresholded[[j]]-1),
                    end = (start(setSequence[i][j])+indexPWMThresholded[[j]]-
                    1 + ncol(PWM)-1)),
                strand = strandLocal[[j]][indexPWMThresholded[[j]]],
                PWMScore=
                AccessibleSequence[[i]][[1]][[j]][indexPWMThresholded[[j]]])
            AllSites <- c(AllSites,GRangesList(GRLocal))
        }
        AccessibleSequence[[i]] <- AllSites
        }
    }

    AccessibleSequence <- lapply(AccessibleSequence,unlist)
    names(AccessibleSequence) <- names(setSequence)
    if(length(grep("-",NoAccess)) == 0){
    warning(NoAccess," do not contain any accessible sites ")
    }

    genomicProfileParameters <-.AllSitesAboveThresholdReplace(
        genomicProfileParameters,GRangesList(
        AccessibleSequence[which(lapply(AccessibleSequence,length) != 0)]))

    genomicProfileParameters <-.NoAccessReplace(genomicProfileParameters,
        NoAccess)
    return(genomicProfileParameters)

}
