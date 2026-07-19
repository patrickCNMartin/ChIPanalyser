#######################################################################
# Functions
#######################################################################


computeOccupancy <- function(genomicProfiles,
                             parameterOptions = NULL, norm = TRUE, verbose = TRUE) {
    # Validity checking
    .validateGenomicProfiles(genomicProfiles, parameterOptions)

    # Generating parameterOptions if not given by user.
    # All Value will be defaut settings
    if (!is.null(parameterOptions)) {
        genomicProfiles <- .updateGenomicProfiles(genomicProfiles, parameterOptions)
    }
    # GenomicProfiles parameter Extraction

    PWMGRList <- GRangesList(profiles(genomicProfiles))
    DNALength <- as.numeric(DNASequenceLength(genomicProfiles))
    averageExpPWMScore <- averageExpPWMScore(genomicProfiles)
    lambda <- lambdaPWM(genomicProfiles)

    ploidy <- as.numeric(ploidy(genomicProfiles))
    boundMolecules <- boundMolecules(genomicProfiles)
    backgroundSignal <- backgroundSignal(genomicProfiles)
    maxSignal <- maxSignal(genomicProfiles)


    # Computing Occupancy at sites higher than threshold
    MultiParam <- vector("list", (length(lambda) * length(boundMolecules)))
    PWMScore <- vector("list", length(PWMGRList))

    ### Essentially making sure that if regions dont have accessible DNA
    ## They will still be kept it will make things easier later down the line
    ## Also you just return a big fat falt line anyway

    dropLoci <- drop(genomicProfiles)

    if (length(grep("No loci dropped", dropLoci)) == 0) {
        widthDisplay <- round(options()$width * 0.5)
        message(
            "No Accessible DNA in: ", paste(rep(" ",
                times = (widthDisplay - nchar("StepSize: ") - nchar(dropLoci[1]))
            ), collapse = ""),
            dropLoci, "\n", "\n"
        )
        chr <- vapply(strsplit(dropLoci, ":"), "[[", character(1), 1)
        start <- vapply(strsplit(vapply(strsplit(dropLoci, ":"), "[[", character(1), 2), "\\.."), "[[", character(1), 1)
        end <- vapply(strsplit(vapply(strsplit(dropLoci, ":"), "[[", character(1), 2), "\\.."), "[[", character(1), 2)
        LocalDrop <- GRanges(
            seqnames = chr, ranges = IRanges(as.numeric(start), as.numeric(end)),
            PWMScore = rep(0, length(chr)),
            DNAaffinity = rep(0, length(chr)),
            Occupancy = rep(0, length(chr))
        )
        name <- c(rep(names(PWMGRList), times = vapply(PWMGRList, length, integer(1))), dropLoci)
    } else {
        name <- rep(names(PWMGRList), times = vapply(PWMGRList, length, integer(1)))
    }


    # chromatinState if not continuous Data

    Occupancy <- vector("list", length(PWMGRList))
    # names(Occupancy) <- names(PWMGRList)
    result <- list()
    emptyGR <- GRanges()
    # Progress message when required
    if (verbose) {
        message("Computing Occupancy at sites higher than threshold. \n")
    }

    counter <- 0
    ParaVal <- c()

    # Computing Occupancy
    for (k in seq_along(lambda)) {
        for (j in seq_along(boundMolecules)) {
            PWMScore <- unlist(PWMGRList)$PWMScore

            Access <- unlist(PWMGRList)$DNAaffinity

            Occupancy <- rep(0, length(PWMScore))
            Occupancy <- (Access * boundMolecules[j] * exp((1 / lambda[k]) * PWMScore)) /
                (Access * boundMolecules[j] * exp((1 / lambda[k]) * PWMScore) +
                    DNALength * ploidy * averageExpPWMScore[k])

            Occupancy <- backgroundSignal + Occupancy * (maxSignal - backgroundSignal)
            # Normalising Ocupancy Signal
            if (norm == TRUE) {
                maxOccupancy <- max(c(maxSignal, max(Occupancy)))
                Occupancy <- Occupancy / maxOccupancy
                ZeroBackground <- backgroundSignal / maxOccupancy
            } else {
                ZeroBackground <- backgroundSignal
            }


            buffer <- unlist(PWMGRList)

            buffer$Occupancy <- Occupancy
            if (length(grep("No loci dropped", dropLoci)) == 0) {
                buffer <- c(buffer, LocalDrop)
            }
            # Extracting and pasting names of Parameters
            names(buffer) <- name
            if (length(name > 1)) {
                level <- factor(names(buffer), levels = unique(name))
                result <- split(buffer, level)
            } else {
                result <- buffer
            }

            counter <- counter + 1
            MultiParam[[counter]] <- result
            ParaVal <- c(ParaVal, paste0(
                "lambda = ", lambda[k], " & ",
                "boundMolecules = ", boundMolecules[j]
            ))
        }
    }

    # ReBuilding GenomicProfileParameters
    names(MultiParam) <- ParaVal

    .ZeroBackground(genomicProfiles) <- ZeroBackground
    .profiles(genomicProfiles) <- MultiParam
    .tags(genomicProfiles) <- "Occupancy"
    .paramTag(genomicProfiles) <- "Occupancy"
    return(genomicProfiles)
}
