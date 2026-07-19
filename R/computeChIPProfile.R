#######################################################################
# Functions
#######################################################################


computeChIPProfile <- function(
  genomicProfiles, loci = NULL,
  parameterOptions = NULL,
  norm = TRUE, method = c("moving_kernel", "truncated_kernel", "exact"),
  peakSignificantThreshold = NULL, cores = 1, verbose = TRUE
) {
    # Validity checking
    if (!is.null(loci)) {
        if (is(loci, "ChIPScore")) {
            chipMean(genomicProfiles) <- chipMean(loci)
            chipSd(genomicProfiles) <- chipSd(loci)
            chipSmooth(genomicProfiles) <- chipSmooth(loci)
            maxSignal(genomicProfiles) <- maxSignal(loci)
            backgroundSignal(genomicProfiles) <- backgroundSignal(loci)
            loci <- loci(loci)
        } else if (is(loci, "Granges")) {
            loci <- loci
        }
    } else {
        stop(
            "Please provide a set of loci to be analysed \n",
            "If you have used processingChIP, you can parse that object \n",
            "Otherwise, you can provide your own loci as a GRanges. "
        )
    }
    .validateGenomicProfiles(genomicProfiles, parameterOptions)

    if (!is.null(parameterOptions)) {
        genomicProfiles <- .updateGenomicProfiles(genomicProfiles, parameterOptions)
    }

    # Extraction of Ocuupancy and associated values
    Occup <- profiles(genomicProfiles)
    ZeroBackground <- .ZeroBackground(genomicProfiles)

    # Extraction of Occupancy Profile Parameters
    stepSize <- stepSize(genomicProfiles)
    backgroundSignal <- backgroundSignal(genomicProfiles)
    removeBackground <- removeBackground(genomicProfiles)
    chipMean <- chipMean(genomicProfiles)
    chipSd <- chipSd(genomicProfiles)
    dropLoci <- drop(genomicProfiles)
    if (!is.null(chipSmooth(genomicProfiles))) {
        chipSmooth <- chipSmooth(genomicProfiles)
    } else {
        chipSmooth <- NULL
    }
    maxSignal <- maxSignal(genomicProfiles)


    # Extracting names of sequences with no accesible DNA

    if (length(grep("No loci dropped", dropLoci)) == 0) {
        widthDisplay <- round(options()$width * 0.5)
        message(
            "No Accessible DNA in: ", paste(rep(" ",
                times = (widthDisplay - nchar("StepSize: ") - nchar(dropLoci[1]))
            ), collapse = ""),
            dropLoci, "\n", "\n"
        )
    }

    # loci fragmentation
    # Spliting loci for Chip PRofile computing
    if (!is.null(loci) & is.null(names(loci))) {
        names(loci) <- paste0(seqnames(loci), ":",
            start(ranges(loci)), "..", end(ranges(loci)),
            sep = ""
        )
    }


    LocalSet <- split(loci, names(loci))
    LocalSet <- LocalSet[match(names(Occup[[1]]), names(LocalSet))]


    # Computing Chip like profile
    if (verbose) {
        message("Computing ChIP Profile \n")
    }

    SplitGRList <- BiocParallel::bplapply(LocalSet, .internalChIPLociSplit,
        stepSize,
        BPPARAM = .bpParam(cores)
    )

    names(SplitGRList) <- names(LocalSet)


    ## method set
    if (length(Occup) > length(Occup[[1]])) {
        OccupancyVals <- BiocParallel::bplapply(Occup,
            .internalChIPOccupValsParam,
            BPPARAM = .bpParam(cores)
        )

        profile <- BiocParallel::bpmapply(.internalChIPParam,
            Occup = Occup,
            OccupancyVals = OccupancyVals,
            MoreArgs = list(
                SplitGRList = SplitGRList, LocalSet = LocalSet,
                chipMean = chipMean, chipSd = chipSd,
                stepSize = stepSize, norm = norm, chipSmooth = chipSmooth,
                peakSignificantThreshold = peakSignificantThreshold,
                ZeroBackground = ZeroBackground,
                removeBackground = removeBackground, method = method
            ),
            BPPARAM = .bpParam(cores), SIMPLIFY = FALSE
        )
    } else {
        profile <- vector("list", length(Occup))
        OccupancyVals <- vector("list", length(Occup))
        for (i in seq_along(Occup)) {
            OccupancyVals[[i]] <- BiocParallel::bplapply(Occup[[i]],
                .internalChIPOccupValsLoci,
                BPPARAM = .bpParam(cores)
            )

            profile[[i]] <- SplitGRList
            profile[[i]] <- BiocParallel::bpmapply(.internalChIPLoci,
                profile = profile[[i]], Occup = Occup[[i]], LocalSet = LocalSet,
                OccupancyVals = OccupancyVals[[i]],
                MoreArgs = list(
                    chipMean = chipMean, chipSd = chipSd,
                    stepSize = stepSize, norm = norm, chipSmooth = chipSmooth,
                    peakSignificantThreshold = peakSignificantThreshold,
                    ZeroBackground = ZeroBackground,
                    removeBackground = removeBackground, method = method
                ), BPPARAM = .bpParam(cores)
            )
        }
    }
    ### GRlist output

    profile <- lapply(profile, GRangesList)
    names(profile) <- names(Occup)

    .profiles(genomicProfiles) <- profile
    .tags(genomicProfiles) <- "ChIPProfile"
    return(genomicProfiles)
}
