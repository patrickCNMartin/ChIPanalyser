#######################################################################
# Functions profile Accuracy
#######################################################################


profileAccuracyEstimate <- function(
  genomicProfiles, ChIPScore,
  parameterOptions = NULL, method = "all", cores = 1
) {
    # Validity checking
    .validateGenomicProfiles(genomicProfiles, parameterOptions)
    if (!inherits(ChIPScore, "ChIPScore") &
        !is.null(ChIPScore)) {
        stop(
            deparse(substitute(ChIPScore)),
            " is not a ChIPScore Object."
        )
    }

    if (!is.null(parameterOptions)) {
        genomicProfiles <- .updateGenomicProfiles(genomicProfiles, parameterOptions)
    }

    # If Sites are not accesible or have no overlapp with chipprofile
    dropLoci <- drop(genomicProfiles)
    if (length(grep("No loci dropped", dropLoci)) == 0) {
        widthDisplay <- round(options()$width * 0.5)
        message(
            "No Accessible DNA in: ", paste(rep(" ",
                times = (widthDisplay - nchar("StepSize: ") - nchar(dropLoci[1]))
            ), collapse = ""),
            dropLoci, "\n", "\n"
        )
    }

    ## Extract and reorder
    predictionSet <- profiles(genomicProfiles)

    ValidationSet <- scores(ChIPScore)

    stepSize <- stepSize(genomicProfiles)

    GoF <- .cleanGoF(predictionSet, ValidationSet, stepSize)
    ## paralle : choosing which one to split over

    if (length(predictionSet) > length(ValidationSet)) {
        GoF <- BiocParallel::bplapply(GoF, .GoFPred, method,
            step = stepSize, BPPARAM = .bpParam(cores)
        )
    } else {
        GoF <- lapply(GoF, .GoFLoci, method = method, step = stepSize, cores = cores)
    }

    ### computing mean over all regions
    GoF <- lapply(GoF, .meanGoFScore)

    .profiles(genomicProfiles) <- GoF
    .tags(genomicProfiles) <- "GoF"

    return(genomicProfiles)
}
