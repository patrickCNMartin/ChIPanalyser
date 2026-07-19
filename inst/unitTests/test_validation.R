## Tests for the shared .validateGenomicProfiles() argument check.
##
## The point of these tests is not that invalid input errors -- it is that the
## error message names the variable *the user passed*, not the formal argument
## name of the internal helper. That property depends on a frame walk in
## .validateGenomicProfiles() and breaks silently if the helper is refactored
## into a nested function, so it is worth pinning down.

library("ChIPanalyser")

.validGP <- function() {
    genomicProfiles(PFM = matrix(1, 4, 6), PFMFormat = "matrix")
}

## Error message reports the caller's variable name for genomicProfiles
test_validate_reportsCallerNameForGenomicProfiles <- function() {
    myProfile <- "not a genomicProfiles object"
    for (f in list(
        function() computeOptimal(myProfile),
        function() computeOccupancy(myProfile),
        function() profileAccuracyEstimate(myProfile)
    )) {
        msg <- tryCatch(f(), error = conditionMessage)
        checkTrue(
            grepl("myProfile", msg, fixed = TRUE),
            paste("error should name the caller's variable, got:", msg)
        )
    }
}

## ... and for parameterOptions, which is the second frame-walked argument
test_validate_reportsCallerNameForParameterOptions <- function() {
    myOpts <- "not a parameterOptions object"
    gp <- .validGP()
    msg <- tryCatch(computeOptimal(gp, parameterOptions = myOpts),
        error = conditionMessage
    )
    checkTrue(
        grepl("myOpts", msg, fixed = TRUE),
        paste("error should name the caller's variable, got:", msg)
    )
}

## Non-symbol arguments deparse to their expression, not the formal name
test_validate_reportsCallerExpression <- function() {
    lst <- list(gp = "not a genomicProfiles object")
    msg <- tryCatch(computeOptimal(lst$gp), error = conditionMessage)
    checkTrue(
        grepl("lst$gp", msg, fixed = TRUE),
        paste("error should name the caller's expression, got:", msg)
    )
}

## parameterOptions is optional package-wide: NULL must never raise the
## parameterOptions error (it should fall through to later validation)
test_validate_acceptsNullParameterOptions <- function() {
    gp <- .validGP()
    msg <- tryCatch(computeOptimal(gp, parameterOptions = NULL),
        error = conditionMessage
    )
    checkTrue(
        !grepl("is not a parameterOptions object", msg, fixed = TRUE),
        paste("NULL parameterOptions must be accepted, got:", msg)
    )
}