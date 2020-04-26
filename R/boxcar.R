##' This function extracts the box start and end positions for all
##' spectra in an `MSnExp` object. The first box coordinates are not
##' returned, as they corresond to either the full spectrum (for
##' non-boxcar MS1 spectra) or the first start and last end box
##' positions for a boxcar MS1 spectrum.
##'
##' @title Extract box positions
##' 
##' @param x An `MSnExp` object
##' 
##' @param offset `numeric(1)` defining the offset to remove/add to
##'     the start/end of the box. Default is `0`, i.e. to leave them
##'     as is.
##' 
##' @param fcol The name of the feature variable containing the box
##'     data. Default is `filterString` and is expected to be
##'     `"FTMS + p NSI SIM msx ms [start_1-end_1, start_2-end_2, ..., start_n-end_n]"`
##'     where `start` and `end` are box start and end positions.
##' 
##' @return A list of data frames of length `length(x)` with start and
##'     end box values.
##'
##' @export
##'
##' @import MSnbase
bc_boxes <- function(x,
                     offset = 0L,
                     fcol = "filterString") {
    stopifnot(inherits(x, "MSnExp"))
    stopifnot(fcol %in% fvarLabels(x))
    if (length(offset) != 1)
        warning("Offset must be a `numeric(1)`. Ignoring other values.")
        offset <- offset[1]
    lapply(seq_len(length(x)),
           function(i) {
               bx <- as.character(fData(x)[i, fcol])
               bx <- gsub(".+\\[(.+)\\]", "\\1", bx)
               bx <- strsplit(strsplit(bx, ", ")[[1]], "-")
               bx <- t(sapply(bx, as.numeric))
               if (offset) {
                   bx[, 1] <- bx[, 1] - offset
                   bx[, 2] <- bx[, 2] + offset
               }
               bx[-1, ]
           })}


##' This function uses the identification of boxcar boxes by the
##' `bc_boxes` function to determine if the spectrum is indeed a
##' boxcar spetrum.
##'
##' @title Identifies boxcar spectra
##' 
##' @param x An `MSnExp` object containing a set of boxcar spectra.
##' 
##' @param fcol The name of the feature variable containing the box
##'     data. Default is `filterString`. See the `bc_boxes` function
##'     for details.
##' 
##' @return A `logical` of length `length(x)` indicating of a spectrum
##'     is a boxcar spectrum.
##' 
##' @author Laurent Gatto
##'
##' @export
bc_is_boxcar <- function(x, fcol = "filterString") {
    stopifnot(inherits(x, "MSnExp"))
    stopifnot(fcol %in% fvarLabels(x))
    bx <- bc_boxes(x, fcol = fcol)
    as.logical(sapply(bx, nrow))
}

##' This function takes an `MSnExp` as input and defines the boxcar
##' groups, i.e. the set of boxcar spectra that match the same full
##' MS1 spectrum. Grouping for full (non-boxcar) spectra are set the
##' `NA`. Note that the presence of full (non-boxcar) spectra is
##' required to define the groups, and it is assumed that boxcar
##' spectra following a full MS1 spectrum below to the same group. The
##' output of this function is best used to create a new feature
##' variable column used when merging boxcar spectra.
##'
##' TODO: make this an endomorphic function that return an MSnExp
##' object
##' 
##' @title Define boxcar groups
##' 
##' @param x An `MSnExp` object containing a set of boxcar spectra.
##' 
##' @param fcol The name of the feature variable containing the box
##'     data. Default is `filterString`. See the `bc_boxes` function
##'     for details.
##' 
##' @return A `numeric` of length `length(x)` with the boxcar groups.
##' 
##' @author Laurent Gatto
##'
##' @export
bc_groups <- function(x, fcol = "filterString") {
    bx <- bc_is_boxcar(x, fcol = fcol)
    if (all(bx))
        stop("No non-boxcar spectra found. Can't assign boxcar groups.")
    grps <- rep(NA_integer_, length(bx))
    k <- 0
    for (i in seq_along(bx)) {
        if (!bx[i]) k <- k + 1
        else grps[i] <- k
    }
    grps
}


##' A function to be passed to `combineSpectra` to combine boxcar
##' spectra. This function simply sums intensities and is provided as
##' an example.
##'
##' @title Combine boxcar spectra
##' 
##' @param x A list of spectra
##'
##' @export
boxcarCombine <- function(x)
    meanMzInts(x, intensityFun = base::sum, mzd = 0)


## Cast back into an MSnExp
.as_MSnExp <- function(x) {
    res <- as(x, "MSnExp")
    MSnbase:::logging(res, "boxcar processed")
}

##' Takes an MS experiment with boxcar spectra and sets any peak
##' outside of the boxes (+/- an optional offset) to zero. If the
##' spectrum is a fill, non-boxcar spectrum, it is left as is.
##'
##' @title Set peaks outside of boxes to zero
##' 
##' @param x x An `MSnExp` object containing a set of boxcar spectra.
##' 
##' @param offset `numeric(1)` defining the offset to remove/add to
##'     the start/end of the box. Default is `0`, i.e. to leave them
##'     as is.
##' 
##' @return An object of class `Spectra` with processed spectra.
##' 
##' @author Laurent Gatto
##'
##' @importFrom S4Vectors DataFrame
##'
##' @export
bc_zero_out_box <- function(x, offset = 0L) {
    stopifnot(inherits(x, "MSnExp"))
    bx <- bc_boxes(x, offset)
    mzs <- mz(x)
    ints <- intensity(x)
    .rt <- rtime(x)
    .centroided <- centroided(x)
    .smoothed <- smoothed(x)
    .fromFile <- fromFile(x)
    .polarity <- polarity(x)
    l <- lapply(seq_along(ints),
                function(i) {
                    ## Only do the following if there are any boxes,
                    ## i.e. it's a boxcar msx scan. Otherwise, the
                    ## output of bc_boxes() is an empty matrix and we
                    ## return the intensities as they are.
                    if (nrow(bx[[i]])) {
                        itvl <- findInterval(mzs[[i]], as.numeric(t(bx[[i]])))
                        ## keep odd intervals
                        out_box <- itvl %% 2 == 0
                        ints[[i]][out_box] <- 0
                    }
                    MSnbase:::Spectrum1(mz = mzs[[i]],
                                        intensity = ints[[i]],
                                        rt = .rt[[i]],
                                        centroided = .centroided[[i]],
                                        smoothed = .smoothed[[i]],
                                        fromFile = .fromFile[[i]],
                                        polarity = .polarity[[i]])
                })
    names(l) <- featureNames(x)
    res <- Spectra(l, elementMetadata = DataFrame(fData(x)))
    .as_MSnExp(res)
}
