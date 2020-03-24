## code for BoxCar processing and visualisation on MSnExp data

sp_plot <- function(x) {
    d <- as(x, "data.frame")
    if (centroided(x))
        p <- ggplot(d, aes(x = mz, xend = mz,
                           y = 0, yend = i)) +
            geom_segment()
    else
        p <- ggplot(d, aes(x = mz, y = i)) +
            geom_line()
    p
}

##' This function extracts the box start and end positions for all
##' spectra in an `MSnExp` object. The first box coordinates are not
##' returned, as they corresond to either the full spectrum (for
##' non-boxcar MS1 spectra) or the first start and last end box
##' positions for a boxcar MS1 spectrum.
##'
##' @title Extract box positions
##' @param x An `MSnExp` object
##' @param offset `numeric(1)` defining the offset to remove/add to
##'     the start/end of the box. Default is `0`, i.e. to leave them
##'     as is.
##' @param fcol The name of the feature variable containing the box
##'     data. Default is `filterString` and is expected to be
##'     `"FTMS + p NSI SIM msx ms [start_1-end_1, start_2-end_2, ..., start_n-end_n]"`
##'     where `start` and `end` are box start and end positions.
##' @md
##' @return A list of data frames of length `length(x)` with start and
##'     end box values.
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


##' Takes an experiment with a set of boxcar spectra and plots them as
##' a single spectrum, colouring the peaks based on their origin.
##'
##' @title Overlay boxcar spectra
##' @param x An `MSnExp` object containing a set of boxcar spectra.
##' @param ... Additional arguments passed to `plot`.
##' @return Used for its side effects.
##' @md
##' @author Laurent Gatto
## bc_plot <- function(x, j, ...) {
##     stopifnot(inherits(x, "MSnExp"))
##     bx <- bc_boxes(x)
##     stopifnot(all(sapply(bx, nrow) > 0))
##     d <- as(x, "data.frame")
##     d$sp <- as.numeric(as.factor(d$rt))
##     mx <- max(d$i)
##     plot(d[, c("mz", "i")], type = "h", ylim = c(0, mx), col = NA, ...)
##     grid()
##     segments(d$mz, rep(0, nrow(d)),
##              d$mz, d$i,
##              col = d$sp)
## }


##' Takes an experiment with a set of boxcar spectra and plots them as
##' a single spectrum, colouring the peaks based on their origin.
##'
##' @title Overlay boxcar spectra
##' @param x An `MSnExp` object containing a set of boxcar spectra.
##' @return Used for its side effects. Invisibly returns a `ggplot` or
##'     `plotly` object.
##' @md
##' @author Laurent Gatto
##' @aliases bc_plotly
bc_plot <- function(x) {
    stopifnot(inherits(x, "MSnExp"))
    d <- as(x, "data.frame")
    d$rt <- factor(d$rt)
    p <- ggplot(d, aes(x = mz, xend = mz,
                       y = 0, yend = i,
                       colour = rt)) +
        geom_segment(show.legend = FALSE)
    bx <- bc_boxes(x)
    bx <- do.call("rbind", bx)
    bx <- bx[order(bx[, 1]), ]
    p + annotate("rect",
                 xmin = bx[, 1],
                 xmax = bx[, 2],
                 ymin = rep(0, nrow(bx)),
                 ymax = rep(max(d$i) * 1.01, nrow(bx)),
                 alpha = 0.2,
                 colour = "black",
                 size = 0.1,
                 fill = "grey")
}

##' @rdname bc_plot
bc_plotly <- function(x) {
    p <- bc_plot(x)
    plotly::ggplotly(p)
}

##' This function uses the identification of boxcar boxes by the
##' `bc_boxes` function to determine if the spectrum is indeed a
##' boxcar spetrum.
##'
##' @title Identifies boxcar spectra
##' @param x An `MSnExp` object containing a set of boxcar spectra.
##' @param fcol The name of the feature variable containing the box
##'     data. Default is `filterString`. See the `bc_boxes` function
##'     for details.
##' @return A `logical` of length `length(x)` indicating of a spectrum
##'     is a boxcar spectrum.
##' @md
##' @author Laurent Gatto
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
##' @param x An `MSnExp` object containing a set of boxcar spectra.
##' @param fcol The name of the feature variable containing the box
##'     data. Default is `filterString`. See the `bc_boxes` function
##'     for details.
##' @return A `numeric` of length `length(x)` with the boxcar groups.
##' @author Laurent Gatto
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


##' Takes an MS experiment with boxcar spectra and sets any peak
##' outside of the boxes (+/- an optional offset) to zero. If the
##' spectrum is a fill, non-boxcar spectrum, it is left as is.
##'
##' @title Set peaks outside of boxes to zero
##' @param x x An `MSnExp` object containing a set of boxcar spectra.
##' @param offset `numeric(1)` defining the offset to remove/add to
##'     the start/end of the box. Default is `0`, i.e. to leave them
##'     as is.
##' @return An object of class `Spectra` with processed spectra.
##' @md
##' @author Laurent Gatto
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
    Spectra(l, elementMetadata = DataFrame(fData(x)))
}


##' A function to be passed to `combineSpectra` to combine boxcar
##' spectra.
##'
##' @title Combine boxcar spectra
##' @param x A list of spectra
boxcarCombine <- function(x)
    meanMzInts(x, intensityFun = base::sum, mzd = 0)


## Cast back into an MSnExp
.as_MSnExp <- function(x, fn) {
    res <- as(x, "MSnExp")
    res@processingData@files <- fn
    MSnbase:::logging(res, "boxcar processed")
}
