##' Takes an experiment with a set of boxcar spectra and plots them as
##' a single spectrum, colouring the peaks based on their origin.
##'
##' @title Overlay boxcar spectra
##' 
##' @param x An `MSnExp` object containing a set of boxcar spectra.
##' 
##' @return Used for its side effects. Invisibly returns a `ggplot` or
##'     `plotly` object.
##' 
##' @author Laurent Gatto
##' 
##' @aliases bc_plotly
##'
##' @import ggplot2
##' 
##' @export
bc_plot <- function(x) {
    stopifnot(inherits(x, "MSnExp"))
    i <- mz <- NULL
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
##'
##' @export
bc_plotly <- function(x) {
    stopifnot(require("plotly"))
    p <- bc_plot(x)
    plotly::ggplotly(p)
}


