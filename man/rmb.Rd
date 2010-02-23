\name{rmb}
\Rdversion{1.1}
\alias{rmb}
\concept{Multiple Barchart for relative frequencies and generalized Spineplots}
\title{
  Multiple Barchart for relative frequencies and generalized Spineplots
}
\description{
  The \code{rmb} function basically produces a Multiple Barchart for the relative frequencies of some target categories within each combination of the explanatory variables.
  The weights corresponding to each of these combinations (that is the absolute frequencies) will be represented by a Multiple Barchart with horizontal bars in the background.
  The width of the barcharts for the target variable can be constricted to the length of the horizontal bars.
 Additionally the \code{rmb} function allows to draw spineplots instead of the barcharts within each combination. On that score it can be seen as a generalization of Spineplots.
}
\usage{rmb(f, dset, hsplit = NULL, spine = FALSE,  hlcat = 1,
            eqwidth = FALSE, tfreq = "Id", max.scale = 1,
            use.na = FALSE, expected = NULL, mod.type = "poisson",
            resid.type = "pearson", use.expected.values = FALSE,
            resid.max = NULL, cut.rv = TRUE, cut.rs = 5,base = 0.2,
            mult = 1.5, colv = NULL, lab = TRUE, yaxis = TRUE,
            label = NULL, min.alpha = 0.1, base.alpha = 0.75,
            boxes = TRUE, lab.tv = FALSE, varnames = TRUE,
            abbr = FALSE, lab.cex = 1.2,\dots)}
\arguments{
  \item{f}{The \code{\link[stats]{formula}} specifying the variables in their given order with the last variable being the target variable. The left hand is either empty or denotes the frequency variable.}
  \item{dset}{The dataset in form of a frequency table (see \code{\link[stats]{ftable}} or \code{\link[extracat]{subtable}} for more information) with a column named \code{"Freq"} or in raw format (that is the rows represent the cases and the columns represent the variables). }
  \item{hsplit}{Logical vector with split directions where \code{TRUE} stands for horizontal splitting. The last (target) variable is always arranged on the x-axis. }
  \item{spine}{If \code{TRUE} a spineplot will be drawn instead of each barchart. This is recommended for binary target variables.}
  \item{hlcat}{A vector specifying the categories of the target variable which shall be highlighted (spineplot version) or plotted as bars(barchart version) in the given order. Specifying less than 2 indices will keep the original order in the barchart version.}
  \item{eqwidth}{If \code{TRUE} the bar length of the multiple barchart in the background no longer restricts the width of the barcharts for the relative frequencies of the target variable.}
  \item{tfreq}{The absolute frequencies used for the multiple barchart in the background can be transformed by means of \code{"log"} or square root(\code{"sqrt"}). No transformation is represented by \code{"Id"}.}
  \item{max.scale}{The maximum value of the probability (y-axis) scale for each combination. Unsurprisingly the default is 1. The axis will be drawn if \code{yaxis} is \code{TRUE}.}
  \item{use.na}{If \code{TRUE} missing values will be changed to a level \code{"N/A"} and else (which is default) the function \link[stats]{na.omit} will be called to reduce the dataset to complete cases only.}
  \item{expected}{A list of integer vectors denoting the interaction terms in the poisson or proportional odds model. If undefinded no residual shadings will be used.}
  \item{mod.type}{Either \code{"poisson"} or \code{"polr"}. See \link[stats]{glm} and \link[MASS]{polr}.}
  \item{resid.type}{ \code{"pearson"}, \code{"deviance"}, \code{working}, \code{partial}  or  \code{"response"}. For \code{polr} models only the latter is available. }
  \item{use.expected.values}{A logical specifying whether or not to use the frequencies predicted by the model instead of the observed frequencies.}
  \item{resid.max}{The maximum scale value for the residuals. If undefined it will be chosen automatically.}
  \item{cut.rv}{A logical. If \code{FALSE} the residual shading alpha values will be chosen exactly and otherwise depending on the specified scale intervals.  }
  \item{cut.rs}{The number of cuttingpoints for the residual scale.}
  \item{base}{The maximum proportion of the total plot width which is used for the gaps.}
  \item{mult}{The incremental multiplier for the gaps of different dimensions. }
  \item{colv}{A vector defining the colors of the bars or \code{NULL} for \link[grDevices]{rainbow} colors. Has no effect if an \code{expected} model is defined.}
  \item{lab}{If \code{TRUE} the plot size will decrease in order to make place for a labeling. }
  \item{yaxis}{If \code{TRUE} a vertical axis will be drawn at both sides of the plot. This is recommended when changing the \code{max.scale} argument.}
  \item{label}{An optional vector of logicals defining which variables shall be labelled. Has no effect if \code{lab} if \code{FALSE}.}
  \item{min.alpha}{In case of \code{eqwidth = T} alpha blending with respect to the corresponding weight is applied to the background color of the bars. In order to save very sparse combinations from disappearing there is a minimum alpha value. }
  \item{base.alpha}{A basic alpha value which will be applied to the bar colors. Does also work with residual shadings but not with colors chosen individually by setting \code{colv}}
  \item{boxes}{Should the labels be surrounded by boxes?}
  \item{lab.tv}{Should the target variable be included in the labeling?}
  \item{varnames}{Should the variable names be shown as labels?}
   \item{abbr}{If \code{TRUE} the labels will automatically be abbreviated (3 characters) using the \link[base]{abbreviate} function.}
   \item{lab.cex}{The fontsize multplier.}
   \item{\dots}{further arguments. Usually not necessary.}
}
\details{
  See \cite{upcoming paper (2009)}.
}
\value{
  No return value. 
}
\references{
  Alexander Pilhoefer \emph{ New approaches in visualization of categorical data:
       R-package extracat}\cr
  Journal of Statistical Software, submitted March 2010
}
\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany
}
\note{

}
\seealso{
  
}
\examples{
    data(housing)
    # simple example
    rmb(f = ~Type+Infl+Cont+Sat, dset = housing, mult = 2,
        hsplit = c(FALSE,TRUE,TRUE,FALSE), abbr = TRUE)
    
    # with sqrt-transformation and horizontal splits only
    rmb(f = ~Type+Infl+Cont+Sat, dset = housing, mult = 2,
        hsplit = c(TRUE,TRUE,TRUE,TRUE), tfreq = "sqrt",  abbr = TRUE)
    
    # a generalized spineplot with the first category highlighted
    rmb(f = ~Type+Infl+Cont+Sat, dset = housing, spine = TRUE, 
        hlcat = 1, mult = 2, hsplit = c(TRUE,FALSE,TRUE,TRUE), 
        tfreq = "sqrt", abbr = TRUE)
    
    # a generalized spineplot with all categories highlighted 
    # in a changed order
    rmb(f = ~Type+Infl+Cont+Sat, dset = housing, spine = TRUE,
        hlcat = c(3,1,2), mult = 2, hsplit = c(TRUE,FALSE,TRUE,TRUE),
        tfreq = "sqrt", abbr = TRUE)
    
    # with equal widths
    rmb(f = ~Type+Infl+Cont+Sat, dset = housing, eqwidth = TRUE,
        mult = 2, hsplit = c(TRUE,FALSE,TRUE,TRUE), lab.tv = TRUE,
        abbr = TRUE)
        
    # with equal widths, residual shadings and expected values
    rmb(f = ~Type+Infl+Cont+Sat, dset = housing, eqwidth = TRUE,
        mult = 2, hsplit = c(TRUE,FALSE,TRUE,TRUE), lab.tv = TRUE,
        abbr = TRUE, expected = list(c(1,2),c(2,3),c(3,4)),
        use.expected.values = TRUE)

}

