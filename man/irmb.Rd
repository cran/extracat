\name{irmb}
\alias{irmb}
\alias{irmb.formula}
\alias{irmb.table}
\alias{irmb.ftable}

\concept{Interactive Multiple Barchart for relative frequencies and generalized Spineplots}
\title{
  Interactive Multiple Barchart for relative frequencies and generalized Spineplots
}
\description{
  The \code{irmb} function is an interactive version of the basic \link{rmb} function. It provides sliders for x- and y-scales, variable inclusion/exclusion, axis changes, switching between barcharts, spineplots and mosaic plots as well as model selection for residual shadings and a few more options.
  It depends on the \pkg{iWidgets} package and therefore only works when called from the \pkg{JGR} console.}
\usage{

\method{irmb}{formula}(formula, data, gap.prop = 0.2,
               gap.mult = 1.5, col = "hcl", col.opt =list(), Z=100, abbrev = FALSE,
               use.na = FALSE, expected = TRUE,  resid.type = "pearson", 
               max.resid.scale = NULL, cut.rs = 7, yaxis = TRUE, min.alpha = 0.1,
               boxes = TRUE, lab.tv = FALSE, varnames = TRUE, lab.cex = 1.2,\dots)
\method{irmb}{ftable}(x , gap.prop = 0.2, gap.mult = 1.5, col = "hcl", col.opt=list(),
			Z=100, abbrev = FALSE, use.na = FALSE, expected = TRUE,
			resid.type = "pearson", max.resid.scale = NULL, cut.rs = 7,
			yaxis = TRUE, min.alpha = 0.1, boxes = TRUE, lab.tv = FALSE,
			varnames = TRUE, lab.cex = 1.2, \dots)
              }

\arguments{
\item{x}{A table either of class \code{table} or of class \code{ftable}. The arguments formula and data will be omitted.}
  \item{formula}{Either a \code{\link[stats]{formula}} a \code{\link[stats]{formula}} or an \code{\link[stats]{formula}}. If \code{x} is a formula it also specifies the order of the variables in the plot. The las variable will be the target variable. The left hand side of the formula is either empty or denotes a frequency variable.}
  \item{data}{The dataset in form of a frequency table (see \code{\link[stats]{ftable}} or \code{\link[extracat]{subtable}} or a \code{data.frame}. Ignored if \code{x} is a table. }
   \item{gap.prop}{The maximum proportion of the total plot width which is used for the gaps.}
  \item{gap.mult}{The incremental multiplier for the gaps of different dimensions. The gaps corresponding to any one variable are \code{gap.mult} times larger than those corresponding to the next variable on the same axis.}
  \item{col}{Either a vector defining the colors of the bars or a name specifying a palette: \code{"hsv"} and \code{"rgb"} for hsv-based rainbow colors, \code{"hcl"} for hcl-based rainbow colors (default), \code{"div"} or \code{"diverge"} for hcl-based diverging colors and finally \code{"seq"} or \code{"sequential"}
  for hcl-based sequential colors.Additional arguments can be specified via the \code{col.opt} argument according to the underlying functions in the colorspace package, e.g. \link[colorspace]{rainbow_hcl}. For the hsv-based colors see \link[grDevices]{rainbow}. Specifying a color or palette has no effect if an \code{expected} model is defined.}
  \item{col.opt}{Further options for the color palettes. See e.g. \link[colorspace]{rainbow_hcl} or \link[grDevices]{rainbow}.}
    \item{Z}{The number of intervals for the sliders.}
  \item{use.na}{Logical, specifying whether or not to drop cases with missing entries. If not, a category \code{"N/A"} will be created in each variable where missings occur.}
  \item{expected}{If \code{TRUE} residual shadings according to loglinear models will be available.}
  \item{resid.type}{ \code{"pearson"}, \code{"deviance"}, \code{"working"}, \code{"partial"}  or  \code{"response"}. The underlying model is a \link[stats]{glm} poisson model. The mosaicplot option will always use Pearson residuals. }
  \item{max.resid.scale}{If undefined the scale of the residuals will update with the model.}
  \item{cut.rs}{The number of cuttingpoints for the residual scale.}
  \item{yaxis}{If \code{TRUE} a vertical axis will be drawn at both sides of the plot. This is recommended when using the y-axis slider.}
  \item{min.alpha}{In case of \code{eqwidth = T} alpha blending with respect to the corresponding weight is applied to the background color of the bars. In order to save very sparse combinations from disappearing there is a minimum alpha value. }
    \item{boxes}{Should the labels be surrounded by boxes?}
  \item{lab.tv}{Should the target variable be included in the labeling?}
  \item{varnames}{Should the variable names be shown as labels?}
  \item{abbrev}{If \code{TRUE} the labels will automatically be abbreviated (3 characters) using the \link[base]{abbreviate} function.}
  \item{lab.cex}{The fontsize multplier.}
  \item{\dots}{Further arguments. Usually not necessary.}
}
\details{

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
  This function is based on \pkg{iWidgets} and therefore needs the \pkg{JGR} console to run. It is not sufficient to load the \pkg{JGR} package.
  For \pkg{iWidgets} please visit \url{http://www.rforge.net} to download the \pkg{iWidgets} source file or install it from the console using:\cr
  \code{install.packages("iWidgets",,'http://rforge.net/',type='source')}\cr
  For \pkg{JGR} please visit \url{http://www.rosuda.org/software} and download the JGR launcher which will start the console automatically and install all other required packages on first run.
  
}
\seealso{
  
}
\examples{
\dontrun{
    data(housing)
    # example 1
    irmb(~Type+Infl+Cont+Sat,data=housing,abbrev=4)
    
    # example 2
    irmb(~Type+Infl+Cont+Sat,data=housing,expected = FALSE,
    	gap.mult=2,abbrev=4)
}
 }

