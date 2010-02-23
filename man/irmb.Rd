\name{irmb}
\Rdversion{1.1}
\alias{irmb}
\concept{Interactive Multiple Barchart for relative frequencies and generalized Spineplots}
\title{
  Interactive Multiple Barchart for relative frequencies and generalized Spineplots
}
\description{
  The \code{irmb} function is an interactive version of the basic \link{rmb} function. It provides sliders for x- and y-scales, variable inclusion/exclusion, axis changes, switching between barcharts and spineplots, model selection for residual shadings as well as a highlighting-like choice of the categories in the spineplot.
  It depends on \pkg{iWidgets} which only works with \pkg{JGR}.}
\usage{irmb(f, dset, Z=100, use.na = FALSE,
 expected = TRUE, resid.type = "pearson",
 use.expected.values = FALSE, max.resid.scale = NULL,
 cut.rv = TRUE, cut.rs = 5, base = 0.2, mult = 1.5, 
 colv = NULL, yaxis = TRUE, min.alpha = 0.1, 
 base.alpha = 0.75, boxes = TRUE, lab.tv = FALSE,
 varnames = TRUE, abbr = FALSE, lab.cex = 1.2,\dots)}

\arguments{
  \item{f}{The \code{\link[stats]{formula}} specifying the variables in their given order with the last variable being the target variable. The left hand is either empty or denotes the frequency variable.}
  \item{dset}{The dataset in form of a frequency table (see \code{\link[stats]{ftable}} or \code{\link[extracat]{subtable}} for more information) with a column named \code{"Freq"} or in raw format (that is the rows represent the cases and the columns represent the variables). }
  \item{Z}{The number of intervals for the sliders.}
  \item{use.na}{Logical, specifying whether or not to drop cases with missing entries. If not, a category \code{"N/A"} will be created in each variable where missings occur.}
  \item{expected}{Should model residuals be used to define the colors of the rectangles or bars?}
  \item{resid.type}{ \code{"pearson"}, \code{"deviance"}, \code{"working"}, \code{"partial"}  or  \code{"response"}. The underlying model is a \link[stats]{glm} poisson model. }
  \item{use.expected.values}{Logical specifying whether or not to use the frequencies expected from the model instead of the observed ones in the dataset.}
  \item{max.resid.scale}{If undefined the scale of the residuals will update with the model.}
  \item{cut.rv}{A logical. If \code{FALSE} the residual shading alpha values will be chosen exactly and otherwise depending on the specified scale intervals.  }
  \item{cut.rs}{The number of cuttingpoints for the residual scale.}
  \item{base}{The maximal proportion of the total plot width which is used for the gaps.}
  \item{mult}{The incremental multiplier for the gaps of different dimensions. }
  \item{colv}{A vector defining the colors of the bars or \code{NULL} for \link[grDevices]{rainbow} colors. Has no effect if \code{expected = TRUE}.}
  \item{yaxis}{If \code{TRUE} a vertical axis will be drawn at both sides of the plot. This is recommended when using the y-axis slider.}
  \item{min.alpha}{In case of \code{eqwidth = T} alpha blending with respect to the corresponding weight is applied to the background color of the bars. In order to save very sparse combinations from disappearing there is a minimum alpha value. }
  \item{base.alpha}{A basic alpha value which will be applied to the bar colors. Does also work with residual shadings but not with colors chosen individually by setting \code{colv}}
  \item{boxes}{Should the labels be surrounded by boxes?}
  \item{lab.tv}{Should the target variable be included in the labeling?}
  \item{varnames}{Should the variable names be shown as labels?}
  \item{abbr}{If \code{TRUE} the labels will automatically be abbreviated (3 characters) using the \link[base]{abbreviate} function.}
  \item{lab.cex}{The fontsize multplier.}
  \item{\dots}{Further arguments. Usually not necessary.}
}
\details{
  The \code{hlcat} editable field allows to choose the categories for the spineplot mode. For example, \code{3,1,4} will put the third category on the bottom and then stack the first and fourth category. Category 2 will be left out.
  The model selection field has to contain the interaction terms to fit the \link{glm} poisson model. The single terms will be included automatically.
  The \code{hsplit} checkboxes can be used to switch variables from one axis to the other but it has no effect if the corresponding variable is the last one.
    
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
  This function is based on \pkg{iWidgets} and therefore needs the \pkg{JGR} console to run. 
  For \pkg{iWidgets} please visit \url{http://www.rforge.net} to download the \pkg{iWidgets} source file or install it from the console using:\cr
  \code{install.packages("iWidgets",,'http://rforge.net/',type='source')}\cr
  For \pkg{JGR} please visit \url{http://www.rosuda.org/software} and download the JGR launcher which will start the console automatically and install all other required packages on its first run.
  
}
\seealso{
  
}
\examples{
\dontrun{
    data(housing)
    # example 1
    irmb(~Type+Infl+Cont+Sat,dset=housing,expected = TRUE,
        eqwidth=FALSE,base=0.2,mult=2,lab.tv=TRUE,abbr=TRUE)
    
    # example 2
    irmb(~Type+Infl+Cont+Sat,dset=housing,expected = FALSE,
        eqwidth = TRUE,base=0.2,mult=2,lab.tv=FALSE,abbr=TRUE)
}
 }

