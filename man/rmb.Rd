\name{rmb}
\Rdversion{1.1}
\alias{rmb}
\alias{rmb.formula}
\alias{rmb.table}
\alias{rmb.ftable}
\concept{Multiple Barchart for relative frequencies and generalized Spineplots}
\title{
  Multiple Barchart for relative frequencies and generalized Spineplots
}
\description{
  The \code{rmb} function basically produces a Multiple Barchart for the relative frequencies of some target categories within each combination of the explanatory variables. The weights of those combinations (that is the absolute frequencies) are represented in the total with of the corresponding barchart.
  The result is a graphic which allows to read the conditional target distributions exactly from the graphic without losing the information about the importance (in the sense of the number of observations) of the different combinations.
 
 Additionally the \code{rmb} function allows to draw spineplots instead of the barcharts within each explanatory combination. On that score it can be seen as a generalization of Spineplots.
}
\usage{
\method{rmb}{formula}(formula, data, col.vars = NULL, spine = FALSE,  cat.ord = NULL,  eqwidth = FALSE,
                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL,
                mod.type = "poisson", resid.type = "pearson",
                use.expected.values = FALSE, resid.max = NULL, cut.rs = 5,                
                gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), yaxis = TRUE, label = TRUE,
                min.alpha = 0.1, boxes = TRUE, lab.tv = FALSE, varnames = TRUE, abbrev = FALSE, lab.cex = 1.2,\dots)
\method{rmb}{ftable}(x, col.vars = NULL, spine = FALSE,  cat.ord = NULL,  eqwidth = FALSE,
                freq.trans = NULL, max.scale = 1,  use.na = FALSE, expected = NULL,
                mod.type = "poisson", resid.type = "pearson",
                use.expected.values = FALSE, resid.max = NULL, cut.rs = 5,                
                gap.prop = 0.2, gap.mult = 1.5, col = "hcl",col.opt = list(), yaxis = TRUE, label = TRUE,
                min.alpha = 0.1, boxes = TRUE, lab.tv = FALSE, varnames = TRUE, abbrev = FALSE, lab.cex = 1.2,\dots)
}
\arguments{
 \item{x}{A table either of class \code{table} or of class \code{ftable}. The latter also implicitly defines the the order in which the variables will be added to the plot. The arguments formula and data will be omitted.}
  \item{formula}{The \code{\link[stats]{formula}} specifying the variables in their given order with the last variable being the target variable. The left hand is either empty or denotes the frequency variable.}
  \item{data}{The dataset in form of a frequency table (see \code{\link[stats]{ftable}} or \code{\link[extracat]{subtable}} for more information) with a column named \code{"Freq"} or in raw format (that is the rows represent the cases and the columns represent the variables). }
  \item{col.vars}{Logical vector with split directions where \code{TRUE} stands for horizontal splitting. The last (target) variable is always arranged on the x-axis. }
  \item{spine}{If \code{TRUE} a spineplot will be drawn instead of each barchart. This is recommended for binary target variables.}
  \item{cat.ord}{A vector specifying the categories of the target variable which will be visualized in the specified order. The default is to use all categories.}
  \item{eqwidth}{If \code{TRUE} the bar length of the multiple barchart in the background no longer restricts the width of the barcharts/spineplots for the relative frequencies of the target variable.}
  \item{freq.trans}{This parameter allows to transform the absolute frequencies used for the underlying multiple barchart. Possible values are \code{"log"}, \code{"sqrt"} or \code{list("krt",k)}. The latter stands for the k-th root transformation.}
  \item{max.scale}{The maximum value of the probability (y-axis) scale for each combination. Unsurprisingly the default is 1. The axis will be drawn if \code{yaxis} is \code{TRUE}.}
  \item{use.na}{If \code{TRUE} missing values will be changed to a level \code{"N/A"} and else (which is the default) the function \link[stats]{na.omit} will be called to reduce the dataset to complete cases only.}
  \item{expected}{A list of integer vectors denoting the interaction terms in the poisson or proportional odds model. If undefinded no residual shadings will be used. Can also be set to \code{TRUE} which stands for the default independence model: the null model in the sense of logistic regression  models. In the Poisson case the log-linear equivalent is to use a saturated model for the explanatory variables with no interaction term including the target variable.}
  \item{mod.type}{Either \code{"poisson"} or \code{"polr"}. See \link[stats]{glm} and \link[MASS]{polr}.}
  \item{resid.type}{ \code{"pearson"}, \code{"deviance"}, \code{working}, \code{partial}  or  \code{"response"}. For \code{polr} models only the latter is available. }
  \item{use.expected.values}{A logical specifying whether or not to use the frequencies predicted by the model instead of the observed frequencies.}
  \item{resid.max}{The maximum scale value for the residuals. If undefined it will be chosen automatically.}
   \item{cut.rs}{The number of cuttingpoints for the residual scale.}
  \item{gap.prop}{The maximum proportion of the total plot width which is used for the gaps.}
  \item{gap.mult}{The incremental multiplier for the gaps of different dimensions. The gaps corresponding to any one variable are \code{gap.mult} times larger than those corresponding to the next variable on the same axis.}
  \item{col}{Either a vector defining the colors of the bars or a name specifying a palette: \code{"hsv"} and \code{"rgb"} for hsv-based rainbow colors, \code{"hcl"} for hcl-based rainbow colors (default), \code{"div"} or \code{"diverge"} for hcl-based diverging colors and finally \code{"seq"} or \code{"sequential"}
  for hcl-based sequential colors.Additional arguments can be specified via the \code{col.opt} argument according to the underlying functions in the colorspace package, e.g. \link[colorspace]{rainbow_hcl}. For the hsv-based colors see \link[grDevices]{rainbow}. Specifying a color or palette has no effect if an \code{expected} model is defined.}
  \item{col.opt}{Further options for the color palettes. See e.g. \link[colorspace]{rainbow_hcl} or \link[grDevices]{rainbow}.}
   \item{yaxis}{If \code{TRUE} a vertical axis will be drawn at both sides of the plot. This is recommended when changing the \code{max.scale} argument.}
  \item{label}{Either a logical specifying whether or not to draw labels or a numeric vector defining which variables shall be labelled. }
  \item{min.alpha}{In case of \code{eqwidth = T} alpha blending with respect to the corresponding weight is applied to the background color of the bars. In order to save very sparse combinations from disappearing there is a minimum alpha value. }
  \item{boxes}{Should the labels be surrounded by boxes?}
  \item{lab.tv}{Should the target variable be included in the labeling?}
  \item{varnames}{Should the variable names be shown as labels?}
   \item{abbrev}{An integer value specifying the number of characters to which the labels will automatically be abbreviated. The logical \code{TRUE} leads to 3-character-abbreviations.}
   \item{lab.cex}{The fontsize multplier.}
   \item{\dots}{further arguments. Usually not necessary.}
}
\details{
  A similar way to regard the graphic is the following: A Multiple Barchart of the explanatory variables is drawn with bars in horizontal direction. Then within each of the resulting bars a barchart of the conditional distribution of the target variable is drawn with bars in vertical direction.
}
\value{
 \code{invisible(TRUE)}
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
  \code{\link[vcd]{mosaic}} for an implementation of mosaicplots.
}
\examples{
    data(housing)
    # simple example
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, gap.mult = 2,
        col.vars = c(FALSE,TRUE,TRUE,FALSE), abbrev = 3)
    
    # with sqrt-transformation and horizontal splits only
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, gap.mult = 2,
        col.vars = c(TRUE,TRUE,TRUE,TRUE), freq.trans = "sqrt",  abbrev = TRUE)
    
    # a generalized spineplot with the first category highlighted
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, spine = TRUE, 
        cat.ord = 1, mult = 2, col.vars = c(1,3,4), 
        freq.trans = list("krt",3), abbrev = 2)
    
    # a generalized spineplot with all categories highlighted 
    # in a changed order
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, spine = TRUE,
        cat.ord = c(3,1,2), gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE),
        freq.trans = "sqrt", abbrev = TRUE)
    
    # with equal widths
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, eqwidth = TRUE,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), lab.tv = TRUE,
        abbrev = 2)
        
    # with equal widths, residual shadings and expected values
    rmb(formula = ~Type+Infl+Cont+Sat, data = housing, eqwidth = TRUE,
        gap.mult = 2, col.vars = c(TRUE,FALSE,TRUE,TRUE), lab.tv = TRUE,
        abbrev = 2, expected = list(c(1,2),c(2,3),c(3,4)),
        use.expected.values = TRUE)
        
     # using an ftable to create the plot
     tt = xtabs(Freq~Type+Cont+Infl+Sat, data = housing)
     ft = ftable(tt, col.vars= c(1,4))
      rmb(tt)
      rmb(ft)

}

