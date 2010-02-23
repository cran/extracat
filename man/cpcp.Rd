\name{cpcp}
\Rdversion{1.1}
\alias{cpcp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parallel Coordinates Plot for categorical data
}
\description{
This function provides a parallel coordinates plot for categorical as well as continuous data based on the \link[iplots]{ipcp} function in the \pkg{iplots} package.
It applies sorted numeric point sequences to the categories which indicate the relative frequencies and allow a sensible interactive highlighting.
There are options to change the rule for the gaps between these sequences and to apply an additional ordering algorithm.
}
\usage{cpcp(V, ord = NULL, freqvar = NULL, numerics = NULL,
        gap.type = "equal.tot", na.rule = "omit", spread = 0.3,
        gap.space = 0.2, psp = FALSE, scr.res = c(1280, 1024),
        sort.individual = FALSE, jitter = FALSE)}
\arguments{
  \item{V}{
data frame containing the variables for the plot. A special choice and order can be made by setting the second argument \code{ord}.If \code{V} is a frequency table (see e.g. \link[stats]{ftable}) the frequency variable has to be named "Freq".
}
\item{freqvar}{
    Optional name of a frequency variable in \code{V}. 
    If \code{V} contains a variable called \code{"Freq"} (see \link[stats]{ftable}) 
    it will be defined as frequency variable if \code{freqvar} is unspecified.
}
\item{numerics}{
    An integer vector specifying the indices of numeric variables.
    Variables of \link[base]{class} \code{"numeric"} are indicated automatically but integer variables are not.
}
  \item{ord}{
An integer vector specifying the indices and order of the variables.
}
  \item{gap.type}{
The rule for the gaps between different categories: \code{"equal.gaps"} will lead to equal gap sizes in the whole plot (and thus depends on the maximum number of categories). 
\code{"equal.tot"} will force the gaps within each variable to sum up to \code{gap.space}.
\code{"spread"} will lead to equal distances between the centerpoints of each category.
In the latter case the \code{spread} value denotes the maximal spread width around each centerpoint.
}
%  \item{num.last}{
%If \code{TRUE} all numerical variables will be arranged at the right end side of the plot (after the categorical ones).
%}
  \item{na.rule}{
\code{"omit"} will omit all incomplete cases from \code{V}. Otherwise all \code{NA}s will be transformed to categories called \code{"n/a"}. 
}
  \item{spread}{
 The maximum width of a single category if \code{gap.type} is set to \code{"spread"}. Should lie between 0 and 1.
}
  \item{gap.space}{
The space parameter for each gap. Should lie in \eqn{[0,1)}.
}
\item{psp}{
If \code{TRUE} parallel spineplots will be drawn below the plot.
}
\item{scr.res}{
The screen resolution used for the plot size.
}
\item{sort.individual}{
If \code{TRUE} all polygones with the same categories in two neighboring variables will be drawn together in the right variable.
}
\item{jitter}{logical: whether or not to jitter integer variables which are treated as a numeric variable.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
This function is based on the \pkg{iplots} package which will be installed automatically on the first run. For unix systems it is also necessary to run it from the \pkg{JGR} console.
To receive \pkg{JGR} please visit \url{http://www.rosuda.org/software} and download the JGR launcher which will start the console automatically and install all other required packages on its first run.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link[iplots]{ipcp}
\link[extracat]{listen}
\link[extracat]{resort}
}
\examples{
\dontrun{
cpcp(V = housing)
}
}
\keyword{ categorical }
\keyword{ parallel coordinates }
\keyword{ interactive }
