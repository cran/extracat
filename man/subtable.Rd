\name{subtable}
\alias{subtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
data.frame reduction
}
\description{
Reduces a dataframe into a frequency table with prespecified entries. Uses a modified version of the \link[plyr]{count} function which also accepts weights.
Zero-entries can be included or excluded and the variables can be coerced into factors if necessary.
}
\usage{subtable(data, cols, freqvar = NULL, 
		keep.zero = FALSE, allfactor = FALSE, return.type=class(data))}

\arguments{
  \item{data}{
The \code{data.frame} to reduce.
}
  \item{cols}{
An ordered integer vector containing the indices of the colums to keep.
}
  \item{keep.zero}{
A logical indicating whether to include zero-cases in the output.
}
  \item{allfactor}{
A logical indicating whether to convert all variables into factor variables. Integer variables will be applied a fitting (non-lexicographic) level order.
}
\item{freqvar}{
    Optional name of a frequency variable in \code{V}. 
    If \code{dset} contains a variable called \code{"Freq"} (see \link[stats]{ftable}) 
    it will be defined as frequency variable if \code{freqvar} is unspecified.
}
\item{return.type}{
\code{"table"} will return a table and anything else will result in a data.frame.
}
}
\details{
This function uses a modified version of \link[plyr]{count} in order to provide the arguments \code{keep.zero} and \code{allfactor} and to be able to handle both raw datasets and datasets with a frequency variable.
}
\value{
A \code{data.frame} including a \code{"Freq"} variable.
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


\examples{
hs2 = subtable(housing,c(3,1))
summary(hs2)
}
