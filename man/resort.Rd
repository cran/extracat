\name{resort}
\alias{resort}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rearrangement of cpcp plots
}
\description{
Rearranging of the cpcp plot after a new selection.
}
\usage{resort(V = iset(which = iset.cur()),
        sel = iset.selected(), listen = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
An iset created by the \code{cpcp} function for which the sorting is to be computed.
Usually it is not necessary to choose \code{V} individually.
}
  \item{sel}{
The indices of the selected cases.
}
  \item{listen}{
If \code{TRUE} the function \code{listen} will be (re-)called after the updating process is done.
}
}
\details{
If called individually after making a new selection the function will simply rearrange the selected cases at the top of each category.
Can also be called indirectly by \link{listen}. 
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


\seealso{
\link[extracat]{listen}
\link[extracat]{cpcp}
}
\examples{
\dontrun{
    cpcp(V = housing, ord = c(4,3,1,2))
    # selection via console
    s = iset()
    iset.select(what = which(ivar.data(s$Type) == "Terrace") )
    resort()
}
}

