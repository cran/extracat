\name{fluctile}
\alias{fluctile}

\title{
fluctuation diagrams
}
\description{
Create a fluctuation diagram from a multidimensional table.
}
\usage{
fluctile(tab, gap.prop = 0.1, hsplit = TRUE, env = NULL, abbrev = 3, border = NULL, label = TRUE,...  )
}

\arguments{
  \item{tab}{
The table which is to be plotted.
}
  \item{gap.prop}{
proportion of the gaps between the rows/columns within each block.
}
  \item{hsplit}{
A logical for alternating columns and rows or a vector of logicals with TRUE for each variable on the x-axis.
}
  \item{env}{
Not necessary.
}
  \item{abbrev}{
abbreviate the labels to a given number of characters.
}
  \item{border}{
The proportion of the borders. The labels will be plotted in this region.
}
  \item{label}{
Whether or not to plot labels.
}
  \item{\dots}{
dots
}
}

\value{
The viewport tree behind the graphic.
}

\author{
Alexander Pilhoefer\cr
Department for Computer Oriented Statistics and Data Analysis\cr
University of Augsburg\cr
Germany
}
\note{
This was part of the Google Summer of Code 2011.
}


\seealso{
\link[vcd]{mosaic}
}

\examples{
M <- bloma(4,4,1.3)
fluctile(M)

M2 <- optile(M)
fluctile(M2)

data(CPScluster)
cpsX = subtable(CPScluster,c(26,34,38,39),allfactor=TRUE)

ssC <- optile(cpsX,presort=FALSE, fun = "csvd", return.data=TRUE, method=NULL)
fluctile(xtabs(Freq~.,data=ssC[,-4]))



}

