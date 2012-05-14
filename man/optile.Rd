\name{optile}
\alias{optile}
\alias{optile.list}

\title{
optimization of category orders for nominal categorical datasets and tables
}
\description{
This function will take a categorical data object (\code{data.frame, table, ftable, matrix, array}) and optimize its category orders graphical representations.
It offers an interface which will by default return the same type of object that has been passed to the function such that it is possible to
write  \code{myplot( optile(x) )} for an optimized version of \code{myplot(x)}. 
The interface is extendible using ones own functions should they fulfill the requirements (see details).
}
\usage{
optile(x, fun = "class", presort = FALSE, foreign = NULL,
 args = list(), perm.cat = TRUE, method = NULL, iter = 1,
 past = NULL, scale = FALSE, freqvar = NULL,
 return.data = TRUE, return.type = "data.frame", vs = 0, tree = NULL,\dots)
 
 \method{optile}{list}(x, fun = "class", presort = FALSE, foreign = NULL,
 args = list(), perm.cat = TRUE, method = NULL, iter = 1,
 past = NULL, scale = FALSE, freqvar = NULL,
 return.data = TRUE, return.type = "table", vs = 0, tree = NULL,
  k = NULL, h = NULL, \dots)
}

\arguments{
  \item{x}{
The categorical data of one of the following classes:
\code{data.frame, table, ftable, matrix, array} 
}
  \item{fun}{
The optimization function. Currently available are:
\code{class, hamming, casort, csvd , rmca} and \code{preclass}. \code{Class, hamming} and \code{preclass} require \code{foreign = ".Call"}
See Details.
}
  \item{presort}{
Whether or not to use \code{preclass} for an initial pre-ordering
}
  \item{foreign}{
Where to find the optimization function \code{fun}. \code{NULL} if it is an R function.
}
  \item{args}{
further arguments which will be passed to \code{fun}.
}
  \item{perm.cat}{
All variables with a fixed category order (e.g. ordinal variables) should have a \code{FALSE} here. Has no effect if \code{fun = "casort"}
}
  \item{method}{
Either \code{NULL, "joint"} or \code{"stepwise"}. \code{NULL} will pass the dataframe to \code{fun, "joint"} will pass the Burt matrix (with all pairwise tables) and \code{"stepwise"} will repeatedly call \code{fun} for 2, 3, 4, and so on variables. See Details.
}
  \item{iter}{
Some optimizations depend on the initial category orders (e.g. \code{"class"} and \code{"preclass"}). This will compute the optimization for \code{iter} random starting points and choose the best.
}
  \item{past}{
Not implemented yet.
}
  \item{scale}{
Not implemented yet. ( Function to scale resulting criteria to make them comparable)
}
  \item{freqvar}{
The name of the frequency variable, if any.
}
  \item{return.data}{
Whether to return the data or just the new orders.
}
  \item{return.type}{
The class of the object which will be returned. Defaults to the input type.
}
  \item{vs}{
An optional version number. \code{"hamming"} is currently equivalent to \code{"class"} and \code{vs =1}}
  \item{tree}{
A list whose entries are either a tree object like they are obtained from \link[stats]{hclust} or the string \code{"hc"}. If the ith entry is a tree object, the ith variable should to be the result from cutting the tree into \code{dim(x)[i]} clusters via \link[extracat]{subtree}. \code{"hc"} will compute a hierarchical clustering for the rows and columns with arguments specified in \code{args}.
}
  \item{k}{
A vector of integers specifying the numbers of clusters into which the tree objects shall be cut. See \link[stats]{hclust}. }
  \item{h}{
Instead of a number of clusters \code{k} the height at which the dendrogram shall be cut can be specified. See \link[stats]{hclust}.}
  \item{\dots}{
dots
}
}
\details{
How to add a custom optimization function:

It is possible to use custom functions for the optile interface as long as they meet the following requirements:

The function should have the form\cr
\code{fun( data, dims, perm.cat, ... )} or \cr
\code{foreign( "fun", data, dims, perm.cat, ...)} \cr
where fun is the name of the function and foreign is \code{".Call", ".C", ...}

The function returns a vector of the new category orders and the resulting criterion, e.g. \code{c( 1,3,5,2,4, 5,4,3,2,1,6,7, 78978 )}
(CURRENTLY: 0..k-1 instead of 1..k)

The function works with integer input weights/frequencies.

\code{dims} is a vector with the number of categories for each variable and \code{perm.cat} is a 0/1 vector which indicates whether or not to change the category order of a variable.

There are three possible types for the \code{data} argument of \code{fun} which can be set via \code{method}:

The argument \code{method} can be one of \code{NULL, "stepwise"} or \code{"joint"}.
The default \code{method = NULL} indicates that \code{fun} accepts a multidimensional table as for instance can be produced via \code{xtabs}.

If \code{method = "joint"} a Burt matrix is computed and passed to \code{fun}. This matrix contains all 2-way tables and is an easy way to pass the joint 2-way-interaction data to a function. For instance \code{"fun=casort"} uses this data representation.

\code{method = "stepwise"} or \code{method = "sw"} passes \code{fun, data, foreign} as well as any \code{args} to a function called \code{steptile} which initially builds a 2-way table of the first pair of variables, passes it to fun and stores the computed category orders. Afterwards the other variables are added one by one. i.e. in a step for the k-th variable the function passes a k-way table to \code{fun} and a new category order for this variable is computed given the (already fixed) category orders of the variables 1 to k-1. E.g. for the \code{"class"} and \code{"hamming"} criteria this is a lot faster than the multivariate (method = NULL) version.


CURRENTLY AVAILABLE SORTING FUNCTIONS:

\code{"class"} and \code{"hamming"}: both maximize the number of pairs of observations which are  consistent with a pseudo diagonal,
i.e. for such a pair all factor levels of one observation should be higher than the corresponding ones of the other observation. The \code{"hamming"} version additionally uses weights according to the hamming distance of the observations.

\code{"casort"}: computes a correspondence analysis (SVD) and sorts by the first coordinate vector of each dimension. For more than two dimensions Multiple CA based on the Burt matrix is used.
It is also useful to compute initial category orders.

\code{"rmca"}: 
Adopts the idea of CA for k > 2 dimensions without dropping information: For each dimension d = 1..k with categories d1...dr compute the scaled average k-1 dimensional profile sdd and perform an SVD of (sd1...sdr)-sdd. Like in correspondence analysis the first coordinate vector is used for the reordering.

\code{"csvd"}:
For each variable d in 1..k (iteratively)
compute the cumulative sums over the multidimensional table for each variable except d. Transform this multidimensional table to an r x s matrix with r being the number of categories of variable d and s being the product of these numbers for all other variables.
Resort the categories of variable d by the first coordinate vector of an SVD of that matrix.
Repeat this procedure for all variables in turn until a stopping criterion is met.
Idea: for any variable h != d we have h1 < h2 < ... < hx due to the cumulative sums. Hence the current order of the categories will (tend to) be the same as in the coordinates of the svd which means that the svd computes coordinates for variable d with respect to the current category orders of the other variables. The algorithm uses \code{casort} for an initial solution to start from.


\code{"distcor"}:
Two-way tables or matrices can also be optimized by means of the distance correlation. See \link[extracat]{dcor}.

\code{"preclass"}:
This function is used for the presorting and is written in C (\code{foreign =".Call"}).
It iteratively sorts the categories of one variable at a time. Therefore it computes the average over the remaining dimensions and scales the profiles of each category as well as the average profile. It then computes the classification criterion between each category profile and the average profile which results in one value per category. The categories are then sorted by this criterion.
The procedure is very quick and yields good results in most cases. It strongly depends on the initial category orders as do the classification or hamming algorithms. It is possible to choose this for the main sorting function by setting \code{fun = "preclass"}.

}
\value{
The function returns the data. The return type is by default the same as the input type but can be redefined via \code{return.type}.
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


\examples{
\dontrun{
############ ------------ EXAMPLE I ------------ ############
# ----- Cluster results from the Current Population Survey ----- #
data(CPScluster)
cpsX = subtable(CPScluster,c(5, 12, 26, 34, 38, 39), allfactor=TRUE)

# joint and stepwise classification criterion. both use a very quick algorithm.
ss <- optile(cpsX,presort=TRUE, return.data=TRUE, method="joint")
ss2 <- optile(cpsX,presort=TRUE, return.data=TRUE, method="sw")

# original cpcp plot
cpcp(cpsX)

# cpcp for joint algorithm
cpcp(ss)

# cpcp and fluctuation for the stepwise algorithm
# (should be good for pcp plots and hierarchical plots)
fluctile(xtabs(Freq~.,data=ss2[,-4]))
cpcp(ss2)

# The multivariate algorithm
ss3 <- optile(cpsX,presort=TRUE, return.data=TRUE, method=NULL)
cpcp(ss3)

# cpcp for casort algorithm
ssca <- optile(cpsX,presort=FALSE, fun = "casort", return.data=TRUE, method="joint")
cpcp(ssca)

# cpcp for rmca algorithm results. works better for the dmc data
ssR <- optile(cpsX,presort=FALSE, fun = "rmca", return.data=TRUE, method=NULL)
cpcp(ssR)


# cpcp for csvd algorithm
ssC <- optile(cpsX,presort=FALSE, fun = "csvd", return.data=TRUE, method=NULL)
fluctile(xtabs(Freq~.,data=ssC[,-4]))
cpcp(ssC)

# cpcp for presort algorithm with 20 iterations
ssP <- optile(cpsX,presort=FALSE, fun = "preclass",
return.data=TRUE, method=NULL, foreign = ".Call",iter=20)
cpcp(ssP)

############ ------------ EXAMPLE II ------------ ############
# ------------------------------- Italian wines ------------------------------ #
library(MMST)
data(wine)

swine <- scale(wine[,1:13])
kmd <- data.frame(wine$class, replicate(9, kmeans(swine, centers = 6)$cluster) )
kmd <- subtable(kmd, 1:10, allfactor = TRUE)

cpcp(kmd)

# there is a good joint order and hence the joint result is better than the stepwise
kmd2 <- optile(kmd, method = "sw")
kmd3 <- optile(kmd, method = "joint")

cpcp(kmd2)
cpcp(kmd3)



############ ------------ EXAMPLE III ------------ ############
# ---------------- The BicatYeast microarray dataset  ---------------- #

# ----- with different clusterings for the genes ----- #
library(biclust)
data(BicatYeast)

Dby <- dist(BicatYeast)

hc1 <- hclust(Dby, method = "ward")
hc2 <- hclust(Dby, method = "average")
hc3 <- hclust(Dby, method = "complete")

hcc1 <- cutree(hc1, k = 6)
hcc2 <- cutree(hc2, k = 6)
hcc3 <- cutree(hc3, k = 6)

km1 <- kmeans(BicatYeast, centers = 6, nstart = 100, iter.max = 30)$cluster

library(mclust)
mc1 <- Mclust(BicatYeast, G = 6)$class

clusterings <- data.frame(hcc1,hcc2,hcc3,km1,mc1)
clusterings <- subtable(clusterings, 1:5, allfactor = TRUE)

clusterings2 <- optile(clusterings, method = "joint")
clusterings3 <- optile(clusterings, fun = "casort")

cpcp(clusterings2)

# a fluctuation diagram of all but the avg. clustering
fluctile(xtabs(Freq~.,data=clusterings2[,-2]))

# compute agreement via Fleiss kappa in irr:
require(irr)
rawdata <- untableSet(clusterings2)
for(i in 1:5) levels(rawdata[,i]) <- 1:6
(kappam.fleiss(rawdata))
(kappam.fleiss(rawdata[,-2]))




############ ------------ EXAMPLE IV ------------ ############
# ------------------------- The Eisen Yeast data ------------------------- #
library(biclust)
data(EisenYeast)
SEY <- scale(EisenYeast)

Dby2 <- dist(SEY)

hc1 <- hclust(Dby2, method = "ward")
hc2 <- hclust(Dby2, method = "complete")

hcc1 <- cutree(hc1, k = 16)
km1 <- kmeans(scale(EisenYeast), centers = 16, nstart = 20, iter.max = 30)$cluster
optile( table(hcc1, km1) )


############ ------------ EXAMPLE V ------------ ############
# ------------------------- The Bicat Yeast data ------------------------- #

# how many clusters are a good choice for kmeans?
# one possible way to find out: 
# compute kmeans for 100 random initial settings, sort the results (clusters) 
# and compute their agreement
# e.g. via Fleiss' Kappa (available in package irr)

require(biclust)
data(BicatYeast)
require(irr)

st <- Sys.time()
fk <- NULL
for(k in 3:8){
	test <- subtable(replicate(100,kmeans(BicatYeast, centers = k)$cluster),1:100)
	test <- optile(test, fun = "casort")
	test <- optile(test, method="joint")
	test <- untableSet(test)
	for(i in 1:100) levels(test[,i]) <- 1:k
	fk <- c(fk,kappam.fleiss(test)$value)
}
Sys.time()-st
plot(x = 3:8, y = fk, type="l", lwd=2)

############ ------------ EXAMPLE VI ------------ ############
# ------------------------- hierarchical clustering ------------------------- #

# A list with hierarchical clustering objects:
require(amap)

hc1 <- hcluster(t(plants[,-1]), method="manhattan", link = "ward")
hc2 <- hcluster(t(plants[,-1]), method="manhattan", link = "complete")

hclist <- list(hc1, hc2)
tfluctile( optile(hclist, k= c(8,8) ) )

# or a table with corresponding tree objects:

tt <- table( subtree(hc1, 12)$data, subtree(hc2, 8)$data )

tfluctile(optile(tt, tree = list(hc1, hc2)))

# only one tree object, the other variable is free:

tt <- table( subtree(hc1, 8)$data, kk <- kmeans(t(plants[,-1]),centers=8)$cluster )

tfluctile(optile(tt, tree = list(hc1, NA)))

}




}

