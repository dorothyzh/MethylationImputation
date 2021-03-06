\name{MetIm}
\alias{MetIm}
\title{MethyImpute packages}
\description{
Methylation imputation function
}
\usage{
MetIm()
}
\arguments{
\item{sequence}{matrix}
\item{microarray}{matrix}
\item{lambda}{numeric}
\item{cut}{numeric}
\item{cvfold}{numeric}
\item{use.mvalue}{logical}
\item{qc_frac}{numeric}

}
\details{
nothing
}
\value{
object
}
\references{
nothing
}
\author{
Degui Zhi
}
\note{
Require packages as following: "plotrix", "coxme", "kinship2", "nlme","MASS"
}
\seealso{
nothing
}
\examples{
library(plotrix)
library(coxme)
library(kinship2)
library(nlme)
library(MASS)
library(MetImpute)
data(package="MetImpute")
data(dat1)
data(dat2)
qc_frac = 1-5/dim(dat2)[2]
Seq <- MetIm(sequence = t(dat1), microarray = t(dat2), lambda=0.3, cut=10, cvfold=0, use.mvalue = T, qc_frac = qc_frac)
}
\keyword{ MethyImpute }
