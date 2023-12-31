% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/codes_from_literature.R
\name{law.func}
\alias{law.func}
\title{Implementing LAWS procedure.}
\usage{
law.func(pvs, pis, q)
}
\arguments{
\item{pvs}{p-values.}

\item{pis}{estimation of non-null propability pi.}

\item{q}{desired FDR level.}
}
\value{
A list with the elements
\item{nr}{number of rejections}
\item{th}{the threshold for rejection of the procedure}
\item{de}{the rejections}
}
\description{
Implementing LAWS procedure.
}
