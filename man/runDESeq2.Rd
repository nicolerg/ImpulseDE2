% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srcImpulseDE2_runDESeq2.R
\name{runDESeq2}
\alias{runDESeq2}
\title{Wrapper function for running DESeq2}
\usage{
runDESeq2(dfAnnotationProc, matCountDataProc, boolCaseCtrl, vecConfounders)
}
\arguments{
\item{dfAnnotationProc}{(data frame samples x covariates) 
{Sample, Condition, Time (numeric), TimeCateg (str)
(and confounding variables if given).}
Processed annotation table with covariates for each sample.}

\item{matCountDataProc}{(matrix genes x samples)
Read count data.}

\item{boolCaseCtrl}{(bool) 
Whether to perform case-control analysis. Does case-only
analysis if FALSE.}

\item{vecConfounders}{(vector of strings number of 
confounding variables)
Factors to correct for during batch correction. Have to 
supply dispersion factors if more than one is supplied.
Names refer to columns in dfAnnotationProc.}
}
\value{
(numeric vector length number of genes)
Dispersion parameter estimates for each gene.
In format of parameter size of \link{dnbinom}
which is 1/dispersion factor of DESeq2.
}
\description{
Run DESeq2 and extract dispersion parameter estimates.
Catch and remove dispersion outlier exception on samples
with zero-count observations.
}
\seealso{
Called by \link{runImpulseDE2}.
}
\author{
David Sebastian Fischer
}
