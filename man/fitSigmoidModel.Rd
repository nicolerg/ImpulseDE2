% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/srcImpulseDE2_fitSigmoid.R
\name{fitSigmoidModel}
\alias{fitSigmoidModel}
\title{Fit a sigmoidal model to data of a gene}
\usage{
fitSigmoidModel(vecSigmoidParamGuess, vecCounts, scaDisp, vecSizeFactors,
  lsvecidxBatch, vecTimepointsUnique, vecidxTimepoint, MAXIT = 1000,
  RELTOL = 10^(-8), trace = 0, REPORT = 10)
}
\arguments{
\item{vecSigmoidParamGuess}{(numeric vector length 4)
\{beta, h0, h1, t\}
Up model initialisations of sigmoidal model parameters.}

\item{vecCounts}{(numeric vector number of samples)
Read count data.}

\item{scaDisp}{(scalar) Gene-wise 
negative binomial dispersion hyper-parameter.}

\item{vecSizeFactors}{(numeric vector number of samples) 
Model scaling factors for each sample which take
sequencing depth into account (size factors).}

\item{lsvecidxBatch}{(list length number of confounding variables)
List of index vectors. 
One vector per confounding variable.
Each vector has one entry per sample with the index batch
within the given confounding variable of the given sample.
Batches are enumerated from 1 to number of batches.}

\item{vecTimepointsUnique}{(numeric vector length number of unique time points)
Unique time points of set of time points of given samples.}

\item{vecidxTimepoint}{(index vector length number of samples)
Index of of time point assigned to each sample in vector
vecTimepointsUnique.}

\item{MAXIT}{(scalar) [Default 1000] 
Maximum number of BFGS iterations for model fitting with \link{optim}.}

\item{RELTOL}{(scalar) [Default 10^(-8)]
Maximum relative change in loglikelihood to reach convergence in
numerical optimisation by BFGS in \link{optim}.}

\item{trace}{(scalar) [Defaul 0]
Reporting parameter of \link{optim}.}

\item{REPORT}{(scalar) [Default 10]
Reporting parameter of \link{optim}.}
}
\value{
(list) List of sigmoid fit parameters and results.
\itemize{
\item vecSigmoidParam (numeric vector length 4)
\{beta, h0, h1, t\}
Maximum likelihood estimators of sigmoidal model parameters.
\item vecSigmoidValue (numeric vector length number of time points)
Values of sigmoid model fit at time points used for fit.
\item lsvecBatchFactors (list length number of confounders)
List of vectors of scalar batch correction factors for each sample.
These are also maximum likelihood estimators.
NULL if no confounders given.
\item scaDispParam (scalar) Dispersion parameter estimate
used in fitting (hyper-parameter).
\item scaLL (scalar) Loglikelihood of data under maximum likelihood
estimator model.
\item scaConvergence (scalar) 
Convergence status of optim for sigmoid model.
}
}
\description{
[Model fitting function hierarchy: 3 out of 3]
This tertiary fitting wrapper performs sigmoidal model fitting:
This function executes numerical optimisaiton and error-handling
thereof.
}
\seealso{
Called by \link{fitSigmoidGene} to fit sigmoidal
model to samples of one condition and one gene.
Calls sigmoidal model cost function 
\link{evalLogLikSigmoid_comp} within \link{optim}.
}
\author{
David Sebastian Fischer
}
