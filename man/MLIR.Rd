% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLIR.R
\name{MLIR}
\alias{MLIR}
\title{MLIR}
\usage{
MLIR(Mplus_Output, endogenous, exogenous, Interaction_factor, interaction_code)
}
\arguments{
\item{Mplus_Output}{the mplus output of the SEM model with one or more two-way interactions}

\item{endogenous}{the endogenous latent factor}

\item{exogenous}{the exogenous latent factors that are regressed on the endogenous latent factor}

\item{Interaction_factor}{the latent factors that are used to define the two-way interactions}

\item{interaction_code}{interactions that are defined by the XWITH function in the mplus}
}
\description{
R square of one or multiple two-way interactions between latent factors in structural equation model
}
\author{
{Lu Qin, Howard University, \email{lu.qin@howard.edu}}

{Jihong Zhang, University of Iowa, \email{jihong-zhang@uiowa.edu}}
}
