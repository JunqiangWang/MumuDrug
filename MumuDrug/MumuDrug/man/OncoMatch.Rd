% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OncoMatchExtended.R
\name{OncoMatch}
\alias{OncoMatch}
\alias{oncoTreat}
\title{DarwinOncoMatch}
\usage{
OncoMatch(
  vpmat,
  vpsig,
  nn = 50,
  adaptive = FALSE,
  modsig = NULL,
  filter = TRUE
)
}
\arguments{
\item{vpmat}{Viper matrix for the phenotypes}

\item{vpsig}{Viper matrix for the signatures or perturbation object: list of 2 components vpmat (viper matrix) and protweight (named vector of weights for the rows)}

\item{nn}{Number indicating the number of MR to consider for the enrichment analysis}

\item{modsig}{Optional protein activity signature for the model. If provided, only MRs conserved in the model are used for the analysis}
}
\value{
Matrix of Normalized Enrichment Scores
}
\description{
This function performs oncoMatch analysis
}
