% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OncoMatchExtended.R
\name{oncoMatch_Darwin}
\alias{oncoMatch_Darwin}
\title{DarwinOncoMatch}
\usage{
oncoMatch_Darwin(
  vpmat,
  vpsig,
  nn = 50,
  balanced = TRUE,
  modsig = NULL,
  reverse = FALSE
)
}
\arguments{
\item{vpmat}{Viper matrix for the phenotypes}

\item{vpsig}{Viper matrix for the signatures or perturbation object: list of 2 components vpmat (viper matrix) and protweight (named vector of weights for the rows)}

\item{nn}{Number indicating the number of MR to consider for the enrichment analysis}

\item{balanced}{Logical, whether to use equal number of activated and inactivated MRs}

\item{modsig}{Optional protein activity signature for the model. If provided, only MRs conserved in the model are used for the analysis}

\item{reverse}{Logical, whether the enrichment of the drug MRs should e computed on the phenotype signature instead of the phenotype MRs on the drug signature}
}
\value{
Matrix of Normalized Enrichment Scores
}
\description{
This function performs oncoMatch analysis
}
