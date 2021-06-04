\name{likelihood_func}
\alias{likelihood_func}
\title{
Functions for negative log likelihood
}
\description{
Calculate negative log likelihood for one of fourteen trait divergence models given a dataset and user-defined parameters. 
}
\details{
This function will typically will not be called directly by the user. Typing 'likelihood_func' may, however, be of interest to those wishing see/verify how we coded up the likelihood functions for each evolutionary model. The model returns the negative log-likelihood of a chosen model based on a dataset and a set of parameters. Likelihood_func is called by 'find_mle', which uses the non-linear search algorithm 'nlminb' to minimize its output. 
}