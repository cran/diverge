\name{find_mle}
\alias{find_mle}
\title{
Find Maximum Likelihood Estimates 
}
\description{
Find maximum likelihood parameter estimates for an evolutionary model given data
}
\usage{
find_mle(model, p_starting = NULL, div, ages, me1=NULL, me2=NULL, GRAD = NULL, cats, 
  breakpoint = NULL, domain = NULL, absolute = TRUE, parallel = FALSE, cores = NULL)
}
\arguments{
  \item{model}{
Character string defining one of 11 models of trait divergence (options: "BM_null", "BM_linear", "OU_null", "OU_linear", "DA_null", "DA_linear", "DA_cat", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear"). See Details for info on each model.
}
  \item{p_starting}{
Optional matrix of customized parameter values from which to launch likelihood searches. Must match the structure required for the chosen model. See model descriptions in the model_select man page for details on providing a custom starting matrix.
}
  \item{div}{
Vector of trait divergences for a set of lineage pairs. Calculated for each pair as abs(trait_val_lineage_2 - trait_val_lineage1). Raw values (i.e. not absolute values) can also be used but must be noted by the user in the argument 'absolute'.  
}
  \item{ages}{
Vector containing the age (i.e. estimated time since divergence) for each pair in the dataset. IMPORTANT: div, ages, GRAD, and breakpoint vectors must be aligned such that div[i] age[i] grad[i] and breakpoint[i] represent values for the same lineage pair.
}
  \item{me1}{
Vector containing the measurement error (standard error of mean) for species 1s of each pair (for one species = variance/No.measurements)
}
  \item{me2}{
Vector containing the measurement error (standard error of mean) for species 2s of each pair (for one species = variance/No.measurements)
}
  \item{GRAD}{
Vector containing the gradient position of each pair. This is the value of a continuous variable such as latitude or body size across which parameters are hypothesized to vary. Required for models with 'linear' suffix.
}
  \item{cats}{
Vector containing the category code (0, 1, or 2) for each pair
}
  \item{breakpoint}{
Vector of breakpoint times for each pair in the dataset. These are the times AFTER divergence at which a shift occurs in the psi parameter of DA model. Required for DA_bp and DA_bp_linear. See find_mle for details on how to calculate.
}
  \item{domain}{
Vector of length 2 defining the low and high ends of the gradient domain. Essentially identical to the 'xlim' argument in plotting functions. Required for models with 'linear' suffix.
}
  \item{absolute}{
Logical indicating whether 'div' represents absolute value of trait divergence or the raw divergence values.
}
  \item{parallel}{
Logical indicating whether likelihood searches should be conducted in parallel across multiple cores. Not available on windows machines. Defaults to FALSE.
}
  \item{cores}{
If parallel=TRUE, the number of cores on which to run the function. Defaults to all virtual cores.
}
}
\details{
This function can be used to find the maximum likelihood parameter estimates for one of ten evolutionary models given a lineage-pair dataset.'find_mle' uses 'nlminb' to optimize the likelihood function of the selected model. Since optimization algorithms can often get stuck on local optima, find_mle works by feeding the optimizer a large set of starting parameter values and ranking the likelihoods of estimates that are returned by each. Many sets of starting parameters should converge on the MLE, and find_mle returns the output of nlminb given one of these optimal starting sets. 

RUN TIME CONSIDERATIONS \cr
The set of starting parameter values grows rapidly with parameter number, so calculating MLEs in the most complex models (esp DA_wt, DA_wt_linear) can take a few minutes. We therefore recommend using the parallel=TRUE option when using any model more complex than DA_null. Users can define their own starting parameter sets but should keep in mind this tradeoff between run time and the breadth of parameter space through which to search. 
}
\value{
A list containing maximum likelihood parameter estimates, the maximumum likelihood, and likelihood search diagnostics. This is the output of nlminb as returned from an optimal starting parameter set.
}
\author{
Sean A.S. Anderson and Jason T. Weir
}
\examples{

## Find the maximum likelihood and parameters for a DA_linear model 
# assume an elevational gradient from 0-1000m

# simulate a dataset
ages = rep(c(0.5, 1, 1.5, 2, 3, 8), 25)
grad_cats = rep(c(0, 250, 500, 750, 1000), 30)
grad=c(rep(grad_cats[1], 30), rep(grad_cats[2],30), rep(grad_cats[3],30), 
  rep(grad_cats[4],30), rep(grad_cats[5],30))
alpha = 0.8
sig2 = 0.2
psi_sl = -0.01
psi_int = 2
sis_div = simulate_div(model="DA_linear", ages=ages, 
  pars=c(alpha, sig2, psi_sl, psi_int), GRAD=grad)

# find max. likelihood and estimate parameters
\donttest{res = find_mle(model="DA_linear", div=sis_div, ages=ages, GRAD=grad, domain=c(0,1000))
	}
}