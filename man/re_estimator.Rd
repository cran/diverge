\name{re_estimator}
\alias{re_estimator}
\title{
Estimate model parameters
}
\description{
Re-estimate known model parameters from replicate datasets simulated under a given model. Facilitates power analyses and the measure of accuracy and precision in parameter estimation. 
}
\usage{
re_estimator(model, div, ages, me1=NULL, me2=NULL, GRAD = NULL, cats=NULL, 
  breakpoint = NULL, domain = NULL, p_starting = NULL, parallel = FALSE, cores = NULL, 
  absolute = TRUE)
}
\arguments{
  \item{model}{
Character string defining one of ten models of trait divergence (options: "BM_null", "BM_linear", "OU_null", "OU_linear", "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear"). See find_mle for model descriptions.
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
Vector containing the category code (0, 1, or 2) for each pair (see model_select help page for details).
}
  \item{breakpoint}{
Vector of breakpoint times for each pair in the dataset. These are the times AFTER divergence at which a shift occurs in the psi parameter of DA model. Required for DA_bp and DA_bp_linear. See find_mle for details on how to calculate.
}
  \item{domain}{
Vector of length 2 defining the low and high ends of the gradient domain. Essentially identical to the 'xlim' argument in plotting functions. Required for models with 'linear' suffix.
}
  \item{p_starting}{
Optional matrix of customized parameter values from which to launch likelihood searches. Must match the structure required for the chosen model. See the model descriptions in find_mle for details. If p_starting=NULL, default starting values are used.
}
  \item{parallel}{
Logical indicating whether likelihood searches should be conducted in parallel across multiple cores. Not available on windows machines. Defaults to FALSE.
}
  \item{cores}{
If parallel=TRUE, the number of cores on which to run the function. Defaults to all virtual cores.
}
  \item{absolute}{
Logical indicating whether 'div' represents absolute value of trait divergence or the raw divergence values.
}
}
\details{
re_estimator is a wrapper to find_mle that returns just the parameter values. The argument 'div' can be a vector of trait divergences or a list containing several vectors of trait divergence. When 'div' is a list, the function can be used to generate a distribution of parameter values. These can in turn be used to evaluate precision and accuracy of likelihood searches when replicate datasets are simulated under a known model. 

TIME CONSIDERATIONS
Warning: re-estimating parametres from a large number of replicate datasets can be extremely time intensive when DA_wt, DA_bp, DA_wt_linear, or DA_bp_linear are the model in question. Re-estimating parameters to assess model performance using a reasonable number of replicates (~1000), from one of these four complex models, can take hours on a multi-core server. 
}
\value{
Returns a matrix of parameter estimates. Each column contains replicate estimates of one parameter. The matrix has one row if div is a vector and N rows if div is a list of length N.
}
\author{
Sean A.S. Anderson
}
\examples{
# simulate a dataset of trait divergence with 150 lineage pairs with a range of ages.
ages = rep(c(0.5, 1, 1.5, 2, 3, 8), 25)
alpha = 0.8
sig2 = 0.2
psi=0.8
N=10
sis_div = simulate_div(model="DA_null", ages=ages, pars=c(alpha, sig2, psi), N=N)

# estimate pars from one trait divergence dataset
par_es = re_estimator(model = "DA_null", div=sis_div[[1]], ages=ages)
par_es

# estimate pars from a list of trait divergence datasets and find the median
# of each estimate
\donttest{
par_es = re_estimator(model = "DA_null", div=sis_div, ages=ages)
head(par_es)
meds = apply(par_es, 2, median)
meds
}
}