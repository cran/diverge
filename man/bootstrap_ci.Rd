\name{bootstrap_ci}
\alias{bootstrap_ci}
\title{
Estimate confidence intervals using bootstrap
}
\description{
Calculates 95 percent confidence intervals of parameter estimates through bootstrap resampling
}
\usage{
bootstrap_ci(div, ages, GRAD = NULL, cats=NULL, breakpoint = NULL, domain = NULL, 
model, N, parallel = FALSE, cores = NULL, starting = NULL)
}
\arguments{
  \item{div}{
Vector of trait divergences for a set of lineage pairs. Calculated for each pair as abs(trait_val_lineage_2 - trait_val_lineage1). Raw values (i.e. not absolute values) can also be used but must be noted by the user in the argument 'absolute'.  
}
  \item{ages}{
Vector containing the age (i.e. estimated time since divergence) for each pair in the dataset. IMPORTANT: div, ages, GRAD, and breakpoint vectors must be aligned such that div[i] age[i] grad[i] and breakpoint[i] represent values for the same lineage pair.
}
  \item{GRAD}{
Vector containing the gradient position of each pair. This is the value of a continuous variable such as latitude or body size across which parameters are hypothesized to vary. Required for models with 'linear' suffix.
}
  \item{cats}{
Vector containing the category code (0, 1, or 2) for each pair (see package help page for details).
}
  \item{breakpoint}{
Vector of breakpoint times for each pair in the dataset. These are the times AFTER divergence at which a shift occurs in the psi parameter of DA model. Required for DA_bp and DA_bp_linear. See find_mle for details on how to calculate.
}
  \item{domain}{
Vector of length 2 defining the low and high ends of the gradient domain. Essentially identical to the 'xlim' argument in plotting functions. Required for models with 'linear' suffix.
}
  \item{model}{
Character string defining one of ten models of trait divergence (options: "BM_null", "BM_linear", "OU_null", "OU_linear", "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear"). See find_mle for model descriptions.
}
  \item{N}{
Number of bootstrap replicates to create (use integers only).
}
  \item{parallel}{
Logical indicating whether likelihood searches should be conducted in parallel across multiple cores. Not available on windows machines. Defaults to FALSE.
}
  \item{cores}{
If parallel=TRUE, the number of cores on which to run the function. Defaults to all virtual cores.
}
  \item{starting}{
Optional matrix of customized parameter values from which to launch likelihood searches. Must match the structure required for the chosen model. See the model descriptions in find_mle for details. If starting=NULL, default starting values are used. 
}
}
\details{
WARNING: estimating parameters for a large number of bootstrap datasets is computationally intensive. Run times can be exceedingly high for the breakpoint models -- in particular DA_bp and especially DA_bp_linear, DA_wt, DA_wt_linear. When running breakpoint models, we strongly recommend using the parallel option where resources allow, ideally on a multi-core server. 
}
\value{
A named list of length 4 containing [[1]] the model chosen, [[2]] all parameter estimates, [[3]], the number of replicates, and [[4]] summary statistics for parameter estimation.
}
\references{
Efron, B. and Tibshirani, R. (1986). The Boostrap Method for standard errors, confidence intervals, and other measures of statistical accuracy. Statistical Science, Vol 1., No. 1, pp 1-35.
}
\author{
Sean A.S. Anderson and Jason T. Weir
}
\examples{
## Estimate confidence intervals for sigma^2 under a single-process BM model 
## 10 bootstrap replicates

# simulate dataset
ages = rep(c(0.5, 1, 1.5, 2, 3, 8), 25)
sig2 = 0.2
sis_div = simulate_div(model="BM_null", ages=ages, pars=sig2)

# Run bootstrap_ci
N=100
res = bootstrap_ci(div=sis_div, ages=ages, model=("BM_null"), N=N)
res

\donttest{
## Estimate confidence intervals for all parameters under a DA_linear model 
## 10 bootstrap replicates. 

# simulate dataset under DA_linear 
# pairs are evenly distributed across a 0-60 degree latitudinal gradient
ages = rep(c(0.5, 1, 1.5, 2, 3, 8), 25)
grad_cats = rep(c(0, 15, 30, 45, 60), 30)
grad=c(rep(grad_cats[1], 30), rep(grad_cats[2],30), rep(grad_cats[3],30), 
  rep(grad_cats[4],30), rep(grad_cats[5],30))
alpha = 0.8
sig2 = 0.2
psi_sl = -0.01
psi_int = 2
sis_div = simulate_div(model="DA_linear", ages=ages, pars=c(alpha, sig2, psi_sl, psi_int), 
  GRAD=grad)

# Run bootstrap_ci
N = 10
res = bootstrap_ci(div=sis_div, ages=ages, GRAD=grad, domain=c(0,60), 
model=("DA_linear"), N=N)
res

## Estimate confidence intervals for psi1 and psi2 under a 2-category DA_cat 
## model given 10 bootstrap replicates.
ages = rep(c(0.5, 1, 1.5, 2, 3, 8), 25)
grad_cats = rep(c(0, 15, 30, 45, 60), 30)
cats = c(rep(0, 75), rep(1, 75))
alpha = 0.8
sig2 = 0.2
psi1 = 0.5
psi2 = 1
N = 2
sis_div = simulate_div(model="DA_cat", ages=ages, pars=c(alpha, sig2, psi1, psi2), cats=cats)
res = bootstrap_ci(div=sis_div, ages=ages, cats=cats, model=("DA_cat"), N=N, parallel=FALSE)
res
}
}