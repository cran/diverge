\name{model_error_rate}
\alias{model_error_rate}
\title{
Calculate error rate in model selection
}
\description{
Estimate error rates in model selection using replicate datasets simulated under a known model
}
\usage{
model_error_rate(datasets, ages, sim_model, alternatives, type, me1=NULL, me2=NULL, 
  GRAD = NULL, cats=NULL, breakpoint = NULL, domain = NULL, threshold = 0, 
  parallel = FALSE, cores = NULL)
}
\arguments{
  \item{datasets}{
List containing replicate datasets of trait divergence simulated under a single model
}
  \item{ages}{
Vector of lineage-pair ages for which trait divergence was simulated to create the 'datasets'
}
  \item{sim_model}{
Character vector naming the model under which trait divergence was simulated (options: "BM_null", "BM_linear", "OU_null", "OU_linear", "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear"). See find_mle for info on each model.
}
  \item{alternatives}{
Character vector naming alternative models to be included in the model comparison
}
  \item{type}{
Numeric "1" or "2". Indicates the whether test is for type 1 error rate (i.e. false rejection of true null) or type 2 error rate (i.e. failure to reject incorrect null).
}
  \item{me1}{
Vector containing the measurement error (standard error of mean) for species 1s of each pair (for one species = variance/No.measurements)
}
  \item{me2}{
Vector containing the measurement error (standard error of mean) for species 2s of each pair (for one species = variance/No.measurements)
}
  \item{GRAD}{
Vector of gradient values. Required when any model with linear gradients is included in the tests (see find_mle). Must match length of 'ages' vector. 
}
  \item{cats}{
Vector containing the category code (0, 1, or 2) for each pair (see package help page for details). 
}
  \item{breakpoint}{
Vector of breakpoints. Required when DA_bp or DA_bp_linear is included in the tests. Must match length of 'ages' vector.
}
  \item{domain}{
Vector of length 2 defining the low and high ends of the gradient domain. Essentially identical to the 'xlim' argument in plotting functions. Required when any models with 'linear' suffix is included in the tests.
}
  \item{threshold}{
Numeric value for the delta_AICc threshold corresponding to an error in model selection (see details). Defaults to zero. 
}
  \item{parallel}{
Logical indicating whether likelihood searches should be conducted in parallel across multiple cores. Not available on windows machines. Defaults to FALSE.
}
  \item{cores}{
If parallel=TRUE, the number of cores on which to run the function. Defaults to all virtual cores.
}
}
\details{
'model_error_rate' calculates the proportion of replicate runs of model_select that return an incorrect model based on a user-defined AICc threshold. A commonly-used threshold is 2, and the default threshold is zero. Since users must know which of the 11 evolutionary models produced the replicate sets of trait divergence, the function is typically used in conjunction with simulate_div. 

Users must designate whether type 1 or type 2 error rate is to be calculated via the 'type' argument. If type = 1, sim_model must be the simplest model against which more general alternatives are compared. If type = 2, sim_model must be a more complex model against which a simpler null alternative is compared. Type 1 error rate is calculated as the proportion of model_select runs for which delta_AICc of the sim model exceeds the threshold. Type 2 error rate is calculated as the proportion of model_select runs for which the difference between AICc of the sim_model and that of the best null model is greater than the threshold*(-1).

TIME CONSIDERATIONS \cr
Warning: this function runs model selection tests on replicate datasets. Model selection can be time intensive in its own right, especially when more complex models are involved. Running model selection trials across a large number of replicate datasets can therefore be extremely time intensive. For example, calculating error rate based on 1000 datasets when the models include DA_wt, DA_bp, DA_wt_linear, and DA_bp_linear will often require parallel runs on a multi-core server over several hours
}
\value{
The error rate in decimal form.
}
\author{
Sean A.S. Anderson
}
\examples{
\donttest{
## Calculate type I error rate for BM_null versus BM_linear and OU_linear
# simulate 10 replicate sets under BM_null
sig2 = 0.2
N = 10
ages = rep(c(0.5, 1, 1.5, 2, 3, 8), 25)
dsets = simulate_div(model="BM_null", ages=ages, pars=sig2, N=N)

# generate a continuous gradient (in this case from 0-60 degrees latitude)
grad_cats = rep(c(0, 15, 30, 45, 60), 30)
grad=c(rep(grad_cats[1], 30), rep(grad_cats[2],30), rep(grad_cats[3],30), 
  rep(grad_cats[4],30), rep(grad_cats[5],30))

# calculate type I error rate given the dataset and gradient
err_rateT1 = model_error_rate(datasets=dsets, ages=ages, sim_model="BM_null", 
  alternatives=c("BM_linear","OU_linear"), type=1, GRAD=grad, 
  domain=c(0,60), threshold=2)
err_rateT1 

## Calculate type II error rate of DA_linear v DA_null
# simulate 10 replicate sets under DA_linear
alpha = 0.8
sig2 = 0.2
psi_sl = -0.01
psi_int = 2
ages = rep(c(0.5, 1, 1.5, 2, 3, 8), 25)
grad_cats = rep(c(0, 15, 30, 45, 60), 30)
grad=c(rep(grad_cats[1], 30), rep(grad_cats[2],30), rep(grad_cats[3],30), 
  rep(grad_cats[4],30), rep(grad_cats[5],30))
dsets = simulate_div(model="DA_linear", ages=ages, pars=c(alpha, sig2, psi_sl, psi_int),
  N=10, GRAD=grad)

# calculate type II error rate
# NOTE: this can take ~1-2 mins
err_rateT2 = model_error_rate(datasets=dsets, ages=ages, sim_model="DA_linear", 
  alternatives="DA_null", type=2, GRAD=grad, domain=c(0,60), threshold=2)
err_rateT2 
  }
}
