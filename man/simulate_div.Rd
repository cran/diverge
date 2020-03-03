\name{simulate_div}
\alias{simulate_div}
\title{
Simulate trait divergence 
}
\description{
Simulate replicate sets of lineage pair trait divergence under one of ten evolutionary models.
}
\usage{
simulate_div(model, pars, ages, GRAD = NULL, cats=NULL, breakpoint = NULL, N = 1)
}
\arguments{
  \item{model}{
Character string defining one of 11 models of trait divergence (options: "BM_null", "BM_linear", "OU_null", "OU_linear", "DA_null", "DA_linear", "DA_cat", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear"). See model_select for model descriptions.
}
  \item{pars}{
Vector of model parameters to use in sumulation. Parameters must be in the correct order in the vector for each model. The parameter orders are \cr
BM_null: sig2, \cr
BM_linear: sig2_slope, sig2_int \cr
OU_null: alpha, sig2 \cr
OU_linear: alpha_int, alpha_slope, sig2 \cr
DA_null: alpha, sig2, psi \cr
DA_linear: alpha, sig2, psi_slope, psi_int \cr
DA_cat: alpha, sig2, psi1, psi2, psi3(optional) \cr
DA_wt: alpha, sig2, psi1, psi2, wt \cr
DA_bp: alpha, sig2, psi1, psi2 \cr
DA_wt_linear: alpha, sig2, psi1_slope, psi1_int, psi2_slope, psi2_int, wt \cr
DA_bp_linear: alpha, sig2, psi1_slope, psi1_int, psi2_slope, psi2_int
}
  \item{ages}{
Vector of 'ages' (i.e. times since divergence) of lineage pairs for which to simulate trait divergence.
}
 \item{GRAD}{
Vector containing the gradient position of each pair to be simulated. This is the value of a continuous variable such as latitude or body size across which parameters will vary. Required for all models with 'linear' suffix.
}
 \item{cats}{
Vector containing the category code (0, 1, or 2) for each pair (see model_select help page for details).
}
  \item{breakpoint}{
Vector of breakpoint times (times after 0 at which a discrete shift occurs in psi) for each lineage pair to be simulated. Required for all "bp" models. Must align with "ages" vector such that ages[i] and breakpoint[i] correspond to the same lineage pair. When using bp models, datasets can and should be simulated with a mix of pairs that have and have not experienced a shift in psi (see model_select for details on bp models). Set bp=0 for pairs that experience no shift.
}
  \item{N}{
Number of datasets to simulate. 
}
}
\details{
Simulates replicate sets of trait divegence under one of ten evolutionary models given a user-defined set of lineage pair ages. Since all evolutionary models are stochastic process models, replicate datasets simulated under the same model parameters are not identical.
}
\value{
Returns a vector (if N=1) or a list of vectors (if N>1) for replicate sets of lineage-pair trait divergences.
}
\author{
Sean A.S. Anderson
}
\examples{

## Simulate datasets of trait divergence for 150 lineage pairs 
# under different evolutionary models

# Define vector of 150 lineage pair ages
ages = rep(c(0.5, 1, 1.5, 2, 3, 8), 25)

## simulate a dataset of trait divergence under a single-process Brownian Motion model
sig2 = 0.2
sis_div_BM = simulate_div(model="BM_null", ages=ages, pars=c(sig2))
summary(sis_div_BM)
sis_div_BM

## simulate 1000 divergence datasets under the null model of divergent adaptation
sig2 = 0.2
alpha = 0.8
psi = 0.8
sis_div_DA = simulate_div(model="DA_null", ages=ages, pars=c(alpha, sig2, psi), N=1000)
length(sis_div_DA)
sis_div_DA[[1]]

## Simulate 1000 divergence datasets under the breakpoint model
# note: a breakpoint vector must be created
# pairs that experience no epoch shift (no shift in the psi parameter) are assigned bp=0
# for all two-epoch pairs, the breakpoint time must be lower than the age of the pair
# here we make half the dataset into two-epoch pairs, half into one-epoch pairs
# we set the breakpoint time equal to half the age of each two-epoch pair
bp = c(ages[1:75]/2, rep(0, 75))

# run the simulation
sig2 = 0.2
alpha = 0.8
psi1 = 0.3
psi2 = 0.9
sis_div_bp = simulate_div(model="DA_bp", ages=ages, breakpoint=bp, 
  pars=c(alpha, sig2, psi1, psi2), N=1000)
length(sis_div_bp)
sis_div_bp[[1]]
}
