\name{param_grid}
\alias{param_grid}
\title{
Generate parameter grid
}
\description{
Generate a grid of starting parameter values for each model of trait divergence
}
\usage{
param_grid(model, domain = NULL, ncats=NULL)
}
\arguments{
  \item{model}{
Character string defining one of ten models of trait divergence (options: "BM_null", "BM_linear", "OU_null", "OU_linear", "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear"). See find_mle for model descriptions.
}
  \item{domain}{
Vector of length 2 defining the low and high ends of the gradient domain. Essentially identical to the 'xlim' argument in plotting functions. Required for models with 'linear' suffix.
}
  \item{ncats}{
A number (either 2 or 3) indicating the number of categories in a DA_cat model.
}
}
\details{
Primarily a utility function but might be useful in some other cases. Non-linear optimizers can often get stuck on local optima when finding the maximum likelihood parameter set, especially when calculating likelihoods with complex models. find_mle solves this problem by feeding the optimizer a grid of parameter values from which to launch its algorithm. While users can determine their own starting parameters, default parameter grids in model_select and find_mle are calculated with this function. To see/measure/assess the default starting parameter grid for a function of interest, users can use this function directly. 
}
\value{
Returns a matrix of starting parameter values. Each olumn contains different values for one parameter and each row is a unique parameter combination in the correct order for likelihood estimation. 
}
\author{
Sean A.S. Anderson and Jason T. Weir
}
\examples{
# Call the default parameter grid for the "DA_linear" model 
# asssume we are testing for a latitudinal gradient over 0-60 degrees.
par_grd = param_grid(model="DA_linear", domain=c(0,60))
dim(par_grd)
head(par_grd)

# Call the default parameter grid for a 3-category "DA_cat" model.
par_grd = param_grid(model="DA_cat", ncat=3)
dim(par_grd)
head(par_grd)
}
