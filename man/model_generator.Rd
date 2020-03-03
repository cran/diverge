\name{model_generator}
\alias{model_generator}
\title{
Generate parameter~gradient linear models
}
\description{
Generate a set of parameter~gradient linear models over the domain of a continuous gradient
}
\usage{
model_generator(domain, par_low, par_high)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{domain}{
Vector of length 2 defining the low and high ends of the gradient domain. Essentially identical to the 'xlim' argument in plotting functions. E.g. domain = c(0,60) may be suitable when investigating a latitudinal gradient.  
}
  \item{par_low}{
Parameter value or vector of parameter values at lowest end of domain 
}
  \item{par_high}{
Parameter value or vector of parameter values at hightest end of domain. NOTE: linear models returned are restricted such that the parameter is never more extreme than the values given in par_high and par_low across the domain.
}
}
\details{
An internal utility function called by param_grid that may nonetheless have value to some users. This function generates a set of linear models for a parameter value over a continuous gradient. The slopes and intercepts generated are ultimately used as starting parameter values in the likelihood search function 'find_mle'.

model_generator works by taking a user-defined parameter value or a set of values at opposite ends of a user-defined domain and calculating the slope of each straight line between them. The slopes thus generated create lines such that parameter values across the domain are never more extreme than the most extreme values supplied by the user. This helps avoid issues like nonsensical negative parameters. The number of models returned = length(par_low)*length(par_high). 
}
\value{
Returns a matrix of slope and intercept values where each row is the slope and intercept of one linear model.
}
\author{
Sean A.S. Anderson
}
\examples{
# Generate and plot 100 linear models for a parameter across a latitudinal gradient. 
# Param values must be between 1 and 10 at the equator and between 0 and 9 at 60 degrees.

mods = model_generator(domain=c(0,60), par_low=seq(1, 10), par_high=seq(0,9))
curve(expr=(y=mods[1,1]*x + mods[1,2]), ylim=c(0,10), xlim=c(0,60), ylab="Parmam Value",
  xlab="Abs. Latitude")
for(i in 2:nrow(mods)) curve(expr=(y=mods[i,1]*x + mods[i,2]), add=TRUE)
}