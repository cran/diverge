\name{random_walks}
\alias{random_walks}
\title{
Simulate and plot trait evolution through time
}
\description{
Simulate and plot continuous trait values as they evolve through time in two diverging lineages under one of four general models.
}
\usage{
random_walks(model, TIME, nsim, sig2, alpha = NULL, psi = NULL, psi2 = NULL, 
  wt = NULL, theta = 0, centre = 0, steps = 100, plot = TRUE, col = c("black", "red"),
  labels = TRUE, ylim = NULL, ... )
}

\arguments{
  \item{model}{
Character vector naming the model under which traits are evolving (one of "BM_null", "OU_null", "DA_null", "DA_wt")
}
  \item{TIME}{
Numeric value in units of millions of years over which to plot trait evolution
}
  \item{nsim}{
Number of replicate walks to simulate. If N=40, then 20 replicates of lineages i and 20 replicates of lineage j will be simulated. 
}/
  \item{sig2}{
Numeric value of sigma squared parameter. Required for all models.
}
  \item{alpha}{
Numeric value of alpha parameter for OU and DA models
}
  \item{psi}{
Numeric value of psi parameter for DA models. If using a breakpoint or wait time model, this is the psi for the first evolutionary epoch
}
  \item{psi2}{
Numeric value of psi for the second evolutionary epoch if using DA_wt
}
  \item{wt}{
Numeric value for wait time to a shift in the adaptive regime (aka epoch shift)
}
  \item{theta}{
Optimum trait value to which both lineages are drawn in a single-process OU model. Defaults to zero. 
}
  \item{centre}{
For DA models, the centre value between the optima of two lineages. Defaults to zero.
}
  \item{steps}{
Numeric indicating the number of steps per unit TIME over which to calculate simulated trait values - greater numbers of steps create more continuous-looking trait curves, but also larger image files
}
  \item{plot}{
Logical indicating whether to plot replicate walks.
}
  \item{col}{
Character vector containing the names of the two colours in which to plot replicate random walks for the two lineages. 
}
  \item{labels}{
Logical indicating whether x and y axes are to be labelled
}
  \item{ylim}{
Vector of length 2 indicating the max and min y values for the plot. Same argument as in plotting functions. Defaults so that trait curves fill the plot area
}
  \item{...}{
Additional parameters for the 'plot' function for custom visualization
}
}
\details{
This function is primarily a visualization aid that will help users intuit the processes under which lineages in a pair evolve after departing from an ancestor. In BM and OU models, lineages in a pair undergo identical evolutionary processes after separating and trait differentiation is generated by random deviations in evolutionary trajectory (i.e. differences are solely generated by the stochastic component of stochastic process models). In the DA models, lineages are pulled deterministically toward alternative adaptive optima, though stochastic noise is still present. The x-label in the default graph is the number of time steps. NOTE: To plot a single random process rather than two diverging lineages (possible for BM and OU models only), simply give the same colour for the lines of both lineages (e.g. col=c("black","black"))
}
\value{
The object returned depends on the model selected. For BM_null and OU_null, random_walks returns an N*M matrix of simulated trait values, where N is the number of replicate walks and M is the number of time steps. For DA_null and DA_bp, random_walks returns a list of length two, where the first element is a K*M matrix (where K=N/2) of simulated trait values for lineage "i" and the second element is a K*M matrix of simulated trait values for lineage "j".The number of rows of the output matrix is equal to the number of time steps.
}
\author{
Sean A.S. Anderson
}
\examples{
# plot replicate runs of two lineages evolving under Brownian motion
rw = random_walks(model="BM_null", TIME=5, steps=200, nsim=40, sig2 = 0.2, plot=TRUE)

# plot replicate runs of two lineages evolving under an OU process
rw = random_walks(model="OU_null", TIME=5, steps=200, nsim=40, sig2 = 0.2, alpha=0.8, 
plot=TRUE)

# same as above but with custom y axis
rw = random_walks(model="OU_null", TIME=5, steps=200, nsim=40, sig2 = 0.2, alpha=0.8, 
plot=TRUE, ylim=c(-2,2))

# plot replicate runs of two lineages evolving under an divergent selection
rw = random_walks(model="DA_null", TIME=5, steps=200, nsim=40, sig2 = 0.2, alpha=0.8, 
psi=1, plot=TRUE)

# plot replicate runs of two lineages evolving under an divergent selection 
# with a breakpoint in the strength of divergent selection at 2.5my
rw = random_walks(model="DA_wt", TIME=5, steps=200, nsim=40, sig2 = 0.2, alpha=0.8, 
psi=0, psi2=1.5, wt=2.5, plot=TRUE)

# play with colours
rw = random_walks(model="DA_null", TIME=5, steps=200, nsim=40, sig2 = 0.2, alpha=0.8, 
psi=1.5, plot=TRUE, col=c("purple","orange"), ylim=c(-2,2))

# use custom axes and show time in units of millions of years
rw = random_walks(model="DA_null", TIME=5, steps=200, nsim=40, sig2 = 0.2, alpha=0.8, 
psi=1.5, plot=TRUE, col=c("purple","orange"), labels = FALSE, ylim=c(-2,2), axes=FALSE)
box()
axis(1, labels=NA)
axis(1, lwd = 0, line = -0.6, at = seq(0, 1000, 200), labels = seq(0, 10, 2))
axis(2, labels = NA)
axis(2, lwd = 0, line = -0.6)
title(line = 1.9, xlab = "Time (Ma)")
title(line = 1.9, ylab = "Trait Value (mm)")

# plot a single BM process
rw = random_walks(model="BM_null", TIME=5, steps=200, nsim=40, sig2 = 0.2, plot=TRUE, 
col=c("black","black"))

# plot a single OU process
rw = random_walks(model="OU_null", TIME=5, steps=200, nsim=40, sig2 = 0.2, alpha=0.9, 
plot=TRUE, ylim=c(-2,2), col=c("black","black"))
}
