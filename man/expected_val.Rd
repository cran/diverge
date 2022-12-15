\name{expected_val}
\alias{expected_val}
\title{
Expected value through time
}
\description{
Calculate and plot expected value of trait divergence through time
}
\usage{
expected_val(model, sig2, alpha = NULL, psi = NULL, time_span = c(0, 10), 
  quantile = FALSE, plot = FALSE, labels = TRUE, bm_col = "darkgoldenrod1", 
  ou_col = "firebrick2", da_col = "navy", exval_lwd = 5, ylim = NULL, ...)
}
\arguments{
  \item{model}{
Name of the model whose expected values are to be calculated. Must be one of "BM_null", "OU_null", "DA_null", "DA_OU", "DA_BM", or "OU_BM".
}
  \item{sig2}{
Value of the sigma-squared parameter. Required for all models.
}
  \item{alpha}{
Value of the alpha parameter. Required for all OU and DA models.
}
  \item{psi}{
Value of psi parameter. Required for DA models. NOTE: if using a breakpoint or wait time model, this is the psi for the first evolutionary epoch
}
  \item{time_span}{
The time after initial divergence at which to calculate the expected value. If time_span is a vector lf length 2, the expectation is calculated continuously between the two times. If time_span is any other length, the expected value is calculated only at the given times. Defaults to c(0,10)
}
  \item{quantile}{
Logical indicating whether to additionally calculate (and optionally plot) the expected quantiles (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95) of trait divergence. Defaults to FALSE.
}
  \item{plot}{
Logical indicating whether to plot the expected values through time. If quantile=TRUE, then quantile curves are also plotted. Defaults to FALSE.
}
  \item{labels}{
Logical indicating whether x and y axes are to be labelled.
}
  \item{bm_col}{
Color of the line (and shaded if quantile=TRUE) area showing expected value through time under a BM model. 
}
  \item{ou_col}{
Color of the line (and shaded if quantile=TRUE) area showing expected value through time under an OU model. 
}
  \item{da_col}{
Color of the line (and shaded if quantile=TRUE) area showing expected value through time under an OU model. 
}
  \item{exval_lwd}{
Width in which the line showing the expected value through time should be plotted. Defaults to 5.
}
  \item{ylim}{
Custom y-axis for graph in the form 'ylim = c(miny, maxy)'
}
  \item{...}{
Additional arguments to be passed to the 'plot' function for custom visualization.
}
}
\details{
Calculates and plots the expected value of continuous trait divergence under different models of trait evolution. Time begins when lineages initially depart form a common ancestor. If users provide upper and lower bounds for time_span, then that span is divided into 10,001 time steps at which expectations are calculated to approximate a continuous tracking of expectation through time. Expectations for the 2.5th-97.5th quantile can also be estimated and 95 percent confidence intervals plotted to visualize how the distribution of trait divergence changes. 

IMPORTANT: Not all models can be selected for this function. The expected value of divergence under one of the "linear" models, for instance, would be different for each position on a continuous gradient and thus can't be calculated or plotted in this way. Acceptable models are ("BM_null", "OU_null", "DA_null", "DA_OU", "DA_BM", "OU_BM").
}
\value{
Returns either (1) a matrix of 2 columns in which col1 = time and col2 = the expected trait divergence at that time, or (2) a matrix of 12 columns that additionally contains the expectation for the 10th-95th quantile of trait divergence at each time (if quantiles=TRUE is indicated in the function call). The number of rows of the output matrix equals the length of time_span EXCEPT if length(time_span) == 2, in which case the output matrix has 10,001 rows.
}
\author{
Sean A.S. Anderson and Jason T. Weir
}
\examples{
# Ex. 1. Calculate the expected value of divergence after 5my under DA_null 
sig2 = 0.2
alpha = 0.8
psi = 0.9
exval = expected_val(model="DA_null", sig2=sig2, alpha=alpha, psi=psi, time_span=5)

# Ex. 2. calculate and plot expected trait divergence and expected quantiles through time over 8my 
# under a DA_null model
sig2 = 0.2
alpha = 0.8
psi = 0.3
\donttest{
exval = expected_val(model="DA_null", sig2=sig2, alpha=alpha, psi=psi, time_span=c(0,8),
  quantile=TRUE, plot=TRUE)
}
  
# Ex. 3. same as above but with customized graphical parameters
sig2 = 0.2
alpha = 0.8
psi = 0.3
\donttest{
exval = expected_val("DA_null", sig2=sig2, alpha=alpha, psi=psi, time_span=c(0,8), quantile=TRUE, 
	plot=TRUE, da_col="green", exval_lwd = 4, ylim = c(0, 3), axes=FALSE, labels=FALSE)
box()
axis(1, labels=NA)
axis(1, lwd = 0, line = -0.6)
axis(2, labels = NA)
axis(2, lwd = 0, line = -0.6)
title(line = 1.9, xlab = "Custom X Axis Title")
title(line = 1.9, ylab = "Custom Y AXis Title")
}

# Ex. 4. calculate and plot expected trait divergence and expected quantiles through time over 8my 
# under a DA-OU mixture model
sig2 = 0.05
alpha = 0.8
psi = 1
\donttest{
exval = expected_val(model="DA_OU", sig2=sig2, alpha=alpha, psi=psi, time_span=c(0,8),
  quantile=TRUE, plot=TRUE)
# add a legend
legend("topleft", legend=c("ExVal DA", "ExVal OU"), lwd=2, col=c("navy", "firebrick2"), 
    bty="n", cex=0.8)
  }
}