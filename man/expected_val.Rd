\name{expected_val}
\alias{expected_val}
\title{
Expected value through time
}
\description{
Calculate and plot expected value of trait divergence through time
}
\usage{
expected_val(sig2, alpha = NULL, psi = NULL, psi2 = NULL, wt = NULL, 
  time_span = c(0, 10), quantile = FALSE, plot = FALSE, labels = TRUE, 
  exval_col = "red", quant_col = "black", exval_lwd = 1, quant_lwd = 1, ...)
}
\arguments{
  \item{sig2}{
Value of the sigma-squared parameter. Required for all models.
}
  \item{alpha}{
Value of the alpha parameter. Required for all OU and DA models.
}
  \item{psi}{
Value of psi parameter. Required for DA models. NOTE: if using a breakpoint or wait time model, this is the psi for the first evolutionary epoch
}
  \item{psi2}{
Value of psi for the second evolutionary epoch in breakpoint or wait time DA models.
}
  \item{wt}{
Value representing the shared wait time to a shift in the adaptive regime (aka an epoch shift)
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
  \item{exval_col}{
Color in which the line showing expected value through time should be plotted. Defaults to red.
}
  \item{quant_col}{
Color in which the lines showing the 10-95th quantiles through time should be plotted. Defaults to black.
}
  \item{exval_lwd}{
Line width in which the line showing the expected value through time should be plotted. Defaults to 1.
}
  \item{quant_lwd}{
Line width in which the lines showing the 10-95th quantiles through time should be plotted. Defaults to 1.
}
  \item{...}{
Additional arguments to be passed to the 'plot' function for custom visualization.
}
}
\details{
Calculates and plots the expected value of continuous trait divergence under different models of trait evolution. Time begins when lineages initially depart form a common ancestor. If users provide upper and lower bounds for time_span, then that span is divided into 10,001 time steps at which expectations are calculated to approximate a continuous tracking of expectation through time. Expectations for the 10th-95th quantile can also be estimated and plotted to visualize how the distribution of trait divergence changes. 

IMPORTANT: The typical suite of 11 evolutionary models does not explicitly apply for this function. The expected value of divergence under one of the "linear" models, for instance, would be different for each position on a continuous gradient. Instead, expected_val calculates the expected value of trait divergence under (1) a single BM process, (2) a single OU process, (3) a DA model with no epoch shift, and (4) DA model with a single epoch shift. These models are not called by name as arguments, but are inferred based on the parameters supplied. If only sig2 is supplied, the function will calculate expected value through time under a single-process BM model. If alpha and sig2 are supplied, divergence is calculated under an single-process OU model. If alpha, sig2, and psi are supplied, divergence is calculated under a constant DA model. If alpha, sig2, psi1, psi2, and wt are supplied, divergence is calculated under a DA model with an epoch shift at wt. 

}
\value{
Returns either (1) a matrix of 2 columns in which col1 = time and col2 = the expected trait divergence at that time, or (2) a matrix of 12 columns that additionally contains the expectation for the 10th-95th quantile of trait divergence at each time. The number of rows of the output matrix equals the length of time_span EXCEPT if length(time_span) == 2, in which case the output matrix has 10,001 rows.
}
\author{
Sean A.S. Anderson and Jason T. Weir
}
\examples{
# Calculate the expected value of divergence after 5my under DA_null 
sig2 = 0.2
alpha = 0.8
psi = 0.9
exval = expected_val(sig2=sig2, alpha=alpha, psi=psi, time_span=5)

# calculate and plot expected trait divergence through time over 10my under a BM_null process
sig2 = 0.3
exval = expected_val(sig2, quantile=FALSE, plot=TRUE)

# same as above but with customized graphical parameters
sig2 = 0.3
exval = expected_val(sig2, quantile=FALSE, plot=TRUE, labels=FALSE, exval_col="blue", 
  ylim = c(0, 6), exval_lwd=4,  axes=FALSE)
box()
axis(1, labels=NA)
axis(1, lwd = 0, line = -0.6)
axis(2, labels = NA)
axis(2, lwd = 0, line = -0.6)
title(line = 1.9, xlab = "Custom X Axis Title")
title(line = 1.9, ylab = "Custom Y AXis Title")

# calculate and plot expected trait divergence and expected quantiles through time over 8my 
# under a DA model with an epoch shift at 5my
sig2 = 0.2
alpha = 0.8
psi = 0.3
psi2 = 1.3
exval = expected_val(sig2=sig2, alpha=alpha, psi=psi, psi2=psi2, wt=5, time_span=c(0,8),
  quantile=TRUE, plot=TRUE)
  
# same as above but with customized graphical parameters
sig2 = 0.2
alpha = 0.8
psi = 0.3
psi2 = 1.3
exval = expected_val(sig2, alpha=alpha, psi=psi, psi2=psi2, wt=5, time_span=c(0,8), quantile=TRUE, 
	plot=TRUE, exval_col="green", exval_lwd = 4, ylim = c(0, 3), axes=FALSE, labels=FALSE)
box()
axis(1, labels=NA)
axis(1, lwd = 0, line = -0.6)
axis(2, labels = NA)
axis(2, lwd = 0, line = -0.6)
title(line = 1.9, xlab = "Custom X Axis Title")
title(line = 1.9, ylab = "Custom Y AXis Title")
}