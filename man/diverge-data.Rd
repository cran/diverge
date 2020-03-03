\name{diverge-data}
\alias{all_birds}
\title{
Example Tree
}
\description{
A phylogenetic tree of birds of the world.
}
\usage{
data(all_birds)
}
\details{
all_birds is the maximum clade credibility tree of 9982 avian species (and 3301 sister pairs) from Pulido-Santacruz and Weir (2016). Branch lengths are time-calibrated in millions of years. 
}
\value{
A phylo object.
}
\references{
Pulido-Santacruz, P. and J. T. Weir. 2016. Extinctioin as a driver of avian latitudinal diversity gradients. Evolution 70: 860-872.
}
\author{
Sean A. S. Anderson
}
\examples{
str(all_birds)
}