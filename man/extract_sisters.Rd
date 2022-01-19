\name{extract_sisters}
\alias{extract_sisters}
\title{
Extract sister pairs
}
\description{
Find the names and divergence times for all sister pairs in phylogeny
}
\usage{
extract_sisters(tree, sis_age=FALSE, mol_clock=NULL, crown_age=NULL)
}
\arguments{
  \item{tree}{
A phylo object. Branch lengths must represent time or genetic distances with known molecular clock. 
}
  \item{sis_age}{
Logical indicating whether or not sister pair ages are to be estimated. 
}
  \item{mol_clock}{
The estimated constant rate of sequence evolution per million years.
}
  \item{crown_age}{
The estimated crown age in millions of years.
}
}
\details{
extract_sisters takes a phylo object and identifies all sister pairs by name. If sis_age=TRUE, then extract_sisters also estimates the age of each pair (in Ma) based on branch lengths and either 1) a given (constant) molecular clock, or 2) an estimate for the crown age. If no molecular clock rate or crown age is supplied by the user, the branch lengths are assumed to represent time in Ma. CAUTION: age estimates will be returned for any tree that contains branch lengths, so it's up to the user to ensure that the tree is ultrametric.
}
\value{
Returns a data frame in which each row corresponds to a sister pair and contains the taxon names for the two lineages, and if sis_age=TRUE, then the estimated divergence time (i.e. age) for that pair.
}
\author{
Sean A.S. Anderson
}
\examples{
# Note: this function automatically generates a warning to remind users that it is
# their responsibility to ensure the tree is ultrametric. This warning can be suppressed
# with suppressWarnings

# extract avian sister pairs and calculate ages based on branch lengths in Ma
bird_sis = extract_sisters(tree=all_birds, sis_age=TRUE)

\donttest{
# extract avian sister pairs and calculate ages based on molecular clock of 1.1
# note: the actual branch lengths are in Ma, this is just a toy example
bird_sis = extract_sisters(tree=all_birds, sis_age=TRUE, mol_clock = 1.1)

# extract avian sister pairs and calculate ages based on a crown age of 30 Ma
# note: the actual branch lengths are in Ma, this is just a toy example
bird_sis = extract_sisters(tree=all_birds, sis_age=TRUE, crown_age=30)
}
}


