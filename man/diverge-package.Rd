\name{diverge-package}
\alias{diverge-package}
\alias{diverge}
\docType{package}
\title{
Evolutionary Trait Divergence Between Sister Species and Other Paired Lineages
}
\description{
This package can be used to study the dynamics of continuous trait divergence in sister species and other paired lineages. The primary aim is to characterize the evolutionary processes responsible for continuous trait differentiation between lineages following their departure from a common ancestor. Alternative mechanisms of trait divergence such as divergent selection, drift, and stabilizing selection are expected to leave unique signatures in the distribution of trait differentiation in datasets comprised of many species pairs. The package uses maximum likelihood to estimate parameters and compare the fits of alternative process-based models to such distributions. Key functions for most users will be 'model_select' and 'find_mle'. See details for a description of each evolutionary model, data requirements, and other considerations. 
}
\details{
\tabular{ll}{
Package: \tab diverge\cr
Type: \tab Package\cr
Version: \tab 1.0.2\cr
Date: \tab 2019-26-07\cr
License: \tab GPL (>= 2) \cr
}

BASIC MODELS: \cr
Alternative mechanisms of trait divergence are approximated here by stochastic process models. We use Brownian Motion (BM) models to approximate divergence under random drift or fluctuating selection; Ornstein-Uhlenbeck (OU) models in which two diverging lineages share an optima to approximate divergence under parallel and/or shared stabilizing selection; and OU models in which two diverging lineages are pulled toward alternative optima to approximate divergent adaptation (this last suite of models we refer to as "DA" models). In the DA models, the difference between optima of paired lineages (psi) measures the extent of divergent selection. A more detailed explanation of the 11 model variants used in this package is provided in the model_select man page.

MODEL EXTENSIONS: \cr
(1) Effect of continuous variables on model parameters. \cr
Rates of evolution and patterns of divergent selection in lineage pairs can vary across continuous biological or ecological gradients such as latitude, elevation, and body mass. We follow the approach used in the EvoRAG package (Weir 2014), whereby model parameters are represented as linear functions of the continuous gradient. Tests of these extended models require as input the gradient value for each pair. \cr
(2) Effect of categorical variables on model parameters. \cr
Just as patterns of divergent selection can vary across continuous ecological gradients, they can also vary across categorical variables such as relative geographic range (i.e. allopatric versus sympatric), and different pollinator types. We model this variation by allowing psi to vary among the user-defined categories. Tests of the categorical model (DA-cat) require as input a vector containing the category of each pair (coded as 0, 1, 2, for a maximum of 3 categories). \cr
(3) Discrete temporal change in extent of divergent selection. \cr
In models of divergent selection (DA models), lineages in a pair are pulled toward alternative adaptive optima following their initial divergence from an ancestor. Often one or both of these lineages will experience a secondary shift in their adaptive environment which changes the position of their optima. This manifests in the DA models as a discrete shift in psi, referred to as an "epoch shift". WAIT-TIME models (DA_wt, DA_wt_linear) estimate the timing of this shift as a wait time shared by all pairs in the dataset. BREAKPOINT models (DA_bp, DA_bp_ilnear) require users to provide a set of breakpoint values, which are the times after divergence at which a shift in psi is hypothesized to occur. These will be unique to pairs of different ages. Pairs in which no shift is expected to have occurred are assigned breakpoint value of zero. See the model_select description of DA_bp for more details. 

DATA REQUIREMENTS: \cr
diverge is designed to work with datasets of lineage-pair contrasts in a continuous trait. At minimum, datasets must contain two components for each lineage pair: (1) a measure of the difference between the two lineages in a trait, and (2) an estimate of time since the two lineages initially diverged. Power and model performance increase with the number of lineage pairs in a dataset. Several functions in this package use the same arguments for required data. These include:

div - a vector of trait differences for each pair in the dataset. Required for all evolutionary models. These are calculated for each pair as abs(trait_val_lineage_2 - trait_val_lineage1). Raw values (i.e. not absolute values) can also be used but this must be noted by the user in the absolute argument of most functions.

ages - a vector containing the age (i.e. estimated time since divergence) for each pair in the dataset.

GRAD - a vector containing the gradient position of each pair. This is the value of a continuous variable such as latitude or body size across which parameters are hypothesized to vary. Required for all models with the 'linear' suffix.

cats - a vector containing the category code (0, 1, or 2) for each pair in the dataset. Required when using the DA-cat model only. 

breakpoint - a vector of breakpoint times for each pair in the dataset. These are the times after divergence at which a shift occurs in the psi parameter of a DA model. Required for DA_bp and DA_bp_linear. See model_select for details on how to calculate.

IMPORTANT: \cr
div, ages, GRAD, cats, and breakpoint vectors must be aligned such that div[i], age[i], grad[i], cat[i], and breakpoint[i] represent values for the same i'th pair.

PARALLELIZATION: \cr
Many functions in 'diverge' can be run in parallel across multiple cores through the 'parallel' logical argument. This functionality uses forking and is therefore unavailable on windows machines. All examples shown in function manuals are run in serial. We encourage users to insert parallel=TRUE in these commands where resources allow.  

CUSTOMIZATION: \cr
The code for several key functions has been lightly commented to aid users wishing to hack functions for custom use.
}
\author{
Sean A. S. Anderson and Jason T. Weir

Maintainer: Sean A.S. Anderson <seanas.anderson@mail.utoronto.ca>
}
