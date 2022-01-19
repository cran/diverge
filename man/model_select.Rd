\name{model_select}
\alias{model_select}
\title{
Select the best fit model of trait divergence
}
\description{
Compare likelihood and AICc support for up to ten evolutionary models given a dataset.
}
\usage{
model_select(div, ages, me1 = NULL, me2 = NULL, GRAD =NULL, cats=NULL, breakpoint=NULL, 
    domain=NULL, models, starting = NULL, absolute=TRUE, parallel=FALSE, cores=NULL)
}
\arguments{
  \item{div}{
Vector of trait divergences for a set of lineage pairs. Calculated for each pair as abs(trait_val_lineage_2 - trait_val_lineage1). Raw values (i.e. not absolute values) can also be used but must be noted by the user in the argument 'absolute'.  
}
  \item{ages}{
Vector containing the age (i.e. estimated time since divergence) for each pair in the dataset. IMPORTANT: div, ages, GRAD, and breakpoint vectors must be aligned such that div[i] age[i] GRAD[i] and breakpoint[i] represent values for the same lineage pair.
}
  \item{me1}{
Vector containing the measurement error (standard error of mean) for species 1s of each pair (for one species = variance/No.measurements)
}
  \item{me2}{
Vector containing the measurement error (standard error of mean) for species 2s of each pair (for one species = variance/No.measurements)
}
  \item{GRAD}{
Vector containing the gradient position of each pair. This is the value of a continuous variable such as latitude or body size across which parameters are hypothesized to vary. Required for models with 'linear' suffix.
}
  \item{cats}{
Vector containing the category code (0, 1, or 2) for each pair (see package help page for details).
}
  \item{breakpoint}{
Vector of breakpoint times for each pair in the dataset. These are the times AFTER divergence at which a shift occurs in the psi parameter of DA model. Required for DA_bp and DA_bp_linear. See find_mle for details on how to calculate.
}
  \item{domain}{
Vector of length 2 defining the low and high ends of the gradient domain. Essentially identical to the 'xlim' argument in plotting functions. Required for models with 'linear' suffix.
}
  \item{models}{
Character vector naming models for which to compare likelihood estimates and AICc values (options: "BM_null", "BM_linear", "BM_cat", "OU_null", "OU_linear", "OU_linear_sig", "OU_cat", "DA_null", "DA_linear", "DA_cat", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_OU", "DA_OU_linear", "DA_OU_cat", "DA_BM", "DA_BM_linear", "DA_BM_cat", "OU_BM", "OU_BM_linear", "OU_BM_cat"). See find_mle for model descriptions.
}
  \item{starting}{
Optional list in which each element is a matrix of customized parameter from which to launch likelihood searches. The list must be ordered such that the n'th element contains the starting parameter values for the n'th model named in the "models" argument. The nth elements in the list must match the structure required for starting parameters for the nth model in "models" (see the model descriptions in find_mle for details). If starting=NULL, default starting values are used. 
}
  \item{absolute}{
Logical indicating whether 'div' represents absolute value of trait divergence or the raw divergence values.
}
  \item{parallel}{
Logical indicating whether likelihood searches should be conducted in parallel across multiple cores. Not available on windows machines. Defaults to FALSE.
}
  \item{cores}{
If parallel=TRUE, the number of cores on which to run the function. Defaults to all virtual cores.
}
}
\details{
OVERVIEW \cr
model_select estimates parameters and compares the fit of up to 20 evolutionary models to lineage-pair datasets based on likelihood, AIC, and AICc. It is the key function in this package for most users. We strongly recommend that users define the models they want to fit to the data (i.e. only those models that are relevant to distinguishing among hypotheses of interest) rather than simply fitting all models (the default option).

In most cases, if you find that a model fails to fit (e.g. NaN is being returned for likelihood for one or more models), then the default starting parameters are not well-suited to your dataset and you will need to make a list of custom starting parameters (see 'starting' argument above).

RUN TIME CONSIDERATIONS \cr
model_select calls the underlying function find_mle, which launches likelihood searches from a large set of starting parameter values to avoid getting stuck on local optima (see find_mle for more details). This set of starting parameter values grows rapidly with parameter number, so model comparisons that include the most complex evolutionary models (esp DA_wt, DA_wt_linear) can take several minutes. We recommend using the parallel=TRUE option where resources allow when including any model more complex than DA_null. Users can define their own starting parameter sets but should keep in mind this tradeoff between run time and the breadth of parameter space through which to search. 

MODELS OF EVOLUTIONARY TRAIT DIVERGENCE \cr
Users can select up to twenty evolutionary models that differ in how a continuous trait in two lineages "i" and "j" evolves after their initial departure from a common ancestor at T=0. All models require div and age vectors, but some models require additional data. The models are:

(1) BM_null
Description: The trait in both lineages evolves under Brownian Motion (BM) processes. 
Parameter: sig2 (i.e. sigma^2, the dispersion parameter of a BM process; this parameter the same in both lineages).
Optional custom starting matrix: A user-defined starting matrix must be a 1 column matrix containing sig2 starting values.

(2) BM_linear 
Description: Same as BM_null but sig2 varies as a linear function of a continuous variable such as latitude, elevation, or body size. 
Additional Data Required: values of the continuous variable for each pair must be provided in the 'GRAD' vector. 
Parameters: sig2_slope (slope of sig2~gradient relationship), sig2_int (intercept of sig2~gradient relationship).
Optional custom starting matrix: A user-defined starting matrix must be 2 column matrix with col1 = sig2_slope and col2 = sig2_int.

(3) BM_cat 
Description: similar to BM_null but sig2 varies among discrete categories to which different pairs belong, such as allopatric versus sympatric ranges, or different pollinator type.
Additional Data Required: The numeric category code of each pair must be provided in the 'cats' vector. A value must be provided for every pair (i.e. the cats vector must be as long as the data vectors). For example, in a 2-category dataset (e.g. a dataset containing allopatric and sympatric pairs), all pairs of one category can be assigned a 'cats' value of 0, and all pairs of the other category are assigned a 'cats' value of 1. If a third category existed, those pairs would be assigned a 'cats' value of 2, and so on with additional categories. 
Parameters = sig2_1, sig2_2, ... sig2_n, where n is the number of categories.
Optional custom starting matrix: A user-defined p_starting matrix must be an n column matrix where col1 = sig2_1, col2 = sig2_2, ..., coln= sig2_n. 

(4) OU_null 
Description: The trait in both lineages evolves under Ornstein-Uhlenbeck (OU) processes. 
Parameters: alpha (the OU constraint parameter), sig2 - both parameters are shared by the two lineages
Optional custom starting matrix: A user-defined p_starting matrix must be 2 column matrix in the order col1 = alpha and col2 = sig2.

(5) OU_linear 
Description: Same as OU_null but alpha varies as a linear function of a continuous variable such as latitude, elevation, or body size. 
Additional Data Required: Values of the continuous variable for each pair must be provided in the 'GRAD' vector.  
Parameters: alpha_int (intercept of alpha~gradient relationship), alpha_slope (slope of alpha~gradient relationship), sig2.
Optional custom starting matrix: A user-defined p_starting matrix must be 3 column matrix where col1 = alpha_int, col2 = alpha_slope, and col3 = sig2.

(6) OU_linear_sig
Description: Same as OU_linear but sig2, not alpha, varies as a linear function of a continuous variable such as latitude, elevation, or body size. 
Additional Data Required: Values of the continuous variable for each pair must be provided in the 'GRAD' vector.  
Parameters: alpha_int (intercept of alpha~gradient relationship), alpha_slope (slope of alpha~gradient relationship), sig2.
Optional custom starting matrix: A user-defined p_starting matrix must be 3 column matrix where col1 = alpha, col2 = sig2_slope, and col3 = sig2_int.

(7) OU_cat 
Description: similar to BM_cat but alpha (not sig2) varies among discrete categories. 
Additional Data Required: The numeric category code of each pair must be provided in the 'cats' vector. A value must be provided for every pair (i.e. the cats vector must be as long as the data vectors). For example, in a 2-category dataset (e.g. a dataset containing allopatric and sympatric pairs), all pairs of one category can be assigned a 'cats' value of 0, and all pairs of the other category are assigned a 'cats' value of 1. If a third category existed, those pairs would be assigned a 'cats' value of 2, and so on with additional categories. 
Parameters = sig2, alpha_1, alpha_2, ... alpha_n, where n is the number of categories.
Optional custom starting matrix: A user-defined p_starting matrix must be an n+1 column matrix where col1 = sig2, col2 = alpha_1, col3 = alpha_2, ..., col(n+1)= alpha_n. 

(8) DA_null 
Description: Traits evolve under independent OU processes in lineages i and j after they depart from an ancestor, and these processes differ only in the value of their optima (see Anderson and Weir, 2020, 'Inferring speciation drivers from functional trait divergence'. Am Nat, 196, 429-442.).
Parameters: alpha, sig2, and psi (the distance between optima of independent OU processes).
Optional custom starting matrix: A user-defined starting matrix must be 3 column matrix where col1 = alpha, col2 = sig2, and col3 = psi.

(9) DA_linear 
Description: Same as DA_null but psi varies as a linear function of a continuous variable such as latitude, elevation, or body size. 
Additional Data Required: Values of the continuous variable for each pair must be provided in the 'GRAD' vector. 
Parameters: alpha, sig2, psi_slope (slope of psi~gradient relationship), and psi_int (intercept of psi~gradient relationship)
Optional custom starting matrix: A user-defined starting matrix must be 4 column matrix where col1 = alpha, col2 = sig2, col3 = psi_sl, and col4=psi_int

(10) DA_cat 
Description: A DA model in which psi varies among discrete categories. 
Additional Data Required: The numeric category code of each pair must be provided in the 'cats' vector. A value must be provided for every pair (i.e. the cats vector must be as long as the data vectors). For example, in a 2-category dataset (e.g. a dataset containing allopatric and sympatric pairs), all pairs of one category can be assigned a 'cats' value of 0, and all pairs of the other category are assigned a 'cats' value of 1. If a third category existed, those pairs would be assigned a 'cats' value of 2, and so on with additional categories. 
Parameters = alpha, sig2, psi1 (psi for variable category 1), psi2 (psi for variable category 2), ..., psi_n (psi for variable category n).
Optional custom starting matrix: A user-defined p_starting matrix must be a n+2 column matrix where col1 = alpha, col2 = sig2, col3 = psi1, col4 = psi2, ..., col(n+2) = psi_n. 

(11) DA_OU
Description: mixture model in which some proportion of pairs are diverging under a DA process while the rest diverge under an OU process. 
Parameters: alpha, sig2, psi, and 'prop' (the proportion of pairs diverging under a DA process, with 1-prop being the proportion of pairs diverging under the OU process)
Optional custom starting matrix: A user-defined starting matrix must be a 4 column matrix where col1=alpha, col2=sig2, col3=psi, and col4=prop

(12) DA_OU_linear
Description: mixture model in which the proportion of pairs diverging under a DA process (i.e. 'prop' parameter) varies as a linear function of a continuous variable such as latitude, elevation, or body size. 
Additional Data Required: Values of the continuous variable for each pair must be provided in the 'GRAD' vector. 
Parameters: alpha, sig2, psi, prop_slope, and prop_int
Optional custom starting matrix: A user-defined starting matrix must be a 5 column matrix where col1=alpha, col2=sig2, col3=psi, col4=prop_slope, and col5=prop_int

(13) DA_OU_cat
Description: similar to other categorical models but the parameter that varies among categories is 'prop', the proportion of pairs diverging under a DA process (where 1-prop is the proportion diverging under the OU process)
Additional Data Required: Values of the category variable for each pair must be provided in the 'cat' vector. For example, in a 2-category dataset (e.g. a dataset containing allopatric and sympatric pairs), all pairs of one category can be assigned a 'cats' value of 0, and all pairs of the other category are assigned a 'cats' value of 1. If a third category existed, those pairs would be assigned a 'cats' value of 2, and so on with additional categories. 
Parameters = alpha, sig2, psi, prop_1, prop_2, ..., prop_n (for variable category n).
Optional custom starting matrix: A user-defined p_starting matrix must be a n+3 column matrix where col1 = alpha, col2 = sig2, col3 = psi, col4 = prop_1, col5=prop_2, ..., col(n+3) = prop_n. 

(14) OU_BM
Description: mixture model in which some proportion of pairs ('prop' parameter) are diverging under an OU process while the rest diverge under a BM process. 
Parameters: alpha, sig2, and prop (the proportion of pairs diverging under an OU process, with 1-prop being the proportion of pairs diverging under the BM process)
Optional custom starting matrix: A user-defined starting matrix must be a 3 column matrix where col1=alpha, col2=sig2, and col3=prop

(15) OU_BM_linear
Description: mixture model in which the proportion of pairs diverging under an OU process (i.e. 'prop' parameter) varies as a linear function of a continuous variable such as latitude, elevation, or body size. 
Additional Data Required: Values of the continuous variable for each pair must be provided in the 'GRAD' vector. 
Parameters: alpha, sig2, prop_slope, and prop_int
Optional custom starting matrix: A user-defined starting matrix must be a 4 column matrix where col1=alpha, col2=sig2, col3=prop_slope, and col4=prop_int

(16) OU_BM_cat
Description: similar to other categorical models but the parameter that varies among categories is 'prop', the proportion of pairs diverging under an OU process (where 1-prop is the proportion diverging under the OU process)
Additional Data Required: Values of the category variable for each pair must be provided in the 'cat' vector. For example, in a 2-category dataset (e.g. a dataset containing allopatric and sympatric pairs), all pairs of one category can be assigned a 'cats' value of 0, and all pairs of the other category are assigned a 'cats' value of 1. If a third category existed, those pairs would be assigned a 'cats' value of 2, and so on with additional categories. 
Parameters = alpha, sig2, prop_1, prop_2, ..., prop_n (for variable category n).
Optional custom starting matrix: A user-defined p_starting matrix must be a n+2 column matrix where col1 = alpha, col2 = sig2, col3 = prop_1, col4=prop_2, ..., col(n+2) = prop_n. 

(17) DA_wt
Description: Like DA_null but with a discrete shift in psi that occurs after a wait time shared by all pairs in the dataset. DA_wt estimates this wait time based on divergence levels of pairs of different ages, so datasets should contain many pairs from a broad range of ages. 
Parameters: alpha, sig2, psi1 (distance between optima prior to the wait time), psi2 (distance between optima after the wait time), and wt (the wait time to a shift in psi; the time after initial divergence from a common ancestor at which a shift occurs).
Optional custom starting matrix: A user-defined p_starting matrix must be 5 column matrix where col1 = alpha, col2 = sig2, col3 = psi1, col4=psi2, and col5=wt.

(18) DA_bp 
Description: Like DA_null but with a discrete shift in psi that occurs after some "breakpoint" time (bp) known to the user. Differs from DA_wt in two ways: (1) the timing of the discrete shift isn't isn't estimated as a parameter; it is instead provided by the user in the breakpoint vector that contains a bp value for each pair, and (2) the timing of the discrete shift is not shared by all pairs.
Additional Data Required: Values of bp for each pair must be provided in the breakpoint vector. 
Parameters: alpha, sig2, psi1, psi2
Optional custom starting matrix: A user-defined p_starting matrix must be 4 column matrix where col1 = alpha, col2 = sig2, col3 = psi1, and col4=psi2
Additional Notes: This model is useful for testing hypotheses related to dated biogeographic events. For example, consider a set of lineage pairs with diverge times between 3 and 7my. If a river formed 1mya and isolated some of these already-diverging pairs on either side, we might hypothesize that this caused a shift in psi for those pairs, and a breakpoint model can test this. Bp is then calculated as bp = AgeOfPair - 1my for each pair divided by the river. For pairs that aren't divided by the river, we set bp = 0. 
IMPORTANT: datasets must contain pairs that HAVE and HAVE NOT experienced the hypothesized shift in psi. In the river example, the dataset must contain pairs whose lineages were NOT separated by the river. BP IS SET TO 0 FOR ALL SUCH PAIRS. Parameter estimation and model_selection is most accurate when ~60-75 percent of the dataset is comprised of bp=0 (aka "single-epoch") pairs. 

(19) DA_wt_linear 
Description: Like DA_wt except psi varies with a continuous gradient both before and after the wait time. 
Additional Data Required: Values of the continuous variable for each pair must be provided in the 'GRAD' vector. 
Parameters: alpha, sig2, psi1_sl (slope of psi1~gradient relationship), psi1_int (intercept of psi1~gradient relationship), psi2_sl (slope of psi2~gradient relationship), psi2_int (intercept of psi2~gradient relationship), wt
Optional custom starting matrix: A user-defined p_starting matrix must be 7 column matrix where col1 = alpha, col2 = sig2, col3 = psi1_slope, col4=psi1_int, col5=psi2_slope, col6=psi1_int, and col7=wt

(20) DA_bp_linear 
Description: Like DA_bp except psi varies with a continuous gradient both before and after the wait time. 
Additional Data Required: Values of the continuous variable for each pair must be provided in the 'GRAD' vector AND values of bp for each pair must be provided in the breakpoint vector. 
Parameters: alpha, sig2, psi1_slope (distance between optima prior to bp), psi1_int, psi2_slope (slope of distance beteen optima after bp), psi2_int (intercept of distance between optima after bp)
Optional custom starting matrix: A user-defined p_starting matrix must be 6 column matrix where col1 = alpha, col2 = sig2, col3 = psi1_slope, col4=psi1_int, col5=psi2_slope, and col6=psi2_int

In general, BM-models may be interpreted as representing random drift or fluctuating selection in trait evolution, OU-models may be interpreted as representing parallel adaptation or shared stabilizing selection (these are mathematically identical in this context), and DA-models may be interpreted as representing divergent adaptation (i.e. evolution in which lineages are generally pulled by selection toward alternative adaptive optima following initial divergence from a common ancestor)

}
\value{
Returns a summary matrix of likelihood results, AIC calculations, and maximum likelihood parameter estimates for each chosen model. 
}
\author{
Sean A.S. Anderson and Jason T. Weir
}
\examples{
# to fit model, need a dataset
# for illustrative purposes, can simulate dataset under DA_linear 
# assume a 0-1000m elevational gradient.
ages = rep(c(0.5, 1, 1.5, 2, 3, 8), 25)
GRAD_cats = rep(c(0, 250, 500, 750, 1000), 30)
GRAD=c(rep(GRAD_cats[1], 30), rep(GRAD_cats[2],30), rep(GRAD_cats[3],30), 
  rep(GRAD_cats[4],30), rep(GRAD_cats[5],30))
alpha = 0.8
sig2 = 0.2
psi_sl = -0.01
psi_int = 2
sis_div = simulate_div(model="DA_linear", ages=ages, GRAD=GRAD, 
  pars=c(alpha, sig2, psi_sl, psi_int))


# run model comparison with DA_null and DA_linear
\donttest{
res <- model_select(div=sis_div, ages=ages, GRAD=GRAD, domain=c(0, 1000), 
  models=c("DA_null", "DA_linear"))
  }
}
