simulate_div <-
function(model, pars, ages, me_prop=NULL, GRAD=NULL, cats=NULL, breakpoint=NULL, Nsets=1) {

  if(model %in% c("BM_null", "OU_null", "BM_linear", "OU_linear", "OU_linear_sig", 
    "DA_null", "DA_linear", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }
  
  if(is.numeric(me_prop) & length(me_prop)!=1 & length(me_prop)!= length(ages)) {
    stop("me_prop should be a single value or a vector of the same length as 'ages'")
  }
  # calculate vectors for variance and expected value of one lineage in each pair given 
  # (a) the model 
  # (b) model paramters, and 
  # (c) the age of the pair, plus gradient and breakpoint vectors where approrpriate
  # For DA models:
  # theta of one lineage in each pair is arbitrarily set to zero
  # theta for the other is set to the psi value 
  if(model == "BM_null") {
    B = pars[1]
    u2 = rep(0, length(ages))
    var = 2*B*ages
  }
  if(model == "BM_linear") {
    B_slope = pars[1]
    B_int = pars[2]
    B = GRAD*B_slope + B_int
    u2 = rep(0, length(ages))
    var = 2*B*ages
  }
  if(model == "OU_linear") {
    A_slope = pars[2]
    A_int = pars[1]
    B = pars[3]
    A = GRAD*A_slope + A_int
    u2 = rep(0, length(ages))
    var = (B/A)*(1-exp(-2*A*ages))
  }
  if(model != "BM_null" | model != "BM_linear" | model != "OU_linear") {
    A=pars[1]
    B=pars[2] 
  }
  if(model == "OU_null") {
    u2 = rep(0, length(ages))
    var = (B/A)*(1-exp(-2*A*ages))
  }
  if(model == "DA_null") {
    if(length(pars) != 3) {
      stop("Wait! You've chosen a model that requires 3 parameters")
    }
    peakdist = pars[3]
    var = (B/A)*(1-exp(-2*A*ages)) # vector of variances
    u2 = peakdist*(1-exp(-A*ages)) # vector of expected values for sister 2
  }
  if(model == "DA_linear") {
    if(length(pars) != 4) {
      stop("You've chosen a model that requires 4 parameters")
    }
    if(is.numeric(GRAD) == FALSE) {
      stop("You've chosen a model that requires a GRAD argument containing the gradient positions for each sister pair")
    }
    grad_slope=pars[3]
    grad_intercept=pars[4]
    psi = GRAD*grad_slope + grad_intercept # vector of psi values based on the position of each pair on the gradient
    var = (B/A)*(1-exp(-2*A*ages)) # vector of variances for sisters 1 and 2
    u2 = psi*(1-exp(-A*ages)) # vector of expected values for sister 2
  }
  if(model == "DA_cat") {
    if(length(pars) != 4 & length(pars) != 5) {
      stop("You've chosen a model that requires 4 or 5 parameters")
    }
    psi1 = pars[3]
    psi2 = pars[4]
    if(length(pars) == 5) psi3 = pars[5]
    var = (B/A)*(1-exp(-2*A*(ages)))
    u2 = rep(0, length(ages))
    for(i in 1:length(ages)) {
      if(cats[i]==0) u2[i] = psi1*(1-exp(-A*ages[i]))
      if(cats[i]==1) u2[i] = psi2*(1-exp(-A*ages[i]))
      if(cats[i]==2) u2[i] = psi3*(1-exp(-A*ages[i]))
    }
  }
  if(model == "DA_wt") {
    if(length(pars) !=5) {
      stop("You've chosen a model that requires 5 parameters")
    }
    psi1=pars[3]
    psi2=pars[4]
    wait_time=pars[5]
    var = (B/A)*(1-exp(-2*A*(ages)))
    u2 = rep(0, length(ages))
    for(i in 1:length(ages)) {
      if(ages[i] > wait_time) {
        u2[i] = exp(-A*ages[i])*(psi1*(exp(A*wait_time)-1)+psi2*(exp(A*ages[i])-exp(A*wait_time)))
      } else {
        u2[i] = psi1*(1-exp(-A*ages[i])) 
      }
    }
  }
  if(model == "DA_bp") {
    if(length(pars) !=4) {
      stop("You've chosen a model that requires 4 parameters")
    }
    if(is.null(breakpoint)){
      stop("You've chosen the breakpoint model and haven't provided a breakpoint vector")
    }
    psi1=pars[3]
    psi2=pars[4]
    var = (B/A)*(1-exp(-2*A*(ages)))
    u2 = rep(0, length(ages))
    for(i in 1:length(breakpoint)) {
      if(breakpoint[i] > 0) {
        u2[i] = exp(-A*ages[i])*(psi1*(exp(A*breakpoint[i])-1)+psi2*(exp(A*ages[i])-exp(A*breakpoint[i])))
      } else {
        u2[i] = psi1*(1-exp(-A*ages[i])) 
      }
    }
  }
  if(model == "DA_wt_linear") {
    if(length(pars) != 7) {
      stop("Hold it! You've chosen a model that requires 7 parameters")
    }
    if(length(GRAD) == 0) {
      stop("Wait! You've chosen a model that requires a GRAD argument containing the gradient positions for each sister pair")
    }
    psi1_slope=pars[3]
    psi1_int=pars[4]
    psi2_slope=pars[5]
    psi2_int=pars[6]
    wt=pars[7]
    var = (B/A)*(1-exp(-2*A*ages))
    u2 = rep(0, length(ages))
    for(i in 1:length(ages)) {
      if(ages[i] > wt) {
        psi1 = GRAD[i]*psi1_slope + psi1_int
        psi2 = GRAD[i]*psi2_slope + psi2_int
        u2[i] = exp(-A*ages[i])*(psi1*(exp(A*wt)-1)+psi2*(exp(A*ages[i])-exp(A*wt)))
      } else {
        psi = GRAD[i]*psi1_slope + psi1_int
        u2[i] = psi*(1-exp(-A*ages[i])) 
      }
    }
  }
  if(model == "DA_bp_linear") {
    if(length(pars) != 6) stop("Hold it! You've chosen a model that requires 6 params")
    if(is.null(GRAD)) stop("Wait! You've chosen a model that requires a GRAD argument")
    if(is.null(breakpoint)) stop("Wait! You've chosen to simulate a model that requires a breakpoint vector")
    psi1_slope = pars[3]
    psi1_int = pars[4]
    psi2_slope = pars[5]
    psi2_int = pars[6]
    var = (B/A)*(1-exp(-2*A*ages))
    u2 = rep(0, length(ages))
    for(i in 1:length(ages)) {
      if(breakpoint[i] > 0) {
        psi1 = GRAD[i]*psi1_slope + psi1_int
        psi2 = GRAD[i]*psi2_slope + psi2_int
        u2[i] = exp(-A*ages[i])*(psi1*(exp(A*breakpoint[i])-1)+psi2*(exp(A*ages[i])-exp(A*breakpoint[i])))
      } else {
        psi = GRAD[i]*psi1_slope+ psi1_int
        u2[i] = psi*(1-exp(-A*ages[i])) 
      }
    }
  }
  # if measurement error is to be added, calculate it as the user-defined
  # proportion of the variance in each lineage
  if(is.null(me_prop)==FALSE) var = var + 2*(me_prop*var/2)^2
  # use variance and expected value vectors to simulate trait divergence
  if(Nsets == 1) EDsisters = truncnorm::rtruncnorm(n=length(ages), a=0, mean=u2, sd=sqrt(var))
  if(Nsets > 1) EDsisters = as.list(as.data.frame(replicate(n=Nsets, 
    expr=truncnorm::rtruncnorm(n=length(ages), a=0, mean=u2, sd=sqrt(var)))))
  return(EDsisters)
}