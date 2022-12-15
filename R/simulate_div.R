simulate_div <-
function(model, pars, ages, me_prop=NULL, GRAD=NULL, cats=NULL, breakpoint=NULL, Nsets=1) {

  if(model %in% c("BM_null", "OU_null", "DA_null", "BM_linear", "OU_linear", "DA_linear", 
    "DA_cat", "DA_wt", "DA_bp", "DA_wt_linear", "DA_bp_linear", "DA_OU", "DA_OU_linear", 
    "DA_OU_cat", "DA_BM", "DA_BM_linear", "DA_BM_cat", "OU_BM", "OU_BM_linear", "OU_BM_cat") == FALSE) {
    stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
  }
  
  if(is.numeric(me_prop) & length(me_prop)!=1 & length(me_prop)!= length(ages)) {
    stop("me_prop should be a single value or alpha vector of the same length as 'ages'")
  }
  if(model == "BM_null") {
    sig2 = pars[1]
    u2 = rep(0, length(ages))
    var = 2*sig2*ages
  }
  if(model == "BM_linear") {
    sig2_slope = pars[1]
    sig2_int = pars[2]
    sig2 = GRAD*sig2_slope + sig2_int
    u2 = rep(0, length(ages))
    var = 2*sig2*ages
  }
  if(model == "OU_linear") {
    A_slope = pars[2]
    A_int = pars[1]
    sig2 = pars[3]
    alpha = GRAD*A_slope + A_int
    u2 = rep(0, length(ages))
    var = (sig2/alpha)*(1-exp(-2*alpha*ages))
  }
  if(model != "BM_null" | model != "BM_linear" | model != "OU_linear") {
    alpha=pars[1]
    sig2=pars[2] 
  }
  if(model == "OU_null") {
    u2 = rep(0, length(ages))
    var = (sig2/alpha)*(1-exp(-2*alpha*ages))
  }
  if(model == "DA_null") {
    if(length(pars) != 3) {
      stop("You need to provide 3 parameters to simulate under the chosen model")
    }
    psi = pars[3]
    var = (sig2/alpha)*(1-exp(-2*alpha*ages)) # vector of variances
    u2 = psi*(1-exp(-alpha*ages)) # vector of expected values for sister 2
  }
  if(model == "DA_linear") {
    if(length(pars) != 4) {
      stop("You need to provide 4 parameters to simulate under the chosen model")
    }
    if(is.numeric(GRAD) == FALSE) {
      stop("You've chosen a model that requires a GRAD argument containing the gradient positions for each sister pair")
    }
    grad_slope=pars[3]
    grad_intercept=pars[4]
    psi = GRAD*grad_slope + grad_intercept # vector of psi values based on the position of each pair on the gradient
    var = (sig2/alpha)*(1-exp(-2*alpha*ages)) # vector of variances for sisters 1 and 2
    u2 = psi*(1-exp(-alpha*ages)) # vector of expected values for sister 2
  }
  if(model == "DA_cat") {
    if(length(pars) != 4 & length(pars) != 5) {
      stop("You need to provide 4 or 5 parameters to simulate under the chosen model")
    }
    psi1 = pars[3]
    psi2 = pars[4]
    if(length(pars) == 5) psi3 = pars[5]
    var = (sig2/alpha)*(1-exp(-2*alpha*(ages)))
    u2 = rep(0, length(ages))
    for(i in 1:length(ages)) {
      if(cats[i]==0) u2[i] = psi1*(1-exp(-alpha*ages[i]))
      if(cats[i]==1) u2[i] = psi2*(1-exp(-alpha*ages[i]))
      if(cats[i]==2) u2[i] = psi3*(1-exp(-alpha*ages[i]))
    }
  }
  if(model == "DA_OU_linear") {
    if(length(pars) != 5) {
      stop("You need to provide 5 parameters to simulate under the chosen model")
    }
    props = pars[4]*unique(GRAD) + pars[5]
    if(all(props<=1 & props>=0)==FALSE) {
      stop("Proportions are not between 0 and 1; You have to pick different slope and intercept parameters to go with the gradient provided")
    }
    if(is.numeric(GRAD) == FALSE) {
      stop("You've chosen a model that requires a GRAD argument containing the gradient positions for each sister pair")
    }
    # note! a requirement for simulating data under this model is that there are an equal number of pairs at each gradient value
    # i.e. it won't work if every pair has its own value for the gradient
    k = length(GRAD)/length(unique(GRAD))
    da_inds = list()
    for(i in 1:length(unique(GRAD))) {
      rng = (i*k - (k-1)):(i*k)
      da_inds[[i]] = rng[1:round(props[i]*k)]
    }
    da_inds = do.call(c,da_inds)
    u2=rep(0, length(ages))
    var = (sig2/alpha)*(1-exp(-2*alpha*ages))
    u2[da_inds] = pars[3]*(1-exp(-alpha*ages[da_inds]))
  }
  if(model == "DA_BM_linear") {
    if(length(pars) != 5) {
      stop("You need to provide 5 parameters to simulate under the chosen model")
    }
    props = pars[4]*unique(GRAD) + pars[5]
    if(all(props<=1 & props>=0)==FALSE) {
      stop("Proportions are not between 0 and 1; You have to pick different slope and intercept parameters to go with the gradient provided")
    }
    # note! a requirement for simulating data under this model is that there are an equal number of pairs at each gradient value
    # i.e. it won't work if every pair has its own value for the gradient
    k = length(GRAD)/length(unique(GRAD))
    da_inds = list()
    for(i in 1:length(unique(GRAD))) {
      rng = (i*k - (k-1)):(i*k)
      da_inds[[i]] = rng[1:round(props[i]*k)]
    }
    da_inds = do.call(c,da_inds)
    bm_inds = which(1:length(ages) %in% da_inds == FALSE)
    var=u2=rep(0, length(ages))
    var[da_inds] = (sig2/alpha)*(1-exp(-2*alpha*ages[da_inds]))
    var[bm_inds] = 2*sig2*ages[bm_inds]
    u2[da_inds] = pars[3]*(1-exp(-alpha*ages[da_inds]))
  }
  if(model == "OU_BM_linear") {
    if(length(pars) != 4) {
      stop("You need to provide 4 parameters to simulate under the chosen model")
    }
    props = pars[3]*unique(GRAD) + pars[4]
    if(all(props<=1 & props>=0)==FALSE) {
      stop("Proportions are not between 0 and 1; You have to pick different slope and intercept parameters to go with the gradient provided")
    }
    # note! a requirement for simulating data under this model is that there are an equal number of pairs at each gradient value
    # i.e. it won't work if every pair has its own value for the gradient
    k = length(GRAD)/length(unique(GRAD))
    ou_inds = list()
    for(i in 1:length(unique(GRAD))) {
      rng = (i*k - (k-1)):(i*k)
      ou_inds[[i]] = rng[1:round(props[i]*k)]
    }
    ou_inds = do.call(c,ou_inds)
    bm_inds = which(1:length(ages) %in% ou_inds == FALSE)
    var=u2=rep(0, length(ages))
    var[ou_inds] = (sig2/alpha)*(1-exp(-2*alpha*ages[ou_inds]))
    var[bm_inds] = 2*sig2*ages[bm_inds]
  }
  if(model == "DA_OU") {
    if(length(pars) != 4) {
      stop("You need to provide 4 parameters to simulate under the chosen model")
    }
    if(abs(pars[4]*length(ages)) - round(pars[4]*length(ages)) != 0) {
      stop("Please provide a proportion parameter that results in the division of pairs into whole number subsets")
    }
    da_inds = 1:round(pars[4]*length(ages))
    var = (sig2/alpha)*(1-exp(-2*alpha*ages))
    u2=rep(0, length(ages))
    u2[da_inds] = pars[3]*(1-exp(-alpha*ages[da_inds]))
  }
  if(model == "DA_BM") {
    if(length(pars) != 4) {
      stop("You need to provide 4 parameters to simulate under the chosen model")
    }
    if(abs(pars[4]*length(ages)) - round(pars[4]*length(ages)) != 0) {
      stop("Please provide a proportion parameter that results in the division of pairs into whole number subsets")
    }
    da_inds = 1:round(pars[4]*length(ages))
    bm_inds = which(1:length(ages) %in% da_inds == FALSE)
    var=u2=rep(0, length(ages))
    var[da_inds] = (sig2/alpha)*(1-exp(-2*alpha*ages[da_inds]))
    var[bm_inds] = 2*sig2*ages[bm_inds]
    u2[da_inds] = pars[3]*(1-exp(-alpha*ages[da_inds]))
  }
  if(model == "OU_BM") {
    if(length(pars) != 3) {
      stop("You need to provide 3 parameters to simulate under the chosen model")
    }
    if(abs(pars[3]*length(ages)) - round(pars[3]*length(ages)) != 0) {
      stop("Please provide a proportion parameter that results in the division of pairs into whole number subsets")
    }
    ou_inds = 1:round(pars[3]*length(ages))
    bm_inds = which(1:length(ages) %in% ou_inds == FALSE)
    var=u2=rep(0, length(ages))
    var[ou_inds] = (sig2/alpha)*(1-exp(-2*alpha*ages[ou_inds]))
    var[bm_inds] = 2*sig2*ages[bm_inds]
  }
  if(model == "DA_OU_cat") {
    if(is.numeric(cats) == FALSE) {
      stop("You've chosen a model that requires a 'cats' argument containing the category code for each pair")
    }
    # note! a requirement for simulating data under this model is that there are an equal number of pairs for each category value
    k = length(cats)/length(unique(cats))
    da_inds = list()
    for(i in 1:length(unique(cats))) {
      rng = (i*k - (k-1)):(i*k)
      da_inds[[i]] = rng[1:round(pars[3+i]*k)]
    }
    da_inds = do.call(c,da_inds)
    u2=rep(0, length(ages))
    var = (sig2/alpha)*(1-exp(-2*alpha*ages))
    u2[da_inds] = pars[3]*(1-exp(-alpha*ages[da_inds]))
  }
  if(model == "DA_BM_cat") {
    if(is.numeric(cats) == FALSE) {
      stop("You've chosen a model that requires a 'cats' argument containing the gradient positions for each sister pair")
    }
    # note! a requirement for simulating data under this model is that there are an equal number of pairs for each category value
    k = length(cats)/length(unique(cats))
    da_inds = list()
    for(i in 1:length(unique(cats))) {
      rng = (i*k - (k-1)):(i*k)
      da_inds[[i]] = rng[1:round(pars[3+i]*k)]
    }
    da_inds = do.call(c,da_inds)
    bm_inds = which(1:length(ages) %in% da_inds == FALSE)
    var=u2=rep(0, length(ages))
    var[da_inds] = (sig2/alpha)*(1-exp(-2*alpha*ages[da_inds]))
    var[bm_inds] = 2*sig2*ages[bm_inds]
    u2[da_inds] = pars[3]*(1-exp(-alpha*ages[da_inds]))
  }
  if(model == "OU_BM_cat") {
    if(is.numeric(cats) == FALSE) {
      stop("You've chosen a model that requires a 'cats' argument containing the gradient positions for each sister pair")
    }
    # note! a requirement for simulating data under this model is that there are an equal number of pairs for each category value
    k = length(cats)/length(unique(cats))
    ou_inds = list()
    for(i in 1:length(unique(cats))) {
      rng = (i*k - (k-1)):(i*k)
      ou_inds[[i]] = rng[1:round(pars[2+i]*k)]
    }
    ou_inds = do.call(c,ou_inds)
    bm_inds = which(1:length(ages) %in% ou_inds == FALSE)
    var=u2=rep(0, length(ages))
    var[ou_inds] = (sig2/alpha)*(1-exp(-2*alpha*ages[ou_inds]))
    var[bm_inds] = 2*sig2*ages[bm_inds]
  }
  if(model == "DA_wt") {
    if(length(pars) !=5) {
      stop("You need to provide 5 parameters to simulate under the chosen model")
    }
    psi1=pars[3]
    psi2=pars[4]
    wait_time=pars[5]
    var = (sig2/alpha)*(1-exp(-2*alpha*(ages)))
    u2 = rep(0, length(ages))
    for(i in 1:length(ages)) {
      if(ages[i] > wait_time) {
        u2[i] = exp(-alpha*ages[i])*(psi1*(exp(alpha*wait_time)-1)+psi2*(exp(alpha*ages[i])-exp(alpha*wait_time)))
      } else {
        u2[i] = psi1*(1-exp(-alpha*ages[i])) 
      }
    }
  }
  if(model == "DA_bp") {
    if(length(pars) !=4) {
      stop("You need to provide 4 parameters to simulate under the chosen model")
    }
    if(is.null(breakpoint)){
      stop("You've chosen the breakpoint model and haven't provided a breakpoint vector")
    }
    psi1=pars[3]
    psi2=pars[4]
    var = (sig2/alpha)*(1-exp(-2*alpha*(ages)))
    u2 = rep(0, length(ages))
    for(i in 1:length(breakpoint)) {
      if(breakpoint[i] > 0) {
        u2[i] = exp(-alpha*ages[i])*(psi1*(exp(alpha*breakpoint[i])-1)+psi2*(exp(alpha*ages[i])-exp(alpha*breakpoint[i])))
      } else {
        u2[i] = psi1*(1-exp(-alpha*ages[i])) 
      }
    }
  }
  if(model == "DA_wt_linear") {
    if(length(pars) != 7) {
      stop("You need to provide 7 parameters to simulate under the chosen model")
    }
    if(length(GRAD) == 0) {
      stop("Wait! You've chosen a model that requires a GRAD argument containing the gradient positions for each sister pair")
    }
    psi1_slope=pars[3]
    psi1_int=pars[4]
    psi2_slope=pars[5]
    psi2_int=pars[6]
    wt=pars[7]
    var = (sig2/alpha)*(1-exp(-2*alpha*ages))
    u2 = rep(0, length(ages))
    for(i in 1:length(ages)) {
      if(ages[i] > wt) {
        psi1 = GRAD[i]*psi1_slope + psi1_int
        psi2 = GRAD[i]*psi2_slope + psi2_int
        u2[i] = exp(-alpha*ages[i])*(psi1*(exp(alpha*wt)-1)+psi2*(exp(alpha*ages[i])-exp(alpha*wt)))
      } else {
        psi = GRAD[i]*psi1_slope + psi1_int
        u2[i] = psi*(1-exp(-alpha*ages[i])) 
      }
    }
  }
  if(model == "DA_bp_linear") {
    if(length(pars) != 6) stop("You need to provide 6 parameters to simulate under the chosen model")
    if(is.null(GRAD)) stop("Wait! You've chosen a model that requires a GRAD argument")
    if(is.null(breakpoint)) stop("Wait! You've chosen to simulate under a model that requires a breakpoint vector")
    psi1_slope = pars[3]
    psi1_int = pars[4]
    psi2_slope = pars[5]
    psi2_int = pars[6]
    var = (sig2/alpha)*(1-exp(-2*alpha*ages))
    u2 = rep(0, length(ages))
    for(i in 1:length(ages)) {
      if(breakpoint[i] > 0) {
        psi1 = GRAD[i]*psi1_slope + psi1_int
        psi2 = GRAD[i]*psi2_slope + psi2_int
        u2[i] = exp(-alpha*ages[i])*(psi1*(exp(alpha*breakpoint[i])-1)+psi2*(exp(alpha*ages[i])-exp(alpha*breakpoint[i])))
      } else {
        psi = GRAD[i]*psi1_slope+ psi1_int
        u2[i] = psi*(1-exp(-alpha*ages[i])) 
      }
    }
  }
  if(is.null(me_prop)==FALSE) var = var + 2*(me_prop*var/2)^2
  if(Nsets == 1) EDsisters = abs(rnorm(n=length(ages),mean=u2, sd=sqrt(var)))
  if(Nsets > 1) EDsisters = as.list(as.data.frame(replicate(n=Nsets, expr=abs(rnorm(n=length(ages), mean=u2, sd=sqrt(var))))))
  return(EDsisters)
}

