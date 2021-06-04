likelihood_func <-
function (parameters, model, div, ages, me1 = NULL, me2 = NULL, cats=NULL, GRAD=NULL, bp=NULL, absolute=TRUE) {
    if (model == "BM_null") {
      sig2=parameters[1]
      var=2*sig2*ages
      if(is.null(me1) == FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE) {
        dens = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute == FALSE) {
        dens = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "BM_linear") {
      sig2_slope=parameters[1]
      sig2_int=parameters[2]
      sig2=GRAD*sig2_slope + sig2_int
      var=2*sig2*ages
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE) {
        dens = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute == FALSE) {
        dens = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "BM_cat") {
      # parameters are in the order sig1, sig2, sig3, etc.
      if(is.null(cats)) {
        stop("Hold up! You haven't included category values in the 'cats' vector")
      }
      if(length(parameters) != length(unique(cats))) {
        stop("The number of paramters provided doesn't match the number of categories")
      }
      sigs=rep(NA, length(ages))
      for(i in 1:length(unique(cats))) sigs[cats==unique(cats)[i]] = parameters[i]
      var = 2*sigs*ages
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE) {
        dens = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute == FALSE) {
        dens = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "OU_null") {
      A=parameters[1]
      sig2=parameters[2]
      var = sig2*(1-exp(-2*A*ages))/A
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE){
        dens = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute==FALSE) {
        dens = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "OU_linear") {
      A_int=parameters[1]
      A_slope=parameters[2]
      sig2=parameters[3]
      A=GRAD*A_slope + A_int
      var = sig2*(1-exp(-2*A*ages))/A
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE){
        dens = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute==FALSE) {
        dens = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "OU_linear_sig") {
      A=parameters[1]
      sig2_slope=parameters[2]
      sig2_int=parameters[3]
      sig2=GRAD*sig2_slope + sig2_int
      var = sig2*(1-exp(-2*A*ages))/A
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE){
        dens = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute==FALSE) {
        dens = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "OU_cat") {
      # parameters are in the order sig2, alpha1, alpha2, ... , alphaN
      if(is.null(cats)) {
        stop("Hold up! You haven't included category values in the 'cats' vector")
      }
      if(length(parameters) != length(unique(cats))+1) {
        stop("The number of paramters provided doesn't match the number of categories")
      }
      sig2 = parameters[1]
      alphas = rep(NA, length(ages))
      for(i in 1:length(unique(cats))) alphas[cats==unique(cats)[i]] = parameters[1+i]
      var = sig2*(1-exp(-2*alphas*ages))/alphas
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE){
        dens = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute==FALSE) {
        dens = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "DA_null") {
      A=parameters[1]
      sig2=parameters[2]
      psi=parameters[3]
      u = psi*(1-exp(-A*ages))
      var = sig2*(1-exp(-2*A*ages))/A
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE) {
        dens = log(dnorm(div, mean=u, sd=sqrt(var)) + dnorm(div, mean=-u, sd=sqrt(var)))
      } 
      if(absolute==FALSE) {
        dens = dnorm(div, mean=u, sd=sqrt(var), log=TRUE)
        }
      negLogL = -sum(dens)
    }
    if (model == "DA_linear") { 
      if(length(GRAD) != length(ages)) {
        stop("You haven't provided a gradient value for each pair in the dataset")
      }
      A=parameters[1]
      sig2=parameters[2]
      psi_slope=parameters[3]
      psi_int=parameters[4]
      psi=psi_slope*GRAD + psi_int # vector of peak distances
      #print(psi)
      if(NaN %in% psi) {
      	negLogL = 1e+20
      	} else {
        if(all(psi > 0)) {
          u = psi*(1-exp(-A*ages)) # u is a vector
          var = sig2*(1-exp(-2*A*ages))/A
          if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
          if(absolute==TRUE){
            dens = log(dnorm(div, mean=u, sd=sqrt(var)) + dnorm(div, mean=-u, sd=sqrt(var)))
          } 
          if(absolute==FALSE){
            dens = dnorm(div, mean=u, sd=sqrt(var), log=TRUE)
          }
          negLogL = -sum(dens)
          } 
        if(all(psi > 0) == FALSE) {
          negLogL = 1e+20
          }
      }
    }
    if (model == "DA_cat") {
      # parameters are in the order alpha, sig2, psi1, psi2, ..., psiN
      if(is.null(cats)) {
        stop("Hold up! You haven't included category values in the 'cats' vector")
      }
      if(length(parameters) != length(unique(cats))+2) {
        stop("The number of paramters provided doesn't match the number of categories")
      }
      A = parameters[1]
      sig2 = parameters[2]
      psis = rep(NA, length(ages))
      for(i in 1:length(unique(cats))) psis[cats==unique(cats)[i]] = parameters[2+i]
      u = psis*(1-exp(-A*ages))
      var = sig2*(1-exp(-2*A*ages))/A
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE){
        dens = log(dnorm(div, mean=u, sd=sqrt(var)) + dnorm(div, mean=-u, sd=sqrt(var)))
      } else {
        dens = dnorm(div, mean=u, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "DA_OU") {
      # parameters are in the order alpha, sig2, psi, proportion
      A = parameters[1]
      sig2 = parameters[2]
      psi = parameters[3]
      prop = parameters[4] # this must be limited to be between zero and one!
      u_DA = psi*(1-exp(-A*ages)) # U_OU is zero 
      var = sig2*(1-exp(-2*A*ages))/A #(var is unrelated to theta or psi and thus same for DA and OU)
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE){
        dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var)) + dnorm(div, mean=-u_DA, sd=sqrt(var))
        dens_OU = 2*dnorm(div, mean=0, sd=sqrt(var))
        dens = log(prop*dens_DA + (1-prop)*dens_OU)
      } else {
        dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var))
        dens_OU = dnorm(div, mean=0, sd=sqrt(var))
        dens = log(prop*dens_DA + (1-prop)*dens_OU)
      }
      negLogL = -sum(dens)
    }
    if (model == "DA_BM") {
      # parameters are in the order alpha, sig2, psi, proportion
      A = parameters[1]
      sig2 = parameters[2]
      psi = parameters[3]
      prop = parameters[4] # this must be limited to be between zero and one!
      u_DA = psi*(1-exp(-A*ages)) # U_BM is zero 
      var_DA = sig2*(1-exp(-2*A*ages))/A #(var is unrelated to theta or psi and thus same for DA and OU)
      var_BM = 2*sig2*ages
      if(is.null(me1)==FALSE) {
      	var_DA = var_DA + me1^2 + me2^2
      	var_BM = var_BM + me1^2 + me2^2
      }
      if(absolute==TRUE){
        dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var_DA)) + dnorm(div, mean=-u_DA, sd=sqrt(var_DA))
        dens_BM = 2*dnorm(div, mean=0, sd=sqrt(var_BM))
        dens = log(prop*dens_DA + (1-prop)*dens_BM)
      } else {
        dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var_DA))
        dens_BM = dnorm(div, mean=0, sd=sqrt(var_BM))
        dens = log(prop*dens_DA + (1-prop)*dens_BM)
      }
      negLogL = -sum(dens)
    }
    if (model == "OU_BM") {
      # parameters are in the order alpha, sig2, proportion
      A = parameters[1]
      sig2 = parameters[2]
      prop = parameters[3] # this must be limited to be between zero and one!
      var_OU = sig2*(1-exp(-2*A*ages))/A
      var_BM = 2*sig2*ages
      if(is.null(me1)==FALSE) {
      	var_BM = var_BM + me1^2 + me2^2
      	var_OU = var_OU + me1^2 + me2^2
      }
      if(absolute==TRUE){
        dens_OU = 2*dnorm(div, mean=0, sd=sqrt(var_OU))
        dens_BM = 2*dnorm(div, mean=0, sd=sqrt(var_BM))
        dens = log(prop*dens_OU + (1-prop)*dens_BM)
      } else {
        dens_OU = dnorm(div, mean=0, sd=sqrt(var_OU))
        dens_BM = dnorm(div, mean=0, sd=sqrt(var_BM))
        dens = log(prop*dens_OU + (1-prop)*dens_BM)
      }
      negLogL = -sum(dens)
    }
    if (model == "DA_OU_linear") {
      if(length(GRAD) != length(ages)) {
        stop("You haven't provided a gradient value for each pair in the dataset")
      }
      # parameters are in the order alpha, sig2, psi, prop_slope, prop_int
      A=parameters[1]
      sig2=parameters[2]
      psi=parameters[3]
      prop = parameters[4]*GRAD + parameters[5]
      if(NaN %in% prop) {
        negLogL = 1e+20
      } else {
        if(all(prop <=1 & prop >=0)==FALSE) {
          negLogL=1e+20
        } else {
          u_DA = psi*(1-exp(-A*ages)) # u_BM is zero
          var = sig2*(1-exp(-2*A*ages))/A
          if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
          if(absolute==TRUE){
            dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var)) + dnorm(div, mean=-u_DA, sd=sqrt(var))
            dens_OU = 2*dnorm(div, mean=0, sd=sqrt(var))
            dens = log(prop*dens_DA + (1-prop)*dens_OU)
          } 
          if(absolute==FALSE){
            dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var))
            dens_OU = dnorm(div, mean=0, sd=sqrt(var))
            dens = log(prop*dens_DA + (1-prop)*dens_OU)
          }
          negLogL = -sum(dens)
        }
      }
    }
    if (model == "DA_BM_linear") {
      if(length(GRAD) != length(ages)) {
        stop("You haven't provided a gradient value for each pair in the dataset")
      }
      # parameters are in the order alpha, sig2, psi, prop_slope, prop_int
      A=parameters[1]
      sig2=parameters[2]
      psi=parameters[3]
      prop = parameters[4]*GRAD + parameters[5]
      if(NaN %in% prop) {
        negLogL = 1e+20
      } else {
        if(all(prop <=1 & prop >=0)==FALSE) {
          negLogL=1e+20
        } else {
          u_DA = psi*(1-exp(-A*ages)) # u_BM is zero
          var_DA = sig2*(1-exp(-2*A*ages))/A
          var_BM = 2*sig2*ages
          if(is.null(me1)==FALSE) {
          	var_DA = var_DA + me1^2 + me2^2
          	var_BM = var_BM + me1^2 + me2^2
          }
          if(absolute==TRUE){
            dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var_DA)) + dnorm(div, mean=-u_DA, sd=sqrt(var_DA))
            dens_BM = 2*dnorm(div, mean=0, sd=sqrt(var_BM))
            dens = log(prop*dens_DA + (1-prop)*dens_BM)
          } 
          if(absolute==FALSE){
            dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var_DA))
            dens_BM = dnorm(div, mean=0, sd=sqrt(var_BM))
            dens = log(prop*dens_DA + (1-prop)*dens_BM)
          }
          negLogL = -sum(dens)
        }
      }
    }
    if (model == "OU_BM_linear") {
      if(length(GRAD) != length(ages)) {
        stop("You haven't provided a gradient value for each pair in the dataset")
      }
      # parameters are in the order alpha, sig2, psi, prop_slope, prop_int
      A=parameters[1]
      sig2=parameters[2]
      prop = parameters[3]*GRAD + parameters[4]
      if(NaN %in% prop) {
        negLogL = 1e+20
      } else {
        if(all(prop <=1 & prop >=0)==FALSE) {
          negLogL=1e+20
        } else {
          var_OU = sig2*(1-exp(-2*A*ages))/A
          var_BM = 2*sig2*ages
          if(is.null(me1)==FALSE) {
          	var_OU = var_OU + me1^2 + me2^2
          	var_BM = var_BM + me1^2 + me2^2
          }
          if(absolute==TRUE){
            dens_OU = 2*dnorm(div, mean=0, sd=sqrt(var_OU))
            dens_BM = 2*dnorm(div, mean=0, sd=sqrt(var_BM))
            dens = log(prop*dens_OU + (1-prop)*dens_BM)
          } 
          if(absolute==FALSE){
            dens_OU = dnorm(div, mean=0, sd=sqrt(var_OU))
            dens_BM = dnorm(div, mean=0, sd=sqrt(var_BM))
            dens = log(prop*dens_OU + (1-prop)*dens_BM)
          }
          negLogL = -sum(dens)
        }
      }
    }
    if (model == "DA_OU_cat") {
      # parameters are in the order alpha, sig2, psi, prop1, prop2, ..., propN
      if(is.null(cats)) {
        stop("Hold up! You haven't included category values in the 'cats' vector")
      }
      if(length(parameters) != length(unique(cats))+3) {
        stop("The number of parameters provided doesn't match the number of categories")
      }
      A = parameters[1]
      sig2 = parameters[2]
      psi = parameters[3]
      prop = rep(NA, length(ages))
      for(i in 1:length(unique(cats))) prop[cats==unique(cats)[i]] = parameters[3+i]
      u_DA = psi*(1-exp(-A*ages)) # u_OU is zero
      var = sig2*(1-exp(-2*A*ages))/A # vars are the same for DA and OU models
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      if(absolute==TRUE){
        dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var)) + dnorm(div, mean=-u_DA, sd=sqrt(var))
        dens_OU = 2*dnorm(div, mean=0, sd=sqrt(var))
        dens = log(prop*dens_DA + (1-prop)*dens_OU) # I think the issue is here
      } else {
        dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var))
        dens_OU = dnorm(div, mean=0, sd=sqrt(var))
        dens = log(prop*dens_DA + (1-prop)*dens_OU)
      }
      negLogL = -sum(dens)
    }
    if (model == "DA_BM_cat") {
      # parameters are in the order alpha, sig2, psi, prop1, prop2, ..., propN
      if(is.null(cats)) {
        stop("Hold up! You haven't included category values in the 'cats' vector")
      }
      if(length(parameters) != length(unique(cats))+3) {
        stop("The number of parameters provided doesn't match the number of categories")
      }
      A = parameters[1]
      sig2 = parameters[2]
      psi = parameters[3]
      prop = rep(NA, length(ages))
      for(i in 1:length(unique(cats))) prop[cats==unique(cats)[i]] = parameters[3+i]
      u_DA = psi*(1-exp(-A*ages)) # u_OU is zero
      var_DA = sig2*(1-exp(-2*A*ages))/A
      var_BM = 2*sig2*ages
      if(is.null(me1)==FALSE) {
      	var_DA = var_DA + me1^2 + me2^2
      	var_BM = var_BM + me1^2 + me2^2
      }
      if(absolute==TRUE){
        dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var_DA)) + dnorm(div, mean=-u_DA, sd=sqrt(var_DA))
        dens_BM = 2*dnorm(div, mean=0, sd=sqrt(var_BM))
        dens = log(prop*dens_DA + (1-prop)*dens_BM)
      } else {
        dens_DA = dnorm(div, mean=u_DA, sd=sqrt(var_DA))
        dens_BM = dnorm(div, mean=0, sd=sqrt(var_BM))
        dens = log(prop*dens_DA + (1-prop)*dens_BM)
      }
      negLogL = -sum(dens)
    }
    if (model == "OU_BM_cat") {
      # parameters are in the order alpha, sig2, prop1, prop2, ..., propN
      if(is.null(cats)) {
        stop("Hold up! You haven't included category values in the 'cats' vector")
      }
      if(length(parameters) != length(unique(cats))+2) {
        stop("The number of parameters provided doesn't match the number of categories")
      }
      A = parameters[1]
      sig2 = parameters[2]
      prop = rep(NA, length(ages))
      for(i in 1:length(unique(cats))) prop[cats==unique(cats)[i]] = parameters[2+i]
      var_OU = sig2*(1-exp(-2*A*ages))/A
      var_BM = 2*sig2*ages
      if(is.null(me1)==FALSE) {
      	var_OU = var_OU + me1^2 + me2^2
      	var_BM = var_BM + me1^2 + me2^2
      }
      if(absolute==TRUE){
        dens_OU = 2*dnorm(div, mean=0, sd=sqrt(var_OU))
        dens_BM = 2*dnorm(div, mean=0, sd=sqrt(var_BM))
        dens = log(prop*dens_OU + (1-prop)*dens_BM)
      } else {
        dens_OU = dnorm(div, mean=0, sd=sqrt(var_OU))
        dens_BM = dnorm(div, mean=0, sd=sqrt(var_BM))
        dens = log(prop*dens_OU + (1-prop)*dens_BM)
      }
      negLogL = -sum(dens)
    }
    if (model == "DA_wt") {
      # assumes a shared wait time in units of M.Y. since speciation
      # note: code for likelihood doesn't look like the piecewise likelihood function, but it is equivalent
      A=parameters[1]
      sig2=parameters[2]
      psi1=parameters[3]
      psi2=parameters[4]
      wt=parameters[5]
      DIST_wt=div[which(ages > wt)] # subset of divergence values = sisters older than the wait time 
      DIST_nwt = div[which(ages <= wt)] # subset of divergence values = sisters younger than wait time
      ED=c(DIST_wt, DIST_nwt) # re-combine divergence value vector
      TIME_wt=ages[ages>wt] # subset of sister pair ages = sisters older than wait time
      TIME_nwt = ages[ages<=wt] # subset of sister pair ages = sisters younger than wait time
      TIME=c(TIME_wt, TIME_nwt) # re-combine ages vector
      var = sig2*(1-exp(-2*A*TIME))/A # vector of variances
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      uwt = exp(-A*TIME_wt)*(psi1*(exp(A*wt)-1)+psi2*(exp(A*TIME_wt)-exp(A*wt)))
      unwt = psi1*(1-exp(-A*TIME_nwt))
      u=c(uwt, unwt)
      if(absolute==TRUE){
        dens = log(dnorm(ED, mean=u, sd=sqrt(var)) + dnorm(ED, mean=-u, sd=sqrt(var)))
      } else {
        dens = dnorm(ED, mean=u, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "DA_bp") {
      if(is.null(bp)) {
          stop("Hold up: You haven't provided a 'breakpoint' vector\nReminder: breakpoint value must == 0 for all single-epoch pairs")
      }
      for(i in 1:length(bp)) {
        if(bp[i] > ages[i]) stop("Pls check input data alignment. The breakpoint time for each pair must be less than the age of the pair")
      }
      # this function requires a vector of breakpoints ('bp')
      # bp has to be set to zero for all sisters that did not experience a breakpoint
      A=parameters[1]
      sig2=parameters[2]
      psi1=parameters[3]
      psi2=parameters[4]
      bpnz=bp[bp>0]
      DIST_bp=div[bp>0] # subset of divergence values = sisters that have experienced an epoch shift
      DIST_nbp=div[bp==0] # subset of divergence values = sisters that have not experienced an epoch shift
      ED=c(DIST_bp, DIST_nbp) # re-combine divergence value vector
      TIME_bp=ages[bp>0] # subset of sister pair ages = sisters that have experienced an epoch shift
      TIME_nbp = ages[bp==0] # subset of sister pair ages = sisters that have not experienced an epoch shift
      TIME=c(TIME_bp, TIME_nbp) # re-combine ages vector
      var = sig2*(1-exp(-2*A*TIME))/A
      if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
      ubp = exp(-A*TIME_bp)*(psi1*(exp(A*bpnz)-1)+psi2*(exp(A*TIME_bp)-exp(A*bpnz))) 
      unbp = psi1*(1-exp(-A*TIME_nbp)) 
      u=c(ubp, unbp)
      if(absolute==TRUE){
        dens = log(dnorm(ED, mean=u, sd=sqrt(var)) + dnorm(ED, mean=-u, sd=sqrt(var)))
      } else {
        dens = dnorm(ED, mean=u, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(dens)
    }
    if (model == "DA_wt_linear") {
      if(length(GRAD) != length(ages)) {
        stop("You haven't provided a gradient value for each pair in the dataset")
      }
      if(length(parameters) != 7) {
        stop("You haven't provided the right number of starting parameters for the likelihood search")
      }
      A=parameters[1]
      sig2=parameters[2]
      wt = parameters[3]
      psi1_slope=parameters[4]
      psi1_int=parameters[5]
      psi2_slope=parameters[6]
      psi2_int=parameters[7]
      if(any(ages < wt)) {
        psi1_wt=GRAD[which(ages > wt)]*psi1_slope + psi1_int
        psi1_nwt=GRAD[which(ages <= wt)]*psi1_slope + psi1_int
        psi2=GRAD[which(ages > wt)]*psi2_slope + psi2_int
      } else {
        negLogL= 1e+20
      }
      if(NaN %in% c(psi1_wt, psi1_nwt, psi2)) {
      	negLogL = 1e+20
      }
      #if(any(c(psi1_wt, psi1_nwt, psi2) <= 0)) {
      #  negLogL = 1e+20
      #}
      if(all(c(psi1_wt, psi1_nwt, psi2) > 0)) {
        dist_wt=div[which(ages > wt)] # subset of divergence values = sisters older than the wait time 
        dist_nwt = div[which(ages <= wt)] # subset of divergence values = sisters younger than wait time
        dist=c(dist_wt, dist_nwt) # re-combine divergence value vector
        time_wt=ages[ages>wt] # subset of sister pair ages = sisters older than wait time
        time_nwt = ages[ages<=wt] # subset of sister pair ages = sisters younger than wait time
        TIME=c(time_wt, time_nwt) # re-combine ages vector
        var = sig2*(1-exp(-2*A*TIME))/A
        if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
        uwt = exp(-A*time_wt)*(psi1_wt*(exp(A*wt)-1)+psi2*(exp(A*time_wt)-exp(A*wt)))
        unwt = psi1_nwt*(1-exp(-A*time_nwt))
        u=c(uwt, unwt)
        if(absolute==TRUE){
          dens = log(dnorm(dist, mean=u, sd=sqrt(var)) + dnorm(dist, mean=-u, sd=sqrt(var)))
        } else {
          dens = dnorm(dist, mean=u, sd=sqrt(var), log=TRUE)
        }
      negLogL = -sum(dens)
      }
    }
    if (model == "DA_bp_linear") {
      if(is.null(bp)) {
          stop("Hold up: You haven't provided a 'breakpoint' vector\nReminder: breakpoint value must == 0 for all single-epoch pairs")
      }
      for(i in 1:length(bp)) {
        if(bp[i] > ages[i]) stop("Pls check input data alignment. The breakpoint time for each pair must be less than the age of the pair")
      }
      # this function requires a vector of breakpoints ('bp')
      # bp has to be set to zero for all single-epoch sisters
      A=parameters[1]
      sig2=parameters[2]
      psi1_slope=parameters[3]
      psi1_int=parameters[4]
      psi2_slope=parameters[5]
      psi2_int=parameters[6]
      psi1_bp=GRAD[bp>0]*psi1_slope + psi1_int
      psi1_nbp=GRAD[bp==0]*psi1_slope +psi1_int
      psi2=GRAD[bp>0]*psi2_slope + psi2_int
      if(NaN %in% c(psi1_bp, psi1_nbp, psi2)) {
      	negLogL = 1e+20
      }
      if(all(c(psi1_bp, psi1_nbp, psi2) > 0) == FALSE) {
      	negLogL = 1e+20
      }
      if(all(c(psi1_bp, psi1_nbp, psi2) > 0)) {
        bpnz=bp[bp>0] # a vector with all the non-zero breakpoints
        dist_bp=div[bp>0] # subset of divergence values = sisters that have experienced an epoch shift
        dist_nbp=div[bp==0] # subset of divergence values = sisters that have not experienced an epoch shift
        dist=c(dist_bp, dist_nbp) # re-combine divergence value vector
        time_bp=ages[bp>0] # subset of sister pair ages = sisters that have experienced an epoch shift
        time_nbp = ages[bp==0] # subset of sister pair ages = sisters that have not experienced an epoch shift
        TIME=c(time_bp, time_nbp) # re-combine ages vector
        var = sig2*(1-exp(-2*A*TIME))/A 
        if(is.null(me1)==FALSE) var = var + me1^2 + me2^2
        ubp = exp(-A*time_bp)*(psi1_bp*(exp(A*bpnz)-1)+psi2*(exp(A*time_bp)-exp(A*bpnz)))
        unbp = psi1_nbp*(1-exp(-A*time_nbp)) 
        u=c(ubp, unbp)
        if(absolute==TRUE) {
          dens = log(dnorm(dist, mean=u, sd=sqrt(var)) + dnorm(dist, mean=-u, sd=sqrt(var)))
        } else {
          dens = dnorm(dist, mean=u, sd=sqrt(var), log=TRUE)
        }
        negLogL = -sum(dens)
      }
    }
  if(is.nan(negLogL)) negLogL=1e+20
return(negLogL)
}

