likelihood_func <-
function (parameters, model, div, ages, cats=NULL, GRAD=NULL, bp=NULL, absolute=TRUE) {
    
    if (model == "BM_null") {
      B=parameters[1]
      var=2*B*ages
      if(absolute==TRUE) {
        kk3 = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute == FALSE) {
        kk3 = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(kk3)
    }
    if (model == "BM_linear") {
      B_slope=parameters[1]
      B_int=parameters[2]
      B=GRAD*B_slope + B_int
      var=2*B*ages
      if(absolute==TRUE) {
        kk3 = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute == FALSE) {
        kk3 = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(kk3)
    }
    if (model == "OU_null") {
      A=parameters[1]
      B=parameters[2]
      var = B*(1-exp(-2*A*ages))/A
      if(absolute==TRUE){
        kk3 = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute==FALSE) {
        kk3 = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(kk3)
    }
    if (model == "OU_linear") {
      A_int=parameters[1]
      A_slope=parameters[2]
      B=parameters[3]
      A=GRAD*A_slope + A_int
      var = B*(1-exp(-2*A*ages))/A
      if(absolute==TRUE){
        kk3 = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute==FALSE) {
        kk3 = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(kk3)
    }
    if (model == "OU_linear_sig") {
      A=parameters[1]
      B_slope=parameters[2]
      B_int=parameters[3]
      B=GRAD*B_slope + B_int
      var = B*(1-exp(-2*A*ages))/A
      if(absolute==TRUE){
        kk3 = log(2*dnorm(div, mean=0, sd=sqrt(var)))
      }
      if(absolute==FALSE) {
        kk3 = dnorm(div, mean=0, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(kk3)
    }
    if (model == "DA_null") {
      A=parameters[1]
      B=parameters[2]
      psi=parameters[3]
      u = psi*(1-exp(-A*ages))
      var = B*(1-exp(-2*A*ages))/A
      if(absolute==TRUE) {
        kk3 = log(dnorm(div, mean=u, sd=sqrt(var)) + dnorm(div, mean=-u, sd=sqrt(var)))
      } 
      if(absolute==FALSE) {
        kk3 = dnorm(div, mean=u, sd=sqrt(var), log=TRUE)
        }
      negLogL = -sum(kk3)
      }
    if (model == "DA_linear") { 
      if(length(GRAD) != length(ages)) {
        stop("You haven't provided a gradient value for each pair in the dataset")
      }
      A=parameters[1]
      B=parameters[2]
      psi_slope=parameters[3]
      psi_int=parameters[4]
      psi=psi_slope*GRAD + psi_int # vector of peak distances
      #print(psi)
      if(NaN %in% psi) {
      	negLogL = 1e+20
      	} else {
        if(all(psi > 0)) {
          u = psi*(1-exp(-A*ages)) # u is a vector
          var = B*(1-exp(-2*A*ages))/A
          if(absolute==TRUE){
            kk3 = log(dnorm(div, mean=u, sd=sqrt(var)) + dnorm(div, mean=-u, sd=sqrt(var)))
          } 
          if(absolute==FALSE){
            kk3 = dnorm(div, mean=u, sd=sqrt(var), log=TRUE)
          }
          negLogL = -sum(kk3)
          } 
        if(all(psi > 0) == FALSE) {
          negLogL = 1e+20
          }
      }
    }
    if (model == "DA_cat") {
      if(is.null(cats)) {
        stop("Hold up! You haven't included category values in the 'cats' vector")
      }
      A = parameters[1]
      B = parameters[2]
      psi1=parameters[3]
      psi2=parameters[4]
      DIST1 = div[cats==0]
      DIST2 = div[cats==1]
      TIME1 = ages[cats==0]
      TIME2 = ages[cats==1]
      TIME = c(TIME1, TIME2)
      ED = c(DIST1, DIST2)
      u1 = psi1*(1-exp(-A*TIME1))
      u2 = psi2*(1-exp(-A*TIME2))
      u = c(u1, u2)
      if(length(parameters==5)) {
        psi3 = parameters[5]
        DIST3 = div[cats==2]
        ED = c(ED, DIST3)
        TIME3 = ages[cats==2]
        TIME = c(TIME, TIME3)
        u3 = psi3*(1-exp(-A*TIME3))
        u = c(u, u3)
      }
      var = B*(1-exp(-2*A*TIME))/A
      if(absolute==TRUE){
        kk3 = log(dnorm(ED, mean=u, sd=sqrt(var)) + dnorm(ED, mean=-u, sd=sqrt(var)))
      } else {
        kk3 = dnorm(ED, mean=u, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(kk3)
      }
    if (model == "DA_wt") {
      # assumes a shared wait time in units of M.Y. since speciation
      # note: code for likelihood doesn't look like the piecewise likelihood function, but it is equivalent
      A=parameters[1]
      B=parameters[2]
      psi1=parameters[3]
      psi2=parameters[4]
      wt=parameters[5]
      DIST_wt=div[which(ages > wt)] # subset of divergence values = sisters older than the wait time 
      DIST_nwt = div[which(ages <= wt)] # subset of divergence values = sisters younger than wait time
      ED=c(DIST_wt, DIST_nwt) # re-combine divergence value vector
      TIME_wt=ages[ages>wt] # subset of sister pair ages = sisters older than wait time
      TIME_nwt = ages[ages<=wt] # subset of sister pair ages = sisters younger than wait time
      TIME=c(TIME_wt, TIME_nwt) # re-combine ages vector
      var = B*(1-exp(-2*A*TIME))/A # vector of variances
      uwt = exp(-A*TIME_wt)*(psi1*(exp(A*wt)-1)+psi2*(exp(A*TIME_wt)-exp(A*wt)))
      unwt = psi1*(1-exp(-A*TIME_nwt))
      u=c(uwt, unwt)
      if(absolute==TRUE){
        kk3 = log(dnorm(ED, mean=u, sd=sqrt(var)) + dnorm(ED, mean=-u, sd=sqrt(var)))
      } else {
        kk3 = dnorm(ED, mean=u, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(kk3)
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
      B=parameters[2]
      psi1=parameters[3]
      psi2=parameters[4]
      bpnz=bp[bp>0]
      DIST_bp=div[bp>0] # subset of divergence values = sisters that have experienced an epoch shift
      DIST_nbp=div[bp==0] # subset of divergence values = sisters that have not experienced an epoch shift
      ED=c(DIST_bp, DIST_nbp) # re-combine divergence value vector
      TIME_bp=ages[bp>0] # subset of sister pair ages = sisters that have experienced an epoch shift
      TIME_nbp = ages[bp==0] # subset of sister pair ages = sisters that have not experienced an epoch shift
      TIME=c(TIME_bp, TIME_nbp) # re-combine ages vector
      var = B*(1-exp(-2*A*TIME))/A
      ubp = exp(-A*TIME_bp)*(psi1*(exp(A*bpnz)-1)+psi2*(exp(A*TIME_bp)-exp(A*bpnz))) 
      unbp = psi1*(1-exp(-A*TIME_nbp)) 
      u=c(ubp, unbp)
      if(absolute==TRUE){
        kk3 = log(dnorm(ED, mean=u, sd=sqrt(var)) + dnorm(ED, mean=-u, sd=sqrt(var)))
      } else {
        kk3 = dnorm(ED, mean=u, sd=sqrt(var), log=TRUE)
      }
      negLogL = -sum(kk3)
    }
    if (model == "DA_wt_linear") {
      if(length(GRAD) != length(ages)) {
        stop("You haven't provided a gradient value for each pair in the dataset")
      }
      if(length(parameters) != 7) {
        stop("You haven't provided the right number of starting parameters for the likelihood search")
      }
      A=parameters[1]
      B=parameters[2]
      wt = parameters[3]
      psi1_slope=parameters[4]
      psi1_int=parameters[5]
      psi2_slope=parameters[6]
      psi2_int=parameters[7]
      psi1_wt=GRAD[which(ages > wt)]*psi1_slope + psi1_int
      psi1_nwt=GRAD[which(ages <= wt)]*psi1_slope + psi1_int
      psi2=GRAD[which(ages > wt)]*psi2_slope + psi2_int
      if(NaN %in% c(psi1_wt, psi1_nwt, psi2)) {
      	negLogL = 1e+20
      }
      if(all(c(psi1_wt, psi1_nwt, psi2) > 0) == FALSE) {
        negLogL = 1e+20
      }
      if(all(c(psi1_wt, psi1_nwt, psi2) > 0)) {
        dist_wt=div[which(ages > wt)] # subset of divergence values = sisters older than the wait time 
        dist_nwt = div[which(ages <= wt)] # subset of divergence values = sisters younger than wait time
        dist=c(dist_wt, dist_nwt) # re-combine divergence value vector
        time_wt=ages[ages>wt] # subset of sister pair ages = sisters older than wait time
        time_nwt = ages[ages<=wt] # subset of sister pair ages = sisters younger than wait time
        TIME=c(time_wt, time_nwt) # re-combine ages vector
        var = B*(1-exp(-2*A*TIME))/A
        uwt = exp(-A*time_wt)*(psi1_wt*(exp(A*wt)-1)+psi2*(exp(A*time_wt)-exp(A*wt)))
        unwt = psi1_nwt*(1-exp(-A*time_nwt))
        u=c(uwt, unwt)
        if(absolute==TRUE){
          kk3 = log(dnorm(dist, mean=u, sd=sqrt(var)) + dnorm(dist, mean=-u, sd=sqrt(var)))
        } else {
          kk3 = dnorm(dist, mean=u, sd=sqrt(var), log=TRUE)
        }
      negLogL = -sum(kk3)
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
      B=parameters[2]
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
        var = B*(1-exp(-2*A*TIME))/A 
        ubp = exp(-A*time_bp)*(psi1_bp*(exp(A*bpnz)-1)+psi2*(exp(A*time_bp)-exp(A*bpnz)))
        unbp = psi1_nbp*(1-exp(-A*time_nbp)) 
        u=c(ubp, unbp)
        if(absolute==TRUE) {
          kk3 = log(dnorm(dist, mean=u, sd=sqrt(var)) + dnorm(dist, mean=-u, sd=sqrt(var)))
        } else {
          kk3 = dnorm(dist, mean=u, sd=sqrt(var), log=TRUE)
        }
        negLogL = -sum(kk3)
      }
    }
  if(is.nan(negLogL)) negLogL=1e+20
return(negLogL)
}

