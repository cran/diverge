expected_val <-
function (model, sig2, alpha = NULL, psi = NULL, time_span = c(0, 10), 
  quantile = FALSE, plot = FALSE, labels = TRUE, bm_col = "darkgoldenrod1", 
  ou_col = "firebrick2", da_col = "navy", exval_lwd = 5, ylim = NULL,...) 
{
    if(model %in% c("BM_null", "OU_null", "DA_null", "DA_OU", "DA_BM", "OU_BM")==FALSE) {
      stop("Spell check: you've entered a model that doesn't match the models accepted by this function")
    }
    if (length(time_span) == 1) {
        time <- time_span
    }
    if (length(time_span) > 2) {
        time <- time_span
    }
    if (length(time_span) == 2) {
        time = seq(0, time_span[2], (time_span[2]/1000))
    }
    if(model == "BM_null") {
        V <- sig2 * time * 2
        u = rep(0, length(time))
        exdiv <- ((2 * V)/pi)^0.5
    }
    if(model == "OU_null") {
        V <- (sig2/(alpha)) * (1 - exp(-2 * alpha * time))
        u = rep(0, length(time))
        exdiv <- ((2 * V)/pi)^0.5
    }
    if(model == "DA_null") {
        V <- (sig2/(alpha)) * (1 - exp(-2 * alpha * time))
        u <- psi * (1 - exp(-alpha * time))
        exdiv <- sqrt((2 * V)/pi) * exp(-(u^2)/(2 * V)) + u * (2 * pnorm((u/sqrt(2 * V)) * sqrt(2)) - 1)
    }
    if(model == "DA_OU") {
        V = V2 = (sig2/(alpha)) * (1 - exp(-2 * alpha * time))
        u = rep(0, length(time))
        u2 = psi * (1 - exp(-alpha * time))
        exdiv <- ((2 * V)/pi)^0.5
        exdiv2 <- sqrt((2 * V2)/pi) * exp(-(u2^2)/(2 * V2)) + u2 * (2 * pnorm((u2/sqrt(2 * V2)) * sqrt(2)) - 1)
    }
    if(model == "DA_BM") {
        V <- sig2 * time * 2
        u <- rep(0, length(time))
        exdiv <- ((2 * V)/pi)^0.5
        V2 <- (sig2/(alpha)) * (1 - exp(-2 * alpha * time))
        u2 <- psi * (1 - exp(-alpha * time))
        exdiv2 <- sqrt((2 * V2)/pi) * exp(-(u2^2)/(2 * V2)) + u2 * (2 * pnorm((u2/sqrt(2 * V2)) * sqrt(2)) - 1)
    }
    if(model == "OU_BM") {
        V <- sig2 * time * 2
        u = u2 = rep(0, length(time))
        V2 <- (sig2/(alpha)) * (1 - exp(-2 * alpha * time))
        exdiv <- ((2 * V)/pi)^0.5
        exdiv2 <- ((2 * V2)/pi)^0.5
    }
    if (quantile == TRUE) {
        qs = c(0.025, seq(0.1, 0.9, 0.1), 0.975) # the ends now give you confidence intervals
        #quantiles = lapply(qs, FUN = truncnorm::qtruncnorm, a = 0, mean = u, sd = sqrt(V))
        #quantiles = lapply(qs, function(x) quantile(abs(rnorm(n=10000, mean=u, sd=sqrt(V))), probs=x))
        # NOTE: FOR SMOOTHER CONFIDENCE INTERVAL EDGES, SET ncol=100000 in the appropriate spots below
        dists = matrix(NA, nrow=length(V), ncol=10000)
        for(i in 1:length(time)) dists[i,] = abs(rnorm(10000,mean=u[i], sd=sqrt(V[i])))
        quantiles = matrix(NA, nrow=length(V), ncol=length(qs))
        for(i in 1:length(time)) quantiles[i,] = quantile(dists[i,], probs=qs)
        res = cbind(time, exdiv, quantiles)
        #res <- matrix(NA, length(time), 13)
        colnames(res) <- c("time", "Expectation", "q025", "q10","q20", 
            "q30", "q40", "q50", "q60", "q70", "q80", "q90", 
            "q975")
        #res[, 1] <- time
        #res[, 2] <- exdiv
        #for (i in 3:13) res[, i] = quantiles[[i - 2]]
        if(model %in% c("DA_OU", "DA_BM", "OU_BM")) {
            #quantiles2 = lapply(qs, FUN = truncnorm::qtruncnorm, a = 0, mean = u2, sd = sqrt(V2))
            #quantiles2 = lapply(qs, function(x) quantile(abs(rnorm(n=100000, mean=u2, sd=sqrt(V2))), probs=x))
            dists2 = matrix(NA, nrow=length(V2), ncol=10000)
            for(i in 1:length(time)) dists2[i,] = abs(rnorm(10000,mean=u2[i], sd=sqrt(V2[i])))
            quantiles2 = matrix(NA, nrow=length(V2), ncol=length(qs))
            for(i in 1:length(time)) quantiles2[i,] = quantile(dists2[i,], probs=qs)
            res2=cbind(exdiv2,quantiles2)
            #res2 = res[,-1]
            colnames(res2) = paste("model2_", colnames(res)[-1], sep="")
            #res2[,1] <- exdiv2
            #for (i in 2:12) res2[, i] = quantiles2[[i - 1]]
            res = cbind(res, res2)
        }
    } else {
        res <- matrix(NA, length(time), 2)
        colnames(res) <- c("time", "Expected_Abs_Diff")
        res[, 1] <- time
        res[, 2] <- exdiv
        if(model %in% c("DA_OU", "DA_BM", "OU_BM")) {
          res = cbind(res, exdiv2)
          colnames(res)[3] = "Expected_Abs_Diff_model2"
        }
    }
    if (plot == TRUE) {
        if (quantile == TRUE) {
            x = length(quantiles)
            time.plot = c(time, sort(time, decreasing=T))
            y.m1 <- c(res[,"q025"], sort(res[,"q975"], decreasing=T))
            if(model %in% c("DA_OU", "DA_BM", "OU_BM")) {
              y.m2 <- c(res[,"model2_q025"], sort(res[,"model2_q975"], decreasing=T))
              if(is.null(ylim)) ylim=c(0, max(c(quantiles[[x]],quantiles2[[x]]),na.rm=TRUE)*1.2)
              if(model=="DA_OU") {
                plot(res[,"model2_Expectation"] ~ time, type="l", ann=FALSE, col=da_col, lwd=exval_lwd, ylim=ylim, ...)
                lines(res[,"Expectation"] ~ time, lty=1, col=ou_col, lwd=exval_lwd)
                polygon(time.plot, y.m2, col = adjustcolor(da_col, alpha.f=0.2), border = NA)
                polygon(time.plot, y.m1, col = adjustcolor(ou_col, alpha.f=0.2), border = NA)
              }
              if(model=="OU_BM") {
                plot(res[,"Expectation"] ~ time, type="l", ann=FALSE, col=bm_col, lwd=exval_lwd, ylim=ylim, ...)
                lines(res[,"model2_Expectation"] ~ time, lty=1, col=ou_col, lwd=exval_lwd)
                polygon(time.plot, y.m2, col = adjustcolor(ou_col, alpha.f=0.2), border = NA)
                polygon(time.plot, y.m1, col = adjustcolor(bm_col, alpha.f=0.2), border = NA)                 
              }
              if(model=="DA_BM") {
                if(max(quantiles[[x]], na.rm=T) > max(quantiles2[[x]], na.rm=T)) {
                  plot(res[,"Expectation"] ~ time, type="l", ann=FALSE, col=bm_col, lwd=exval_lwd, ylim=ylim, ...)
                  lines(res[,"model2_Expectation"] ~ time, lty=1, col=da_col, lwd=exval_lwd) 
                } else {
                  plot(res[,"model2_Expectation"] ~ time, type="l", ann=FALSE, col=da_col, lwd=exval_lwd, ylim=ylim, ...)
                  lines(res[,"Expectation"] ~ time, lty=1, col=bm_col, lwd=exval_lwd) 
                }
                polygon(time.plot, y.m1, col = adjustcolor(bm_col, alpha.f=0.2), border = NA) 
                polygon(time.plot, y.m2, col = adjustcolor(da_col, alpha.f=0.2), border = NA)     
              }
            } else {
              if(is.null(ylim)) ylim=c(0, max(quantiles[[x]],na.rm=TRUE)*1.2)
              if(model=="BM_null") {
                plot(exdiv ~ time, type="l", ann=FALSE, col=bm_col, lwd=exval_lwd, ylim=ylim, ...)
                polygon(time.plot, y.m1, col = adjustcolor(bm_col, alpha.f=0.2), border = NA)
              } 
              if(model=="OU_null") {
                plot(exdiv ~ time, type="l", ann=FALSE, col=ou_col, lwd=exval_lwd, ylim=ylim, ...)
                polygon(time.plot, y.m1, col = adjustcolor(ou_col, alpha.f=0.2), border = NA)	
              }
              if(model=="DA_null") {
                plot(exdiv ~ time, type="l", ann=FALSE, col=da_col, lwd=exval_lwd, ylim=ylim, ...)
                polygon(time.plot, y.m1, col = adjustcolor(da_col, alpha.f=0.2), border = NA)	
              }
            }
        }
        if (quantile == FALSE) {
          if(model %in% c("DA_OU", "OU_BM", "DA_BM")) {
            if(is.null(ylim)) ylim=c(0, max(c(exdiv,exdiv2),na.rm=TRUE)*1.2)
            } else {
            if(is.null(ylim)) ylim=c(0, max(exdiv,na.rm=TRUE)*1.2)
          }
          if(model=="DA_OU") {
            plot(exdiv2 ~ time, type="l", ann=FALSE, col=da_col, lwd=exval_lwd, ylim=ylim, ...)
            lines(exdiv ~ time, col=ou_col, lty=1, lwd=exval_lwd)
          }
          if(model=="OU_BM") {
            plot(exdiv ~ time, type="l", ann=FALSE, col=bm_col, lwd=exval_lwd, ylim=ylim, ...)
            lines(exdiv2 ~ time, col=ou_col, lty=1, lwd=exval_lwd)
          }
          if(model=="DA_BM") {
            if(max(exdiv,na.rm=TRUE) > max(exdiv2, na.rm=TRUE)) {
              plot(exdiv ~ time, type="l", ann=FALSE, col=bm_col, lwd=exval_lwd, ylim=ylim, ...)
              lines(exdiv2 ~ time, col=da_col, lty=1, lwd=exval_lwd)
            } else {
              plot(exdiv2 ~ time, type="l", ann=FALSE, col=da_col, lwd=exval_lwd, ylim=ylim, ...)
              lines(exdiv ~ time, col=bm_col, lty=1, lwd=exval_lwd)
            }
          }
          if(model=="DA_null") plot(exdiv ~ time, type="l", ann=FALSE, col=da_col, lwd=exval_lwd, ylim=ylim, ...)
          if(model=='OU_null') plot(exdiv ~ time, type="l", ann=FALSE, col=ou_col, lwd=exval_lwd, ylim=ylim, ...)
          if(model=="BM_null") plot(exdiv ~ time, type="l", ann=FALSE, col=bm_col, lwd=exval_lwd, ylim=ylim, ...)
        }
        if (labels == TRUE) {
            title(xlab = "Age of Lineage Pair")
            title(ylab = "Trait Divergence")
        }
    }
    return(res)
}