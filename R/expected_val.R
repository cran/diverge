expected_val <-
function (model, sig2, alpha = NULL, psi = NULL, psi2 = NULL, wt = NULL, 
    time_span = c(0, 10), quantile = FALSE, plot = FALSE, labels = TRUE, 
    exval_col = "blue", quant_col = "black", exval_col2 = "springgreen2", quant_col2 = "black",
    exval_lwd = 2, quant_lwd = 1, ...) 
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
        time = seq(0, time_span[2], (time_span[2]/10000))
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
    if(model == "DA_wt") {
        if (is.null(wt)) stop("You need to enter a time (wt) at which the shift in epoch occurs")
        if (wt > max(time)) stop("You've chosen an epoch shift that occurs after the max time in this calculation")
        TIME1 = time[0:which(time == wt)]
        TIME2 = time[(which(time == wt) + 1):length(time)]
        u1 = psi * (1 - exp(-alpha * TIME1))
        u2 = exp(-alpha * TIME2) * (psi * (exp(alpha * wt) - 1) + psi2 * (exp(alpha * TIME2) - exp(alpha *  wt)))
        u = c(u1, u2)
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
        qs = c(seq(0.1, 0.9, 0.1), 0.95)
        quantiles = lapply(qs, FUN = truncnorm::qtruncnorm, a = 0, mean = u, sd = sqrt(V))
        res <- matrix(NA, length(time), 12)
        colnames(res) <- c("time", "Expectation", "q10", "q20", 
            "q30", "q40", "q50", "q60", "q70", "q80", "q90", 
            "q95")
        res[, 1] <- time
        res[, 2] <- exdiv
        for (i in 3:12) res[, i] = quantiles[[i - 2]]
        if(model %in% c("DA_OU", "DA_BM", "OU_BM")) {
            quantiles2 = lapply(qs, FUN = truncnorm::qtruncnorm, a = 0, mean = u2, sd = sqrt(V2))
            res2 = res[,-1]
            colnames(res2) = paste("m2_", colnames(res)[-1], sep="")
            res2[,1] <- exdiv2
            for (i in 2:11) res2[, i] = quantiles2[[i - 1]]
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
            if(model %in% c("DA_OU", "DA_BM", "OU_BM")) {
                if(max(quantiles[[x]], na.rm=TRUE) > max(quantiles2[[x]], na.rm=TRUE)) {
                  plot(quantiles[[x]] ~ time, type = "l", ann = FALSE, col = quant_col, lwd = quant_lwd, ...)
                  lines(quantiles2[[x]] ~ time, lty = 3, col = quant_col2, lwd = quant_lwd)
                } else {
                  plot(quantiles2[[x]] ~ time, type = "l", lty=3, ann = FALSE, col = quant_col2, lwd = quant_lwd, ...)
                  lines(quantiles[[x]] ~ time, lty = 1, col = quant_col, lwd = quant_lwd)
                }
                for (i in 1:(length(quantiles) - 1)) {
                  lines(quantiles[[i]] ~ time, lty = 1, col = quant_col, lwd = quant_lwd)
                  lines(quantiles2[[i]] ~ time, lty = 3, col = quant_col2, lwd = quant_lwd)
                  lines(exdiv ~ time, col = exval_col, lwd = exval_lwd)
                  lines(exdiv2 ~ time, col = exval_col2, lwd = exval_lwd)
                }
            } else {
                plot(quantiles[[x]] ~ time, type = "l", ann = FALSE, col = quant_col, lwd = quant_lwd, ...)
                for (i in 1:(length(quantiles) - 1)) {
                  lines(quantiles[[i]] ~ time, lty = 1, col = quant_col, lwd = quant_lwd)
                  lines(exdiv ~ time, col = exval_col, lwd = exval_lwd)
                }
            }
        }
        if (quantile == FALSE) {
            if(model %in% c("DA_OU", "DA_BM", "OU_BM")) {
              if(max(exdiv, na.rm=TRUE) > max(exdiv2, na.rm=TRUE)) {
                plot(exdiv ~ time, type = "l", ann = FALSE, col = exval_col, lwd = exval_lwd, ...)
                lines(exdiv2 ~ time, col = exval_col2, lty=1, lwd = exval_lwd)
              } else {
                plot(exdiv2 ~ time, type = "l", ann = FALSE, col = exval_col2, lwd = exval_lwd, ...)
                lines(exdiv ~ time, col = exval_col, lty=1, lwd = exval_lwd)
              }          
            } else {
              plot(exdiv ~ time, type = "l", ann = FALSE, col = exval_col, lwd = exval_lwd, ...)
            }
        }
        if (labels == TRUE) {
            title(xlab = "Age of Lineage Pair")
            title(ylab = "Trait Divergence")
        }
    }
    return(res)
}