expected_val =
function (sig2, alpha = NULL, psi = NULL, psi2 = NULL, wt = NULL, 
    time_span = c(0, 10), quantile = FALSE, plot = FALSE, 
    labels = TRUE, exval_col = "red", quant_col = "black", 
    exval_lwd = 1, quant_lwd = 1, ...) 
{
    if (length(time_span) == 1) {
        time <- time_span
    }
    if (length(time_span) > 2) {
        time <- time_span
    }
    if (length(time_span) == 2) {
        time = seq(0, time_span[2], (time_span[2]/10000))
    }
    if (is.null(alpha)) {
        V <- sig2 * time * 2
        u = rep(0, length(time))
        exdiv <- ((2 * V)/pi)^0.5 # expected value of absolute divergence
    }
    if (is.null(alpha) == FALSE) {
        V <- (sig2/(alpha)) * (1 - exp(-2 * alpha * time))
        if (is.null(psi)) {
            u = rep(0, length(time))
            exdiv <- ((2 * V)/pi)^0.5
        }
        if (is.null(psi) == FALSE & is.null(psi2)) {
            u <- psi * (1 - exp(-alpha * time))
            exdiv <- sqrt((2 * V)/pi) * exp(-(u^2)/(2 * V)) + 
               u * (2*pnorm((u/sqrt(2 * V))*sqrt(2)) - 1)
        }
        if (is.null(psi) == FALSE & is.null(psi2) == FALSE) {
            if (is.null(wt)) 
                stop("You need to enter a time (wt) at which the shift in epoch occurs")
            if (wt > max(time)) 
                stop("You've chosen an epoch shift that occurs after the max time in this calculation")
            TIME1 = time[0:which(time == wt)]
            TIME2 = time[(which(time == wt) + 1):length(time)]
            u1 = psi * (1 - exp(-alpha * TIME1))
            u2 = exp(-alpha * TIME2) * (psi * (exp(alpha * wt) - 
                1) + psi2 * (exp(alpha * TIME2) - exp(alpha * 
                wt)))
            u = c(u1, u2)
            exdiv <- sqrt((2 * V)/pi) * exp(-(u^2)/(2 * V)) + 
                u * (2*pnorm((u/sqrt(2 * V))*sqrt(2)) - 1)
        }
    }
    if (quantile == TRUE) {
        qs = c(seq(0.1, 0.9, 0.1), 0.95)
        quantiles = lapply(qs, FUN = truncnorm::qtruncnorm, a = 0, 
            mean = u, sd = sqrt(V))
        res <- matrix(NA, length(time), 12)
        colnames(res) <- c("time", "Expectation", "q10", "q20", 
            "q30", "q40", "q50", "q60", "q70", "q80", "q90", 
            "q95")
        res[, 1] <- time
        res[, 2] <- exdiv
        for (i in 3:12) res[, i] = quantiles[[i - 2]]
    }
    if (quantile == FALSE) {
        res <- matrix(NA, length(time), 2)
        colnames(res) <- c("time", "Expected_Abs_Diff")
        res[, 1] <- time
        res[, 2] <- exdiv
    }
    if (plot == TRUE) {
        if (quantile == TRUE) {
        	x = length(quantiles)
        	plot(quantiles[[x]] ~ time, type = "l", ann=FALSE, col=quant_col, lwd = quant_lwd, ...)
            for (i in 1:(length(quantiles)-1)) {
                lines(quantiles[[i]] ~ time, lty = 1, col = quant_col, lwd = quant_lwd)
                lines(exdiv ~ time, col=exval_col, lwd = exval_lwd)
            }
        }
        if (quantile == FALSE) {
            plot(exdiv ~ time, type = "l", ann=FALSE, col=exval_col, lwd = exval_lwd, ...)
        }
        if (labels == TRUE) {
            title(xlab = "Age of Lineage Pair")
            title(ylab = "Trait Divergence")
        }
    }
    return(res)
} 