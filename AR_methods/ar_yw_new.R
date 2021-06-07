set.seed(10)

sourceCpp('eureka.cpp')

# Correct for predicted order = 0
ar.yw.threshold <-
  function (x, aic = TRUE, order.max = NULL, na.action = na.fail,
            demean = TRUE, series = NULL, ...)
  {
    if(is.null(series)) series <- deparse1(substitute(x))
    ists <- is.ts(x)
    x <- na.action(as.ts(x))
    if(ists) xtsp <- tsp(x)
    xfreq <- frequency(x)
    x <- as.matrix(x)
    if(!is.numeric(x))
      stop("'x' must be numeric")
    if(any(is.na(x) != is.na(x[,1]))) stop("NAs in 'x' must be the same row-wise")
    nser <- ncol(x)
    if (demean) {
      xm <- colMeans(x, na.rm=TRUE)
      x <- sweep(x, 2L, xm, check.margin=FALSE)
    } else xm <- rep.int(0, nser)
    n.used <- nrow(x)
    n.obs <- sum(!is.na(x[,1])) # number of non-missing rows
    order.max <- if (is.null(order.max))
      min(n.obs - 1L, floor(10 * log10(n.obs))) else floor(order.max)
    if (order.max < 1L) stop("'order.max' must be >= 1")
    else if (order.max >= n.obs) stop("'order.max' must be < 'n.obs'")
    xacf <- acf(x, type = "covariance", lag.max = order.max, plot = FALSE,
                demean=demean, na.action = na.pass)$acf
    if (xacf[1L] == 0) stop("zero-variance series")
    r <- as.double(drop(xacf))
    #print("yoohoo")
    z <- eureka(as.integer(order.max), r, r, matrix(double(order.max^2), nrow = order.max), 
                double(order.max), double(order.max), 0.001)
    # z <- .Fortran(C_eureka, as.integer(order.max), r, r,
    #               coefs = double(order.max^2),
    #               vars = double(order.max),
    #               double(order.max))
    order = z$order
    ar = if (order) z$coefs[order, seq_len(order)] else numeric()
    var.pred <- c(r[1L], z$vars[1:order.max])
    var.pred <- var.pred[order + 1L]
    var.pred <- var.pred * n.obs/(n.obs - (order + 1L))
    # if(z$order) {
    #   coefs <- matrix(z$coefs[1:order.max, 1:order.max], nrow = order.max)
    #   #coefs <- matrix(z$coefs, order.max, order.max)
    #   partialacf <- array(diag(coefs), dim = c(order.max, 1L, 1L))
    #   var.pred <- c(r[1L], z$vars[1:order.max])
    #   #var.pred <- c(r[1L], z$vars)
    #   xaic <- n.obs * log(var.pred) + 2 * (0L:order.max) + 2 * demean
    #   maic <- min(aic)
    #   xaic <- setNames(if(is.finite(maic)) xaic - min(xaic) else
    #     ifelse(xaic == maic, 0, Inf),
    #     0L:order.max)
    #   order <- if (aic) (0L:order.max)[xaic == 0L] else order.max
    #   ar <- if (order) coefs[order, seq_len(order)] else numeric()
    #   var.pred <- var.pred[order + 1L]
    #   ## Splus compatibility fix
    #   var.pred <- var.pred * n.obs/(n.obs - (order + 1L))
    # }
    # else  {
    #   order = 0
    #   ar = 0
    #   var.pred = 0
    #   xaic = 0
    #   partialacf = 0
    # }
    
    # resid <- if(order) c(rep.int(NA, order), embed(x, order + 1L) %*% c(1, -ar))
    # else as.vector(x) # we had as.matrix() above
    # if(ists) {
    #   attr(resid, "tsp") <- xtsp
    #   attr(resid, "class") <- "ts"
    # }
    res <- list(order = order, ar = ar, var.pred = var.pred, x.mean  =  drop(xm),
                n.used = n.used, n.obs = n.obs, order.max = order.max, method = "Yule-Walker",
                series = series, frequency = xfreq, call = match.call())
    # res <- list(order = order, ar = ar, var.pred = var.pred, x.mean  =  drop(xm),
    #             aic  =  xaic, n.used = n.used, n.obs = n.obs, order.max = order.max,
    #             partialacf = partialacf, resid = resid, method = "Yule-Walker",
    #             series = series, frequency = xfreq, call = match.call())
    if(nser == 1L && order)
      res$asy.var.coef <- var.pred/n.obs *
      solve(toeplitz(drop(xacf)[seq_len(order)]))
    class(res) <- "ar"
    res
  }
