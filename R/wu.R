
#' Fit spike and slab Cox PHMs
#'
#' @param Y survival times
#' @param delta censoring indicator, 0: uncensored, 1: censored
#' @param X design matrix
#' @returns list with coefficient values and number of non-zero coefficients
#' 
#' @examples
#' # Generate data
#' set.seed(1)
#' n <- 250; omega <- 1; censoring_lvl <- 0.4
#' 
#' b <- c(1, 1, rep(0, 250))
#' p <- length(b)
#' X <- matrix(rnorm(n * p), nrow=n)
#' y <- runif(nrow(X))
#' Y <- log(1 - y) / - (exp(X %*% b) * omega)
#' 
#' delta  <- runif(n) > censoring_lvl   # 0: censored, 1: uncensored
#' Y[!delta] <- Y[!delta] * runif(sum(!delta))
#' 
#' fit <- wu.fit(Y, delta, X)
#'
#' @export
wu.fit <- function(Y, delta, X) {
    y <- survival::Surv(as.matrix(Y), as.matrix(as.integer(delta)))

    # find CV lambda
    lambda <- glmnet::cv.glmnet(X, y, alpha=0, standardize=FALSE,
	    family="cox")$lambda.1se
    
    # fit l2-penalised cox phm
    fit <- glmnet::glmnet(X, y, alpha=0, standardize=FALSE,family="cox", 
	    lambda=lambda)

    beta0 <- as.matrix(unlist(coef(fit)))

    aL <- lik_linear(Y, delta, X, beta0)
    linear_onto_fisher_col <- X[aL[[4]], ] %*% beta0 + 
	solve(diag(aL[[3]][aL[[4]], ])) %*% aL[[2]][aL[[4]], ]

    w.star0 <- diag(aL[[3]][aL[[4]], ])^0.5 
    y.star0 <- w.star0 %*% linear_onto_fisher_col
    x.star0 <- w.star0 %*% X[aL[[4]],]

    # variable selection via spikeslab
    x <- spikeslab::spikeslab(x=x.star0, y=y.star0, intercept=FALSE,
	    n.iter1=100, n.iter2=2500, center=FALSE, bigp.smalln=ncol(X) > nrow(X),
	    screen=FALSE)

    return(list(x.gnet=x$gnet, phat=x$phat))
}


#' @keywords internal
risk_set <- function(Y, delta)
{
    unique.Y <- unique(Y[delta == 1])
    order.Y <- sort(unique.Y)

    l <- length(unique.Y)
    n <- length(Y)
    risk <- rep(0,n) %*% t(rep(0, l))

    for(i in 1:l){
	risk[, i] <- (Y >= order.Y[i])
    }
    return(risk)
}


#' @keywords internal
lik_linear <- function(Y, delta, X , b)
{
    n <- nrow(X); p <- ncol(X); b <- as.matrix(b)

    unique.Y <- unique(Y[delta == 1])
    event.t <- sort(unique.Y)
    l <- length(event.t)

    risk <- risk_set(Y, delta)

    sub.hess <- array(0, c(n, 1, l))
    sub.grad <- array(0, c(n, 1, l))
    sub.devi <- rep(0,l)

    linear <- as.vector(X %*% b) 
    eta <- as.vector(exp(X %*% b))
    d.hazard <- diag(eta)
    sum.d <- 1/apply(d.hazard %*% risk, 2, sum)

    for(i in 1:l){
	u.stat <- eta * (Y >= event.t[i])*sum.d[i]   
	tie.event.set <- (Y == event.t[i] & delta == 1)

	if (sum(tie.event.set) != 1) {
	    tie=rep(1,sum(tie.event.set))       
	    sub.grad[,,i] <- as.matrix(tie.event.set) - sum(tie.event.set)*u.stat                     
	    sub.hess[,,i] <- sum(tie.event.set)*(u.stat - u.stat^2)
	    sub.devi[i] <- as.vector(linear[tie.event.set]) %*% tie - 
		sum(tie.event.set) * log(1/sum.d[i])               
	}
	if (sum(tie.event.set) == 1) {
	    sub.grad[,,i] <- as.matrix( tie.event.set ) - u.stat
	    sub.hess[,,i] <- u.stat - u.stat^2
	    sub.devi[i] <- as.numeric(linear[tie.event.set])-log(1/sum.d[i])
	}      
    }

    devi <- -2 * sum(sub.devi)
    score <- apply(sub.grad, c(1,2), sum)
    information <- apply(sub.hess, c(1,2), sum)
    approx.information <- round(information, 5)
    active <- (approx.information != 0)

    return(list(deviance = devi, score = score, 
       information = information, active.set = active))
}

