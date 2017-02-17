


bc <- read.table("bc.csv",sep=",",header=T)

head(bc)
str(bc)

n.records <- nrow(bc)


##
## (a)
##

y <- rep(NA,n.records)

y[bc$Diagnosis=="M"] <- 1
y[bc$Diagnosis=="B"] <- 0
str(y)

y
X <- bc[,-c(1,2)]
X <- as.matrix(X)
str(X)


##
## (b)
##


nll <- function(mu,y) {

    n <- length(y)
    ## calculate mu
    mu.vec <- rep(mu,n)
    ## calculate p
    p.vec <- 1/(1+exp(-mu.vec))
    ## calculate negative log-likelihood
    loglik.vals <- y*log(p.vec)+(1-y)*log(1-p.vec)
    nll <- -sum(loglik.vals)
    ## return value
    return(nll)
}

nll(2,y)


##
## (c)
##


##
## (d)
##


fit <- optim(0, nll, y=y, hessian=TRUE)
fit

mu.hat <- fit$par
mu.se <- sqrt(1/fit$hessian)
CI <- c(mu.hat - 1.96*mu.se, mu.hat + 1.96*mu.se)

print(CI)


##
## (e)
##

nll.beta <- function(beta,y,X) {

    ## note that "beta" should have 11 values in it
    beta.0 <- beta[1]
    beta.reg <- beta[-1]
    ## calculate mu
    mu.vec <- beta.0+X%*%beta.reg
    ## calculate p
    p.vec <- 1/(1+exp(-mu.vec))
    ## calculate negative log-likelihood
    loglik.vals <- y*log(p.vec)+(1-y)*log(1-p.vec)
    nll <- -sum(loglik.vals)
    ## return value
    return(nll)
}

beta.start <- rep(0, 11)

nll.beta(beta.start, y, X)

fit <- optim(beta.start, nll.beta, y=y, X=X, hessian=TRUE)
fit

beta.symmetry <- fit$par[10]
se.beta <- sqrt(solve(fit$hessian)[10,10])
CI <- c(beta.symmetry - se.beta, beta.symmetry + se.beta)
print(CI)


