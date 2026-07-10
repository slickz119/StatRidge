#' Logistic Ridge Regression
#'
#' Fits a logistic ridge regression model using maximum likelihood
#' estimation with L2 (ridge) regularization.
#'
#' @param X Numeric predictor matrix.
#' @param y Binary response vector (0/1).
#' @param lambda Ridge penalty parameter.
#'
#' @return A list containing the estimated regression coefficients and
#' optimization details.
#'
#' @examples
#' # Example (to be added in future releases)
#'
#' @export

logisticRidge <- function(X, y, lambda){

    loss_function <- function(beta){

        p <- 1/(1+exp(-X %*% beta))

        p <- pmin(pmax(p,1e-10),1-1e-10)

        negLL <- -sum(
            y*log(p) +
            (1-y)*log(1-p)
        )

        penalty <- lambda*sum(beta[-1]^2)

        negLL + penalty
    }

    beta0 <- rep(0,ncol(X))

    fit <- optim(
        par=beta0,
        fn=loss_function,
        method="BFGS",
        control=list(maxit=1000)
    )

    result <- list(

        coefficients=fit$par,

        lambda=lambda,

        convergence=fit$convergence,

        objective=fit$value,

        method="BFGS"

    )

    class(result) <- "StatRidge"

    return(result)

}
