#' Poisson Ridge Regression
#'
#' Fits a Poisson ridge regression model using gradient descent with
#' L2 (ridge) regularization.
#'
#' @param X Predictor matrix.
#' @param y Response vector.
#' @param lambda Ridge parameter.
#' @param alpha Learning rate.
#' @param iterations Number of iterations.
#'
#' @return A fitted StatRidge model.
#'
#' @export

poissonRidge <- function(
    X,
    y,
    lambda,
    alpha = 0.01,
    iterations = 1000
){

    gradient <- function(beta){

        mu <- exp(X %*% beta)

        mu <- pmin(mu,1e10)

        t(X) %*% (y-mu) -
            2*lambda*beta

    }

    beta <- matrix(
        0,
        ncol(X),
        1
    )

    for(i in seq_len(iterations)){

        beta <- beta +
            alpha*gradient(beta)

    }

    result <- list(

        coefficients=beta,

        lambda=lambda,

        iterations=iterations,

        method="Gradient Descent"

    )

    class(result) <- "StatRidge"

    return(result)

}
