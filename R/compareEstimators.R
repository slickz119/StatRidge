#' Compare Ridge Regression Estimators
#'
#' Fits logistic ridge regression models using multiple ridge
#' parameter estimators.
#'
#' @param X Predictor matrix.
#' @param y Binary response vector.
#' @param ridge A list of ridge parameter estimators.
#'
#' @return A list of fitted ridge regression models.
#'
#' @export

compareEstimators <- function(X, y, ridge){

    fits <- list()

    for(name in names(ridge)){

        if(!is.na(ridge[[name]])){

            fits[[name]] <-
                logisticRidge(
                    X,
                    y,
                    lambda = ridge[[name]]
                )

        }

    }

    return(fits)

}
