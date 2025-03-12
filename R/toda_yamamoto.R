#' Applies the Toda-Yamamoto Causality Test
#'
#' This function takes a list of data frames with time series and applies the Toda-Yamamoto test
#' to assess causality between variable pairs in VAR models. The Toda-Yamamoto test allows testing
#' for Granger causality even in the presence of non-stationary series by estimating a VAR model
#' with an additional lag. This approach avoids the need for pre-tests of stationarity,
#' preserving the validity of statistical tests.
#'
#' @param data_list A list of data frames containing time series. In each data frame, the first column
#'   should represent the "effect" variable, while the remaining columns are the potential "cause" variables
#'   to be tested.
#' @param criterion A character string indicating the criterion for selecting the optimal lag length 
#'   (e.g., "AIC", "BIC", etc.).
#'
#' @return A list of data frames, where each data frame contains the causality test results for the
#'   corresponding data frame in \code{data_list}. Each result data frame includes the following columns:
#'   \describe{
#'     \item{cause}{The variable considered as the cause.}
#'     \item{effect}{The variable considered as the effect.}
#'     \item{chi_square}{The chi-square statistic obtained from the Wald test.}
#'     \item{p_value}{The p-value associated with the Wald test.}
#'   }
#'
#' @details The Toda-Yamamoto test consists of estimating a VAR model with \code{p + 1} lags, where \code{p}
#'   is the optimal lag length selected by the specified criterion. The inclusion of the additional lag 
#'   allows applying the Granger causality test without requiring the series to be transformed into stationary ones,
#'   avoiding problems arising from stationarity pre-tests and loss of information.
#'
#' @examples
#' \dontrun{
#'   # Suppose df1 and df2 are data frames containing time series
#'   results <- toda_yamamoto(list(df1 = df1, df2 = df2), criterion = "AIC")
#' }
#'
#' @importFrom vars VARselect VAR
#' @importFrom aod wald.test
#' @importFrom dplyr bind_rows tibble
#' @importFrom purrr map
#' @export
toda_yamamoto <- function(data_list, criterion) {
  require(tidyverse)
  require(vars)
  require(aod)
  
  results <- map(data_list, function(df) {
    result_list <- list()
    
    for (j in 2:(ncol(df))) {
      # Select optimal lag length
      var_selection <- df[, c(1, j)] |> 
        VARselect(lag.max = 15, type = "none")
      
      selected_lag <- tryCatch(
        as.numeric(var_selection$selection[str_detect(
          names(var_selection$selection), toupper(criterion)
        )]),
        error = function(e) NA
      )
      
      if (is.na(selected_lag)) next  # Skip if criterion is not found
      
      # Fit the VAR model with an extra lag (p + 1)
      var_model <- df[, c(1, j)] |> VAR(p = selected_lag + 1, type = "none")
      
      # Wald tests to check causality
      test_forward <- wald.test(Sigma = vcov(var_model$varresult[[1]]),
                                b = coef(var_model$varresult[[1]]),
                                Terms = seq(2, (var_model$p * 2), by = 2))
      
      test_reverse <- wald.test(Sigma = vcov(var_model$varresult[[2]]),
                                b = coef(var_model$varresult[[2]]),
                                Terms = seq(1, (var_model$p * 2), 2))
      
      # Store results
      result_list <- append(result_list, list(
        tibble(
          cause = names(var_model$varresult)[2],
          effect = names(var_model$varresult)[1],
          chi_square = round(test_forward$result$chi2["chi2"], 3),
          p_value = round(test_forward$result$chi2["P"], 3)
        ),
        tibble(
          cause = names(var_model$varresult)[1],
          effect = names(var_model$varresult)[2],
          chi_square = round(test_reverse$result$chi2["chi2"], 3),
          p_value = round(test_reverse$result$chi2["P"], 3)
        )
      ))
    }
    
    return(bind_rows(result_list)) # Combine results into a data frame
  })
  
  return(set_names(results, names(data_list)))
}
