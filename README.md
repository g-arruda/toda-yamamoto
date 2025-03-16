# Toda-Yamamoto Causality Test

## Overview

This code provides an implementation of the Toda-Yamamoto causality test for time series analysis in R. The Toda-Yamamoto procedure is an extension of the Granger causality test that allows for testing causality between variables even when the series are non-stationary or have different orders of integration.

## Why Toda-Yamamoto?

Traditional Granger causality tests require stationary time series data, which often means differencing or transforming data. This can lead to:
- Loss of information contained in the original levels
- Pre-testing bias from unit root and cointegration tests
- Incorrect model specification

The Toda-Yamamoto approach addresses these issues by allowing for the testing of Granger causality in level VAR models regardless of the integration or cointegration properties of the data.


## Dependencies

This package requires:
- tidyverse
- vars
- aod
- purrr
- dplyr

## Usage

```r
# Example with two data frames containing time series
results <- toda_yamamoto(
  list(df1 = df1, df2 = df2),
  criterion = "AIC"
)
```

### Input Format

- `data_list`: A list of data frames where each data frame contains your time series variables
  - The first column should be the "effect" variable
  - Remaining columns should be potential "cause" variables
- `criterion`: Selection criterion for determining optimal lag length (e.g., "AIC", "BIC", "HQ")

### Output

The function returns a list of data frames (one for each input data frame) with the following columns:
- `cause`: The variable being tested as the causal factor
- `effect`: The outcome variable
- `chi_square`: Chi-square statistic from the Wald test
- `p_value`: P-value for the causality hypothesis test

## Method Details

The Toda-Yamamoto procedure works as follows:

1. Determine the maximum order of integration (d) for the variables in the system
2. Select the optimal lag length (p) for the VAR model using information criteria
3. Estimate a VAR(p+d) model (add d extra lags beyond the optimal lag length)
4. Perform Wald tests on the first p lags only, ignoring the extra d lags

This implementation automatically selects the optimal lag length using the specified criterion and adds one additional lag, which is suitable for most economic and financial time series (which are typically I(1)).

## Example

```r
# Create example time series data
set.seed(123)
n <- 200
x <- cumsum(rnorm(n))
y <- cumsum(0.7*x + rnorm(n))
z <- cumsum(0.3*y + rnorm(n))

df1 <- data.frame(y, x)
df2 <- data.frame(z, x, y)

# Apply Toda-Yamamoto test
results <- toda_yamamoto(
  list(df1 = df1, df2 = df2),
  criterion = "AIC"
)

# View results
results$df1
results$df2
```

## References

```
Toda, H. Y., & Yamamoto, T. (1995). Statistical inference in vector autoregressions with possibly integrated processes. Journal of Econometrics, 66(1-2), 225-250.
```


## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
