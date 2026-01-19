
<!-- README.md is generated from README.Rmd. Please edit that file -->

# iProPath

<!-- badges: start -->
<!-- badges: end -->

## Installation

``` r
install.packages("iProPath_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Example

### 1. Set basic parameters

``` r
library(iProPath)
file_name <- "melanoma_final_v1"
num_cores <- 6
y_name <- "trans_glyco"
x_name <- "mutation"
```

### 2. Load input data

``` r
load(file.path("Data", paste0(file_name, ".RData")))
```

Loads:

x_df_list, fixed_cov, mediator_list, y_covdf_list, Y_list

### 3. Regression step

``` r
reg_res <- prepare_ipropath_regression(
  x_df_list,
  fixed_cov,
  mediator_list,
  y_covdf_list,
  Y_list,
  x_name = x_name,
  y_name = y_name,
  num_cores = num_cores,
  adjust_methods = c("BH", "bonferroni", "BY")
)
```

### 4. Estimate non-null proportions

``` r
pi_res <- estimate_ipropath_pi(
reg_res$merged_df$p_value_MX,
reg_res$merged_df$p_value_YM,
method = "qvalue.smoother",
run_all = FALSE
)
```

### 5. Run iProPath inference

``` r
ipropath_res <- run_ipropath(
  reg_res$merged_df,
  pi_res$pi1_MX,
  pi_res$pi1_YM,
  num_cores = num_cores,
  y_name = y_name
)
```

### 6. Assemble final results

``` r
res_all <- assemble_ipropath_results(
  pairwise_res = reg_res$pairwise_res,
  y_transreg = reg_res$y_transreg,
  ipropath_res_p = ipropath_res$res_p
)
```

### 7. Write outputs

``` r
file_stub <- paste(file_name, x_name, y_name, sep = ".")
write_ipropath_outputs(
  res_all = res_all,
  real_data_pi_df = pi_res$pi_df,
  ipropath_pi_df = ipropath_res$pi_df,
  file_stub = file_stub,
  output_dir = "Results"
)
```
