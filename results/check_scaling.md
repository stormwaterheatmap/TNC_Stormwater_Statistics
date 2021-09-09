TSS predictor transformation
================
September 8, 2021

This notebook compares scaling and centering of predictors between TSS
regression results and earth engine data layers.

``` r
library(dplyr)
load(here::here("results","TSS_models.RData"))
```

## Regression results

``` r
results <- TSS_models$Landscape_Predictor_Model
coeffs <- results$coefficients$fixed
TSS.data <- results$data
params <- coeffs %>% names()

traffic_data <- TSS.data$traffic %>% as.numeric()
paved_data <- TSS.data$paved


traffic_stats <- data.frame(
  min = traffic_data %>% min(),
 max= traffic_data %>% max(),
 mean= traffic_data %>% mean(), 
 sd = traffic_data %>% sd()
 )


traffic_stats
```

    ##          min      max       mean        sd
    ## 1 -0.9767572 1.592913 0.03113167 0.9707054

## Earth engine results

Code is hosted here:
<https://code.earthengine.google.com/c35bb989a54478723a0c5d006f72aec5>

EE output shown below

      "max": 1.4611211927889414,
      "mean": -6.050715484207103e-15,
      "min": -0.983404831055255,
      "sample_sd": 0.9999999999999934,

These look different - I think that the scaling and centering in the
regression model was done before removal of Pierce County watersheds?
