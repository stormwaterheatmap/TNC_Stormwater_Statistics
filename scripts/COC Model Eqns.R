# This script assembles the model equations for each COC.

# Eva Dusek Jennings
# April 27, 2021 
#--------------------------------------------------------------

#-----------------#
#  Copper Models  #
#-----------------#

source("Copper Model_v5.R")
Cu.coc2 <- coc2
Cu.r1X <- r1X
Cu.vf1X <- vf1X
Cu.Form4 <- Form4
Cu.Form5 <- Form5

Cu.M1 <- gls(data=Cu.coc2, result~1, method="REML")  #here, we use the actual concentration, rather than the transformed one
Cu.M3 <- lme(data=Cu.coc2, result~landuse+rain, random = Cu.r1X, method="REML", weights=Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
Cu.M4 <- lme(data=Cu.coc2, Cu.Form4, random = Cu.r1X, method="REML", weights=Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
Cu.M5 <- lme(data=Cu.coc2, Cu.Form5, random = Cu.r1X, method="REML", weights=Cu.vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))



#---------------#
#  Zinc Models  #
#---------------#



#--------------#
#  TSS Models  #
#--------------#



#--------------------------#
#  Nitrite-Nitrate Models  #
#--------------------------#



