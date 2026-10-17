# This script takes model validation predictors and plots them atop best fit lines from S8 data to 
#  facilitate selection of the best model for each COC.


# Author: Eva Dusek Jennings
# Revised: Sep 25, 2026
#--------------------------------------------------------------------------------------------------

#install.packages("patchwork")
library(patchwork)





### first run up through line 294 from Model Validation v4.R


Cu.mv <- mv.full_res_pred %>%
  filter(Analyte=="Total Copper") %>%
  mutate(rain=daymet_14day_std,
         summer=as.factor(summer),
         rawResult = Result)

TSS.mv <- mv.full_res_pred %>%
  filter(Analyte=="Total Suspended Solids") %>%
  mutate(rain=daymet_precip_std,
         spring=as.factor(spring),
         rawResult = Result)

P.mv <- mv.full_res_pred %>%
  filter(Analyte=="Total Phosphorus") %>%
  mutate(rain=daymet_14day_std,
         rain1=daymet_precip_std,
         summer = as.factor(summer),
         rawResult = Result)

Zn.mv <- mv.full_res_pred %>%
  filter(Analyte=="Total Zinc") %>%
  mutate(rain=daymet_14day_std,
         summer = as.factor(summer),
         rawResult = Result)

TKN.mv <- mv.full_res_pred %>%
  filter(Analyte=="Total Kjeldahl Nitrogen") %>%
  mutate(rain=daymet_14day_std,
         rain1=daymet_precip_std,
         summer = as.factor(summer),
         rawResult = Result)



#-------------#
#  Functions  #
#-------------#
# 
# theModel <- Model
# myCOC <- "Cu"
# myModel <- Model
# aa <- aa[1]
# 

#function to plot best fit LME model line to a single predictor
plot.single.preds.MV <- function(theModel, myCOC) {
  aa <- summary(theModel)$coefficients$fixed
  aa <- names(aa[names(aa) %in% thePreds])
  ppp1 <- Plot.One.Predictor.wAgency.MV(aa[1], theModel, myCOC)
  if(length(aa) > 1)  ppp2 <- Plot.One.Predictor.wAgency.MV(aa[2], theModel, myCOC, hideLegend=FALSE) 
  #arrange plots on the grid
  if(length(aa)==1)  grid.arrange(ppp1, nrow=1, ncol=2)
  # if(length(aa)==2)  grid.arrange(ppp1, ppp2, nrow=1, ncol=2)
  if(length(aa)==2) ppp1 + ppp2 + plot_layout(guides="auto")
}


#Function to plot result vs one predictor, with location shown in color
Plot.One.Predictor.wAgency.MV <- function(aa, myModel, myCOC, hideLegend=TRUE, xEqn=0.15, yEqn=max(s8$result*1.02, mv$result*1.02) ) {
  coef.int <- myModel$coefficients[[1]]["(Intercept)"]
  coef.slope <- myModel$coefficients[[1]][aa]
  
  #linear equation for the primary predictor, to be printed in the plot
  eqnText <- paste("ln(", myCOC, ") = ", 
                   round(coef.int, 2), " + ",
                   round(coef.slope, 2), "*", aa )
  #plot the primary predictor on x-axis, with secondary predictor shown in color
  ggplot() + 
    geom_point(data=s8, aes(x=s8[,aa], y=result), color="darkgray", size=3) + 
#    scale_fill_manual(values = "darkgray") +
    #    geom_point(data=s8, aes(x=s8[,aa], y=result, fill=s8[,"agency"]), size=3, shape=21, stroke=0) + 
    #    scale_fill_manual(values = colors_agency) +
    geom_point(data=mv, aes(x=mv[,aa], y=result, color=Project), size=3) +
#    scale_color_manual(values = colors_agency) +
#    scale_fill_viridis() +
    geom_abline(intercept=coef.int, slope=coef.slope, col="black", lwd=1.5) +  #, aes(x = aa, y = y.hat), col="darkred", lwd = 1.5)
    ylab(paste("ln", myCOC)) + xlab(aa) +
    geom_text(aes(x=xEqn, y=yEqn, label=eqnText), cex=4.5, color="black")+  #add equation to the top of the plot
    #    labs(color="agency", fill="agency") +
    theme_gray() +
    if (hideLegend) {
      theme(legend.position="none")
    }
}    




#------------#
#  All COCs  #
#------------#

r1X <- formula( ~ 1 | agency/location)  #random structure for best fit model
vf1X <- varIdent(form= ~1|location)   #variance structure for best fit model

thePreds <- c("totRES", "roof_totRES", "sqrt_popn", "ln_traffic_k50", "ln_traffic_k100", "oldDev", "devAge", "trees", "impervious", "paved",
              "sqrt_CO2_road", "sqrt_CO2_res", "sqrt_CO2_transport")

colors_agency <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7")
colors_location <- colors_agency[c(1,1,1,2,2,2,3,4,4,4,5,5,5,6,6,6)]


#----------#
#  Copper  #
#----------#

# setwd("/Users/edj/Documents/Consulting/TNC/Statistical Analysis")
load("results/S8 Data & Model Structure_Copper.RData")

s8 <- Cu.coc2
mv <- Cu.mv
my.COC <- "Cu"

# Form4a <- formula(result ~ rain + summer + dry + totRES + ln_traffic_k50)  
# Form4b <- formula(result ~ rain + summer + dry + sqrt_CO2_road + oldDev)
# Form4c <- formula(result ~ rain + summer + dry + impervious + ln_traffic_k100)

Form <- formula(result ~ rain + summer + totRES + ln_traffic_k50)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + sqrt_CO2_road + oldDev)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + impervious + ln_traffic_k100)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)


#-------#
#  TSS  #
#-------#

load("results/S8 Data & Model Structure_TSS.RData")

s8 <- TSS.coc2
mv <- TSS.mv
my.COC <- "TSS"

# Form4a <- formula(result ~ rain + spring + rain21 + oldDev + sqrt_CO2_road)  
# Form4b <- formula(result ~ rain + spring + rain21 + oldDev)

Form <- formula(result ~ rain + spring + oldDev + sqrt_CO2_road)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + spring + oldDev)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)


#--------------#
#  Phosphorus  #
#--------------#

load("results/S8 Data & Model Structure_Phosphorus.RData")

s8 <- P.coc2
mv <- P.mv
my.COC <- "P"

# Form4a <- formula(result ~ rain + summer + rain1 + trees)  
# Form4b <- formula(result ~ rain + summer + rain1 + oldDev)
# Form4c <- formula(result ~ rain + summer + rain1 + impervious)

Form <- formula(result ~ rain + rain1 + summer + trees)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + rain1 + summer + oldDev)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + rain1 + summer + impervious)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + rain1 + summer + devAge)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)


Form <- formula(result ~ rain + rain1 + summer + sqrt_CO2_road + oldDev)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)


#--------#
#  Zinc  #
#--------#

load("results/S8 Data & Model Structure_Zinc.RData")

s8 <- Zn.coc2
mv <- Zn.mv
my.COC <- "Zn"

# Form4a <- formula(result ~ rain + summer + trees + sqrt_CO2_transport)
# Form4b <- formula(result ~ rain + summer + impervious + sqrt_CO2_transport)  
# Form4c <- formula(result ~ rain + summer + totRES + sqrt_CO2_transport)
# Form4d <- formula(result ~ rain + summer + devAge + roof_totRES)
# Form4e <- formula(result ~ rain + summer + oldDev + sqrt_CO2_transport)
# Form4f <- formula(result ~ rain + summer + totRES)

Form <- formula(result ~ rain + summer + trees + sqrt_CO2_transport)    #we don't have sqrt_CO2_total, so here is transport to make the point of paved being not good
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + impervious + sqrt_CO2_transport)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + totRES + sqrt_CO2_transport)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + devAge + roof_totRES)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + oldDev + sqrt_CO2_transport)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + totRES)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)



#-------#
#  TKN  #
#-------#

load("results/S8 Data & Model Structure_TKN.RData")

s8 <- TKN.coc2
mv <- TKN.mv
my.COC <- "TKN"

# Form4a <- formula(result ~ rain + summer + oldDev + sqrt_CO2_road)
# Form4b <- formula(result ~ rain + summer + oldDev + sqrt_CO2_res)  
# Form4c <- formula(result ~ rain + summer + oldDev + sqrt_popn)  
# Form4d <- formula(result ~ rain + summer + sqrt_popn + ln_traffic_k50)  
# Form4e <- formula(result ~ rain + summer + sqrt_CO2_road)  
# Form4f <- formula(result ~ rain + summer + paved + sqrt_CO2_res)

Form <- formula(result ~ rain + summer + sqrt_CO2_road + oldDev)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + sqrt_CO2_res + oldDev)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + sqrt_popn + oldDev)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + ln_traffic_k50 + sqrt_popn)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + sqrt_CO2_road)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)

Form <- formula(result ~ rain + summer + paved + sqrt_CO2_res)  
Model <- lme(data=s8, Form, method="ML", random=r1X, weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
plot.single.preds.MV(Model, my.COC)


