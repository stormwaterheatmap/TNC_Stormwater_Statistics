# Functions for Generalized COC Analysis and COC analysis scripts


#-----------------------------------#
#  Functions for Exploratory Plots  #
#-----------------------------------#

#obtain slope for grob function
slope_grob <- function(predictor) {
  my.lm <- lm(result~predictor, data=coc)
  grob <- grobTree(textGrob(paste("slope p-value =", round(summary(my.lm)$coef[2,4], 4)), x=0.1,  y=0.95, hjust=0,
                            gp=gpar(col="red", fontsize=15, fontface="italic")))
  return(grob)
}

my.ggplot <- function(pred_num) {
  ggplot(coc, aes(coc[,predictors[pred_num]], result)) + geom_point() + xlab(predictors[pred_num]) + 
    geom_smooth(method = "lm") + annotation_custom(slope_grob(coc[, predictors[pred_num]]))  #add a smooth slope and p-value for the slope line
}


#-------------------------------------#
#  Functions for Plotting Model Fits  #
#-------------------------------------#

# myModel <- M.gls1X
# myE <- E.1X
# myForm <- "X"

# #this function plots residuals from models, to examine for heterogeneity
# plot.resids <- function(myModel, myE, myForm) {
#   op <- par(mfrow=c(4,4), mar=c(4,4,1,1))
#   if (myForm=="A") {
#     mydf <- data.frame(roofs=coc2$roofs, grass=coc2$grass, trees=coc2$trees)
#   } else if (myForm=="B") {
#     mydf <- data.frame(paved=coc2$paved, roofs=coc2$roofs, grass=coc2$grass)
#   } else if (myForm=="C") {
#     mydf <- data.frame(impervious=coc2$impervious, grass=coc2$grass)
#   } else if (myForm=="X") {
#     mydf <- data.frame(paved=coc2$paved, roofs=coc2$roofs, grass=coc2$grass, trees=coc2$trees)
#   }
#   plot(fitted(myModel), myE, xlab="Fitted Values", ylab="Residuals", col="gray", pch=16)
#   plot(coc2$location, myE, xlab="location", ylab="Residuals", col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)])
#   plot(coc2$agency, myE, xlab="agency", ylab="Residuals", col=colors_agency)
#   plot(coc2$month, myE, xlab="month", ylab="Residuals", col=c(rep("light blue", 3), rep("light green", 3), rep("yellow", 3), rep("orange", 3)))
#   plot(coc2$season, myE, xlab="season", ylab="Residuals", col=c("light blue", "light green", "yellow", "orange"))
#   plot(coc2$landuse, myE, xlab="land use", ylab="Residuals", col=c("light blue", "light green", "yellow", "orange"))
#   plot(mydf[,1], myE, xlab=names(mydf)[1], ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   plot(mydf[,2], myE, xlab=names(mydf)[2], ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   if (ncol(mydf)==3) {
#     plot(mydf[,3], myE, xlab=names(mydf)[3], ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   }
#   if (ncol(mydf)==4) {
#     plot(mydf[,4], myE, xlab=names(mydf)[4], ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   }
#   plot(coc2$nodev, myE, xlab="nodev", ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   plot(coc2$pm25, myE, xlab="pm25", ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   plot(coc2$traffic, myE, xlab="traffic", ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   plot(coc2$slope, myE, xlab="slope", ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   plot(coc2$rain, myE, xlab="rain", ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   plot(coc2$mPrecip, myE, xlab="monthly precip", ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   plot(coc2$dry, myE, xlab="ant. dry days", ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
# }

# #this function generates a figure of boxplots showing model fit to all possible parameters
# boxplots.resids <- function(myModel, myE, myForm) {
#   aa <- coc2[, c("agency", "impervious", "paved", "roofs", "grass", "trees", "nodev", "pm25", "traffic", "slope")]
#   bb <- unique(aa)
#   op <- par(mfrow=c(4,5), mar=c(4,4,2,1))
#   plot(fitted(myModel), myE, xlab="Fitted Values", ylab="Residuals", col="gray", pch=16)
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$location, col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)], xlab="locations", ylab="Residuals")
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$landuse, col=colors_agency, xlab="land use", ylab="Residuals")
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$agency, col=colors_agency, xlab="agency", ylab="Residuals")
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$year, xlab="year", ylab="Residuals")
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$month, xlab="month", ylab="Residuals", col=c(rep("light blue", 3), rep("light green", 3), rep("yellow", 3), rep("orange", 3)))
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$season, xlab="season", ylab="Residuals", col=c("light blue", "light green", "yellow", "orange"))
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$impervious, xlab="impervious", ylab="Residuals", col=colors_agency[as.numeric(bb[order(bb[,"impervious"]), "agency"])])
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$paved, xlab="paved", ylab="Residuals", col=colors_agency[as.numeric(bb[order(bb[,"paved"]), "agency"])])
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$roofs, xlab="roofs", ylab="Residuals", col=colors_agency[as.numeric(bb[order(bb[,"roofs"]), "agency"])])
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$grass, xlab="grass", ylab="Residuals", col=colors_agency[as.numeric(bb[order(bb[,"grass"]), "agency"])])
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$trees, xlab="trees", ylab="Residuals", col=colors_agency[as.numeric(bb[order(bb[,"trees"]), "agency"])])
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$nodev, xlab="nodev", ylab="Residuals", col=colors_agency[as.numeric(bb[order(bb[,"nodev"]), "agency"])])
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$pm25, xlab="pm25", ylab="Residuals", col=colors_agency[as.numeric(bb[order(bb[,"pm25"]), "agency"])])
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$traffic, xlab="traffic", ylab="Residuals", col=colors_agency[as.numeric(bb[order(bb[,"traffic"]), "agency"])])
#   abline(0,0, col="gray")
#   boxplot(myE~coc2$slope, xlab="slope", ylab="Residuals", col=colors_agency[as.numeric(bb[order(bb[,"slope"]), "agency"])])
#   abline(0,0, col="gray")
#   plot(coc2$rain, myE, xlab="rain", ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
#   abline(0,0, col="gray")
#   #  plot(coc2$dry, myE, xlab="ant. dry days", ylab="Residuals", col=colors_agency[as.numeric(coc2$agency)])
# }
# 


#a more presentable version of the plotting function above
boxplots.resids2 <- function(myModel, myE, myForm) {
  aa <- coc2[, c("agency", all_of(best_predictors))]
  bb <- unique(aa)
  op <- par(mfrow=c(4,5), mar=c(0.5,1,2.5,1))
  # plot(fitted(myModel), myE, xlab="Fitted Values", ylab="Residuals", col=colors_agency[c(1,1,1,2,2,2,3,4,4,4,5,5,5,6,6,6)])
  # abline(0,0, col="gray")
  plot(fitted(myModel), myE, main="residuals", xlab="", ylab="", xaxt="n", yaxt="n", col="gray", pch=16)
  boxplot(myE~coc2$location, col=colors_agency[c(1,1,1,2,3,4,4,4,5,5,5,6,6,6)], main="locations", ylab="", xaxt="n", yaxt="n")
  abline(0,0, col="gray")
  boxplot(myE~coc2$landuse, col=c("green", "dark green", "pink", "gray"), main="land use", ylab="", xaxt="n", yaxt="n")
  abline(0,0, col="gray")
  boxplot(myE~coc2$agency, col=colors_agency, main="agency", ylab="", xaxt="n", yaxt="n")
  abline(0,0, col="gray")
  # boxplot(myE~coc2$year, main="year", ylab="", xaxt="n", yaxt="n")
  # abline(0,0, col="gray")
  boxplot(myE~coc2$month, main="month", ylab="", xaxt="n", yaxt="n", col=c(rep("light blue", 3), rep("light green", 3), rep("yellow", 3), rep("orange", 3)))
  abline(0,0, col="gray")
  boxplot(myE~coc2$season, main="season", ylab="", xaxt="n", yaxt="n", col=c("light blue", "light green", "yellow", "orange"))
  abline(0,0, col="gray")
  
  for (i in 1:length(best_predictors)) {
    boxplot(myE ~ dplyr::pull(coc2, best_predictors[i]), main=best_predictors[i], ylab="", xaxt="n", yaxt="n",
            col=colors_agency[as.numeric(bb[order(pull(bb, best_predictors[i])), "agency"])]) 
    abline(0,0, col="gray")
  }

  plot(coc2$rain, myE, main="rain", xlab="", ylab="", xaxt="n", yaxt="n", col="gray", pch=16)
  abline(0,0, col="gray")
  plot(coc2$dry, myE, main="ant. dry days", xlab="", ylab="", xaxt="n", yaxt="n", col="gray", pch=16)
  abline(0,0, col="gray")
}



#--------------------------------------------#
#  Functions for Plotting Model Predictions  #
#--------------------------------------------#

plot.preds.vs.results <- function(model.pred) {
  ylimits <- c(0.75*min(coc2$result), 1.05*max(coc2$result))
  
  par(mfrow=c(3,3), mar=c(2,2,2,1), oma=c(0,0,0,0))
  plot(model.pred ~ coc2$paved, col="pink", pch=19, cex=2, ylim=ylimits, xlab="", ylab="", main="paved", xaxt="n", yaxt="n", cex.main=2)
  points(coc2$result ~ coc2$paved, col="black", pch="*")
  
  plot(model.pred ~ coc2$trees, col="light green", pch=19, cex=2, ylim=ylimits, xlab="", ylab="", main="trees", xaxt="n", yaxt="n", cex.main=2)
  points(coc2$result ~ coc2$trees, col="black", pch="*")
  
  plot(model.pred ~ coc2$grass, col="goldenrod", pch=19, cex=2, ylim=ylimits, xlab="", ylab="", main="grass", xaxt="n", yaxt="n", cex.main=2)
  points(coc2$result ~ coc2$grass, col="black", pch="*")
  
  plot(model.pred ~ coc2$roofs, col="orange", pch=19, cex=2, ylim=ylimits, xlab="", ylab="", main="roofs", xaxt="n", yaxt="n", cex.main=2)
  points(coc2$result ~ coc2$roofs, col="black", pch="*")
  
  plot(model.pred ~ coc2$pm25, col="yellow", pch=19, cex=2, ylim=ylimits, xlab="", ylab="", main="pm25", xaxt="n", yaxt="n", cex.main=2)
  points(coc2$result ~ coc2$pm25, col="black", pch="*")
  
  plot(model.pred ~ coc2$traffic, col="gray", pch=19, cex=2, ylim=ylimits, xlab="", ylab="", main="traffic", xaxt="n", yaxt="n", cex.main=2)
  points(coc2$result ~ coc2$traffic, col="black", pch="*")
  
  plot(model.pred ~ coc2$rain, col="light blue", pch=19, cex=2, ylim=ylimits, xlab="", ylab="", main="rain", xaxt="n", yaxt="n", cex.main=2)
  points(coc2$result ~ coc2$rain, col="black", pch="*")
  
  boxplot(model.pred ~ coc2$landuse, col="light blue", at=c(2,4.5,7,9.5), xlim=c(0.75,9.75), ylim=ylimits, xlab="", ylab="", main="landuse", xaxt="n", yaxt="n", cex.main=2)
  boxplot(coc2$result ~ coc2$landuse, col="dark blue", at=c(1,3.5,6,8.5), add=TRUE, xaxt="n", yaxt="n")
  axis(side=1, at=c(1.5,4,6.5,9), labels=c("COM", "HDR", "IND", "LDR"))
  
  boxplot(model.pred ~ coc2$agency, col="light blue", at=c(2,4.5,7,9.5,12,14.5), xlim=c(0.75,14.75), ylim=ylimits, xlab="", ylab="", main="agency", xaxt="n", yaxt="n", cex.main=2)
  boxplot(coc2$result ~ coc2$agency, col="dark blue", at=c(1,3.5,6,8.5,11,13.5), add=TRUE, xaxt="n", yaxt="n")
  axis(side=1, at=c(1.5,4,6.5,9,11.5,14), labels=c("King", "Pierce", "POT", "Sea", "Sno", "Tac"))
}



#function to generate interaction plots by inter.quantile (IQ), levels of X, and grouping
Plot.Quantile <- function(aa, bb, myModel, xEqn=0.15, yEqn=max(coc$result)*0.9) {

  #create a list for quantiles of the two landscape predictors aa and bb; this is used to generate IQ
  cc <- list(unlist(qList[[aa]]), unlist(qList[[bb]]))
  names(cc) <- c(aa, bb)
  
  #generate dataframe of interaction fits for landscape predictors aa and bb
  IQ <- effect(paste(aa, "*", bb, sep=""), myModel, xlevels= cc, se=TRUE, confidence.level=0.95, typical=mean)
  IQ <- as.data.frame(IQ)
  IQ[bb] <- factor(as.numeric(unlist(IQ[bb])), levels=qList[[bb]], labels=c("0%", "25%", "50%", "75%", "100%"))
  
  #linear equation for the interaction, to be printed in the plot
  eqnText <- paste("ln(", this_param_short, ") = ", 
                  round(myModel$coefficients[[1]]["(Intercept)"], 2), " + ",
                  round(myModel$coefficients[[1]][aa], 2), "*", aa, " + ",
                  round(myModel$coefficients[[1]][bb], 2), "*", bb, " + ", 
                  round(myModel$coefficients[[1]][paste(aa,":",bb, sep="")], 2), "*", aa, "*", bb, sep="")

  #make the interaction plot  
  ggplot() + 
    geom_line(data=IQ, size=1.5, aes(x=IQ[,aa], y=fit, group=IQ[,bb], color=IQ[,bb]))+
    ylab("COC")+ xlab(names(IQ)[1])+
    ggtitle(paste("Interaction between", names(IQ)[1], "and", names(IQ)[2]))+
    scale_color_manual(values=c("blue", "purple", "yellow", "orange", "red"))+  #colors for lines (color)
    geom_point(data=coc2, aes(x=coc2[,aa], y=result, fill=coc2[,bb]), size=3, shape=21, stroke=0) + 
    scale_fill_gradient2(low = "blue", mid="yellow", high = "red")+  #colors for points (fill)
    geom_text(aes(x=xEqn, y=yEqn, label=eqnText), cex=4.5, color="black")+
    labs(color=names(IQ)[2], fill=names(IQ)[2])
}



#Function to plot result vs one predictor, with second predictor shown in color (no interaction)
Plot.Two.Predictors <- function(aa, bb, myModel, xEqn=0.15, yEqn=max(coc2$result)*0.9) {

  #linear equation for the primary predictor, to be printed in the plot
  eqnText <- paste("ln(", this_param_short, ") = ", 
                   round(myModel$coefficients[[1]]["(Intercept)"], 2), " + ",
                   round(myModel$coefficients[[1]][aa], 2), "*", aa)
  #plot the primary predictor on x-axis, with secondary predictor shown in color
  ggplot() + 
    geom_point(data=coc2, aes(x=coc2[,aa], y=result, fill=coc2[,bb]), size=3, shape=21, stroke=0) + 
    ylab("COC")+ xlab(aa) +
    scale_fill_gradient2(low = "blue", mid="yellow", high = "red")+  #colors for points (fill)
    geom_text(aes(x=xEqn, y=yEqn, label=eqnText), cex=4.5, color="black")+
    labs(color=bb, fill=bb) +
    geom_abline(intercept=myModel$coefficients[[1]]["(Intercept)"], slope=myModel$coefficients[[1]][aa], col="black", lwd=1.5)   #, aes(x = aa, y = y.hat), col="darkred", lwd = 1.5)
}    


#Function to plot result vs one predictor, with second predictor shown in color (no interaction)
Plot.Two.Predictors.Landuse <- function(aa, bb, myModel, xEqn=0.15, yEqn=max(coc2$result)*0.9) {
  
  #linear equation for the primary predictor, to be printed in the plot
  eqnText <- paste("ln(", this_param_short, ") = ", 
                   round(myModel$coefficients[[1]]["(Intercept)"], 2), " + ",
                   round(myModel$coefficients[[1]][aa], 2), "*", aa)
  #plot the primary predictor on x-axis, with secondary predictor shown in color
  ggplot() + 
    geom_point(data=coc2, aes(x=coc2[,aa], y=result, fill=coc2[,bb]), size=3, shape=21, stroke=0) + 
    ylab("COC")+ xlab(aa) +
    scale_fill_manual(values = c("blue", "purple", "yellow", "orange"))+
    geom_text(aes(x=xEqn, y=yEqn, label=eqnText), cex=4.5, color="black")+
    labs(color=bb, fill=bb) +
    geom_abline(intercept=myModel$coefficients[[1]]["(Intercept)"], slope=myModel$coefficients[[1]][aa], col="black", lwd=1.5)   #, aes(x = aa, y = y.hat), col="darkred", lwd = 1.5)
}

#Function to plot result vs one predictor, with location shown in color
Plot.One.Predictor.wAgency <- function(aa, bb, myModel, xEqn=0.15, yEqn=max(coc2$result)*0.9) {
  
  if(class(myModel)=="gls") {
    coef.int <- myModel$coefficients["(Intercept)"]
    coef.slope <- myModel$coefficients[aa]
  } else {
    coef.int <- myModel$coefficients[[1]]["(Intercept)"]
    coef.slope <- myModel$coefficients[[1]][aa]
  }
  #linear equation for the primary predictor, to be printed in the plot
  eqnText <- paste("ln(", this_param_short, ") = ", 
                   round(coef.int, 2), " + ",
                   round(coef.slope, 2), "*", aa )
  #plot the primary predictor on x-axis, with secondary predictor shown in color
  ggplot() + 
    geom_point(data=coc2, aes(x=coc2[,aa], y=result, fill=coc2[,bb]), size=3, shape=21, stroke=0) + 
    ylab("COC")+ xlab(aa) +
    scale_fill_manual(values = c("red", "orange", "yellow", "green", "blue", "purple"))+
    geom_text(aes(x=xEqn, y=yEqn, label=eqnText), cex=4.5, color="black")+
    labs(color=bb, fill=bb) +
    geom_abline(intercept=coef.int, slope=coef.slope, col="black", lwd=1.5)   #, aes(x = aa, y = y.hat), col="darkred", lwd = 1.5)
}    


#------------------------------------------------------#
#  Functions to make choosing a best fit model easier  #
#------------------------------------------------------#

#function to check correlations for a model's predictors
check.cor <- function(theModel) {
  aa <- summary(theModel)$coefficients$fixed
  aa <- aa[names(aa) %in% preds.1]
  pairs(select(coc2, names(aa)),
        lower.panel=panel.smooth2, upper.panel=panel.cor, diag.panel=panel.hist) 
}

#function to plot best fit GLS or LME model line to a single predictor
plot.single.preds <- function(theModel) {
  if(class(theModel)=="gls")  aa <- summary(theModel)$coefficients
  if(class(theModel)=="lme")  aa <- summary(theModel)$coefficients$fixed
  aa <- names(aa[names(aa) %in% preds.1])
  ppp1 <- Plot.One.Predictor.wAgency(aa[1], "agency", theModel, yEqn=max(coc2$result)*1.02, xEqn=0)
  if(length(aa) > 1)  ppp2 <- Plot.One.Predictor.wAgency(aa[2], "agency", theModel, yEqn=max(coc2$result)*1.02, xEqn=0)
  if(length(aa) > 2)  ppp3 <- Plot.One.Predictor.wAgency(aa[3], "agency", theModel, yEqn=max(coc2$result)*1.02, xEqn=0)
  if(length(aa) > 3)  ppp4 <- Plot.One.Predictor.wAgency(aa[4], "agency", theModel, yEqn=max(coc2$result)*1.02, xEqn=0)
  
  if(length(aa)==1)  grid.arrange(ppp1, nrow=2, ncol=2)
  if(length(aa)==2)  grid.arrange(ppp1, ppp2, nrow=2, ncol=2)
  if(length(aa)==3)  grid.arrange(ppp1, ppp2, ppp3, nrow=2, ncol=2)
  if(length(aa)==4)  grid.arrange(ppp1, ppp2, ppp3, ppp4, nrow=2, ncol=2)
}


#function to check if coefficients are in the direction they should be

my.aics <- rep(0, length(ee))
for (i in 1:length(ee)) {
  bb <- lme(data=coc2, ee[[i]], random = r1X, method="ML", weights=vf1X, control = lmeControl(maxIter = 1e8, msMaxIter = 1e8))
  my.aics[i] <- AIC(bb)
}

aa <- summary(bb)$coefficients$fixed
aa <- aa[names(aa) %in% preds.1]
sign(aa)


