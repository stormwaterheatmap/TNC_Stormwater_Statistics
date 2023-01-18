library(raster)

fn <- ("~/Documents/repos/TNC_Stormwater_Statistics/data/predictors_30m.tif")
x <- raster(fn,band=1)
y <- raster(fn,band=2)
z <- raster(fn,band=3)

r1 <- raster::stack(c(x,y,z))

fun=function(x,y,z){return((x + y) * z)}


b <- overlay(r1, fun=function(x,y,z){return(x*y*z)} )
