library(tidyverse)
load("/Users/christiannilsen/Downloads/Copper_diff_U80.Rdata")
sample <- sample_n(
  Cu.d.U80,100000)
head(sample)

dev <- sample$devAge2 %>% sort(FALSE)

devages <- c(
   # -3.186779260635376,
  #  -2.924546241760254,
  #  -2.137847423553467,
  -1.9999,
  0.9999,
    -0.8266825675964355,
  0.4
  #  1.008948206901550
)

devages <- seq(from=-1.75,to=.75,length.out = 5)


pos <- vector()

for (i in 1:length(devages)) {
  x <- devages[i]
  print(x)
  val <- min(which(
    dev>=x
  ))
  pos[i]=val
    
}
truevals <- sample %>% filter(devAge2 %in% dev[pos])
plot(d.U80~sqrt_traffic,truevals %>% sample_n(100))

forplot <- truevals %>% sample_n(1000) %>% 
  mutate(DevAge2 = devAge2 %>% round(2) %>%   as.factor()) %>% #select(c(d.U80,devAge2,sqrt_traffic)) %>% 
  pivot_longer(cols = -c(d.U80,DevAge2))

ggplot(forplot,aes(x=value,y=d.U80,color=DevAge2,group=DevAge2))+
  geom_point(alpha=0.5,shape=16,size=0.6)+
  geom_smooth(aes(group=DevAge2),size=0.1)+
  #scale_y_log10()+#limits = c(0.8,1.2))+
  facet_wrap(~name)

#try a linear model 
mod1 <- lm((d.U80)~sqrt_traffic+poly(pm25_na, 3)+devAge2,data=truevals %>% sample_n(1000))
mod1
summary(mod1)
plot(mod1)

subset <- truevals %>% sample_n(10000)
## An example of polynomial regression
plot(d.U80~sqrt_traffic, subset,cex=0.1,pch=20,col = "cyan")
d <- seq(-1, 1, length.out = 200)
for(degree in 1:4) {
  fm <- lm(d.U80 ~ devAge2+poly(sqrt_traffic, degree), data = subset)
  assign(paste("subset", degree, sep = "."), fm)
 lines(d, predict(fm, data.frame(sqrt_traffic = d,devAge2=-1), col = degree))
}

