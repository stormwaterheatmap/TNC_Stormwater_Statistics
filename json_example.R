library(tidyverse)
load("/Users/christiannilsen/Downloads/Copper_diff_U80.Rdata")
sample <- sample_n(
  Cu.d.U80,100)


#plot(sample,pch=19,cex=0.1)
#Round cu 



# 
# #make list of values 
# df.1 <- sample %>% group_by(devAge2,sqrt_traffic) %>% 
#   nest() %>% 
#   distinct(devAge2,sqrt_traffic,.keep_all=TRUE) 
# df.2 <- df.1 %>% group_by(devAge2) %>% nest()
# %>% nest() %>% group_by(devAge2) %>% nest()
# library(jsonlite)
# 
# 
# df.1[[3]][[1]] %>% pull(1,2)
# 
# temp <- stringr::str_replace_all(js_out, '"(\\w+)"\\s*:', '\\1:')
# 
# jsonlite::write_json(temp,'~/Documents/repos/TNC_Stormwater_Statistics/test_json.json')
# sample
l = list()
for (i in 1:nrow(sample)) {
  l[i]=sample[i,1]
  
}
names(l)<-df.2$devAge2

l %>% toJSON()
