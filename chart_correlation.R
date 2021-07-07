library(tidyverse)
s8data <- read.csv("~/repos/TNC_Stormwater_Statistics/data/s8data.csv")
df_predictors <- read.csv(here::here("data","spatial_predictors_raw.csv"))
# join data ---------------------------------------------------------------
joined <- df_predictors %>% mutate(
  location_id = Location_N
) %>% left_join(s8data, by = "location_id") %>% as.data.frame()



# Correlation -------------------------------------------------------------
library(PerformanceAnalytics)
params = joined$parameter %>% unique() %>% sort()





cortable <- function(i){
  

p = params[i]
df.coc <- joined %>%
  filter(parameter == p) %>%
  mutate(log_concentration = log(result)) %>%
  #select(log_concentration,everything()) %>%
  mutate_at(vars(starts_with('CO_emissions')), log) 
#replace nan with 0 
#df.coc[is.nan(df.coc)] = 0
#df.coc[is.na(df.coc)]=0

M = df.coc %>% select_if(is.numeric) %>% dplyr::select(-access_id,-result)%>%
  dplyr::select(c(log_concentration, everything())) %>% 
  #select(-c(starts_with("roof"),"OPEN")) %>% 
  dplyr::mutate_if(is.numeric, scale)
  
#cols = c(1,3,11:21)
cols = c(1:43)
M[,cols] %>% 
  cor(.) %>% 
  as.data.frame() %>% 
  arrange(desc(log_concentration)) %>% 
  kbl(digits = 2, caption = paste("Correlation table",p)) %>%
  kable_classic(full_width = F,html_font = "serif")
}
cortable(1)
