Correlation Tables
================

## List of Parameters

``` r
params %>% 
kbl(digits = 2, col.names = "Parameter",) %>%
  kable_classic(full_width = F,html_font = "serif")
```

<table class=" lightable-classic" style="font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Parameter
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
Copper - Water - Total
</td>
</tr>
<tr>
<td style="text-align:left;">
CPAH - Water - Total
</td>
</tr>
<tr>
<td style="text-align:left;">
HPAH - Water - Total
</td>
</tr>
<tr>
<td style="text-align:left;">
Lead - Water - Total
</td>
</tr>
<tr>
<td style="text-align:left;">
Nitrite-Nitrate - Water - Dissolved
</td>
</tr>
<tr>
<td style="text-align:left;">
Total Kjeldahl Nitrogen - Water - Total
</td>
</tr>
<tr>
<td style="text-align:left;">
Total PAH - Water - Total
</td>
</tr>
<tr>
<td style="text-align:left;">
Total Phosphorus - Water - Total
</td>
</tr>
<tr>
<td style="text-align:left;">
Total Phthalate - Water - Total
</td>
</tr>
<tr>
<td style="text-align:left;">
Total Suspended Solids - Water - Total
</td>
</tr>
<tr>
<td style="text-align:left;">
Zinc - Water - Dissolved
</td>
</tr>
<tr>
<td style="text-align:left;">
Zinc - Water - Total
</td>
</tr>
</tbody>
</table>

``` r
cortable <- function(i){
  

p = params[i]
df.coc <- joined %>%
  filter(parameter == p) %>%
  mutate(log_concentration = log(result)) 
  #select(log_concentration,everything()) %>%
 # mutate_at(vars(starts_with('CO_emissions')), log) 
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
  kbl(digits = 2, escape = FALSE, caption = paste("Correlation table",p)) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>% 
  kable_classic(full_width = F,html_font = "serif")


}
for(i in 1:length(params)){
  
print(htmltools::h2(params[i])) 
       
print(cortable(i))
}
```

<h2>
Copper - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Copper - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.53
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.60
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.72
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
CPAH - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table CPAH - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.51
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.58
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.71
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
HPAH - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table HPAH - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.51
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.58
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.71
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
Lead - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Lead - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.71
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.58
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.51
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
Nitrite-Nitrate - Water - Dissolved
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Nitrite-Nitrate - Water - Dissolved
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.73
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.61
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.54
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
Total Kjeldahl Nitrogen - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Total Kjeldahl Nitrogen - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.55
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.62
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.73
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
Total PAH - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Total PAH - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.51
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.58
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.71
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
Total Phosphorus - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Total Phosphorus - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.57
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.63
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.74
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.90
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
Total Phthalate - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Total Phthalate - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.22
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.51
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.58
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.71
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
Total Suspended Solids - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Total Suspended Solids - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.61
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.54
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.15
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.73
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.73
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
Zinc - Water - Dissolved
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Zinc - Water - Dissolved
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.89
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.28
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.83
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.21
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1975\_1990
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.20
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
RES
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_1990\_2000
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.55
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.08
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_2000\_2014
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.31
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
grass\_low\_veg
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.62
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.22
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.54
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.70
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.25
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.57
</td>
</tr>
<tr>
<td style="text-align:left;">
slope
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
no\_dev
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.92
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
RURES
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.20
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.80
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
percent\_tree\_cover
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.71
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.89
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.09
</td>
</tr>
<tr>
<td style="text-align:left;">
precip
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
inches\_rain\_per\_hour
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
</tr>
</tbody>
</table>
<h2>
Zinc - Water - Total
</h2>
<table class="table table-striped table-hover table-condensed lightable-classic" style="margin-left: auto; margin-right: auto; font-family: serif; width: auto !important; margin-left: auto; margin-right: auto;">
<caption>
Correlation table Zinc - Water - Total
</caption>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
log\_concentration
</th>
<th style="text-align:right;">
COM
</th>
<th style="text-align:right;">
CO\_emissions\_commercial
</th>
<th style="text-align:right;">
CO\_emissions\_nonroad
</th>
<th style="text-align:right;">
CO\_emissions\_onroad
</th>
<th style="text-align:right;">
CO\_emissions\_residential
</th>
<th style="text-align:right;">
CO\_emissions\_total
</th>
<th style="text-align:right;">
IND
</th>
<th style="text-align:right;">
NO\_2
</th>
<th style="text-align:right;">
PM25\_NA
</th>
<th style="text-align:right;">
RES
</th>
<th style="text-align:right;">
ROW
</th>
<th style="text-align:right;">
RURES
</th>
<th style="text-align:right;">
dev\_1975\_1990
</th>
<th style="text-align:right;">
dev\_1990\_2000
</th>
<th style="text-align:right;">
dev\_2000\_2014
</th>
<th style="text-align:right;">
dev\_pre\_1975
</th>
<th style="text-align:right;">
grass\_low\_veg
</th>
<th style="text-align:right;">
imperv\_ground
</th>
<th style="text-align:right;">
imperv\_roofs
</th>
<th style="text-align:right;">
no\_dev
</th>
<th style="text-align:right;">
particulate\_surface\_area
</th>
<th style="text-align:right;">
percent\_tree\_cover
</th>
<th style="text-align:right;">
pm25
</th>
<th style="text-align:right;">
pop\_per\_ha
</th>
<th style="text-align:right;">
roof\_COM
</th>
<th style="text-align:right;">
roof\_IND
</th>
<th style="text-align:right;">
roof\_RES
</th>
<th style="text-align:right;">
slope
</th>
<th style="text-align:right;">
traffic
</th>
<th style="text-align:right;">
precip
</th>
<th style="text-align:right;">
inches\_rain\_per\_hour
</th>
<th style="text-align:right;">
yday
</th>
<th style="text-align:right;">
month
</th>
<th style="text-align:right;">
year
</th>
<th style="text-align:right;">
season
</th>
<th style="text-align:right;">
latitude
</th>
<th style="text-align:right;">
longitude
</th>
<th style="text-align:right;">
daymet\_precip
</th>
<th style="text-align:right;">
daymet\_2day
</th>
<th style="text-align:right;">
daymet\_3day
</th>
<th style="text-align:right;">
daymet\_5day
</th>
<th style="text-align:right;">
daymet\_7day
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
log\_concentration
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.14
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_total
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.84
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.86
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
traffic
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
dev\_pre\_1975
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.93
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
</tr>
<tr>
<td style="text-align:left;">
IND
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.61
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_IND
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.96
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.07
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_nonroad
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.54
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.94
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.82
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_onroad
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.81
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_ground
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.76
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_commercial
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
-0.40
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
</tr>
<tr>
<td style="text-align:left;">
particulate\_surface\_area
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.77
</td>
<td style="text-align:right;">
0.38
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.64
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.78
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.12
</td>
</tr>
<tr>
<td style="text-align:left;">
pm25
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.87
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.37
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.03
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_COM
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
-0.79
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.46
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.77
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.13
</td>
</tr>
<tr>
<td style="text-align:left;">
COM
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.72
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.30
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.70
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
PM25\_NA
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.32
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
-0.39
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.60
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
ROW
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.88
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
-0.75
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
</tr>
<tr>
<td style="text-align:left;">
imperv\_roofs
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.93
</td>
<td style="text-align:right;">
-0.85
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.83
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
</tr>
<tr>
<td style="text-align:left;">
pop\_per\_ha
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.86
</td>
<td style="text-align:right;">
-0.74
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
-0.72
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.87
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.18
</td>
</tr>
<tr>
<td style="text-align:left;">
CO\_emissions\_residential
</td>
<td style="text-align:right;">
0.29
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.78
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
-0.65
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
NO\_2
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.60
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.28
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.57
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.50
</td>
<td style="text-align:right;">
-0.58
</td>
<td style="text-align:right;">
0.47
</td>
<td style="text-align:right;">
-0.49
</td>
<td style="text-align:right;">
0.53
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.37
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.52
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
month
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
yday
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
season
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.02
</td>
</tr>
<tr>
<td style="text-align:left;">
roof\_RES
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.61
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.64
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.80
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.48
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
-0.59
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
-0.35
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.36
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.19
</td>
</tr>
<tr>
<td style="text-align:left;">
year
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.43
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.06
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_5day
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
-0.01
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.88
</td>
</tr>
<tr>
<td style="text-align:left;">
latitude
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.29
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.34
</td>
<td style="text-align:right;">
-0.51
</td>
<td style="text-align:right;">
-0.27
</td>
<td style="text-align:right;">
-0.66
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
-0.48
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
-0.44
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.00
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
-0.16
</td>
</tr>
<tr>
<td style="text-align:left;">
longitude
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.57
</td>
<td style="text-align:right;">
-0.33
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.47
</td>
<td style="text-align:right;">
-0.42
</td>
<td style="text-align:right;">
-0.63
</td>
<td style="text-align:right;">
-0.26
</td>
<td style="text-align:right;">
-0.56
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
0.34
</td>
<td style="text-align:right;">
-0.67
</td>
<td style="text-align:right;">
0.42
</td>
<td style="text-align:right;">
-0.68
</td>
<td style="text-align:right;">
-0.50
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
-0.41
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
-0.69
</td>
<td style="text-align:right;">
-0.38
</td>
<td style="text-align:right;">
-0.45
</td>
<td style="text-align:right;">
-0.52
</td>
<td style="text-align:right;">
-0.32
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
-0.53
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.05
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_precip
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.21
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.18
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
0.13
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.12
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.17
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.76
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.51
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_2day
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.22
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.23
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.24
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.20
</td>
<td style="text-align:right;">
0.21
</td>
<td style="text-align:right;">
-0.24
</td>
<td style="text-align:right;">
0.23
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.09
</td>
<td style="text-align:right;">
-0.09
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
0.15
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.58
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_7day
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.10
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
-0.07
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.02
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.08
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.04
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
0.03
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
0.10
</td>
<td style="text-align:right;">
-0.12
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
-0.03
</td>
<td style="text-align:right;">
-0.18
</td>
<td style="text-align:right;">
-0.13
</td>
<td style="text-align:right;">
0.07
</td>
<td style="text-align:right;">
-0.19
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:right;">
-0.06
</td>
<td style="text-align:right;">
-0.05
</td>
<td style="text-align:right;">
0.06
</td>
<td style="text-align:right;">
-0.02
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.70
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
1.00
</td>
</tr>
<tr>
<td style="text-align:left;">
daymet\_3day
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.11
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.17
</td>
<td style="text-align:right;">
-0.14
</td>
<td style="text-align:right;">
-0.16
</td>
<td style="text-align:right;">
-0.15
</td>
<td style="text-align:right;">
-0.06