# Explore how to determine% change in PVB from a 1% change in disc rate by applying the ratio of AAL to PVB to the corresponding
# % change in AAL(based on TPL) 

library("ggplot2")
library("scales") # so we can use scales with ggplot2
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("magrittr")
library("lubridate")
library("tibble")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
library("tidyr")
library("dplyr") # always load AFTER plyr
options(dplyr.print_min = 60) # default is 10
options(dplyr.print_max = 60) # default is 20
library("knitr")
library("stringr")
library("grDevices")
library("readr")
library("readxl")
library(XLConnect)




# Load mortality data
data_raw_healthyMale   <- readWorksheetFromFile("RP2000&ScaleBB/RP2000-male-CombinedHealthy.xlsx", 
                                                sheet = "Sheet1", header=FALSE, region="a25:b144", colTypes="character") %>% 
                                                mutate_each(funs(as.numeric))

data_raw_healthyMale



