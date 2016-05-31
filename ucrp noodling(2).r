

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

library("btools") # library that I created (install from github) - not needed here
library("bdata") # library that I created (install from github) - not needed here
library("pdata") # library that I created (install from github) - not needed here
library("fof") # library that I created (install from github) - not needed here
library("apitools") # library that I created (install from github) - not needed here
library("ipoptr") # get from dropbox


df <- read_excel(paste0("./data-raw/", "ucrp_djb.xlsx"), sheet="2015UCRP_djb")
glimpse(df)
ht(df)

dr <- 0.0725
df2 <- df %>% mutate(year=year - min(year) + 1) %>%
  gather(variable, value, -year) %>%
  mutate(value=ifelse(is.na(value), 0, value),
         value=value / 1e6,
         pv=value / (1 + dr)^(year - 1),
         pv.m1=value / (1 + dr - .01)^(year - 1),
         pv.p1=value / (1 + dr + .01)^(year - 1))
glimpse(df2)

pvsums <- df2 %>% group_by(variable) %>%
  summarise_each(funs(sum(., na.rm=TRUE)), contains("pv")) %>%
  mutate(pch.m1=pv.m1 / pv * 100 - 100,
         pch.p1=pv.p1 / pv * 100 - 100,
         duration=(pv.m1 - pv.p1) / (2*pv) * 100)
pvsums


df2 %>% filter(str_detect(variable, "active")) %>% qplot(year, value, data=., colour=variable, geom=c("point", "line"))
df2 %>% filter(str_detect(variable, "termv")) %>% qplot(year, value, data=., colour=variable, geom=c("point", "line"))
df2 %>% filter(str_detect(variable, "retire")) %>% qplot(year, value, data=., colour=variable, geom=c("point", "line"))

# retired looks easiest, so let's start with that
# step 1: shift the yy curve so that it has the same year1 value as known retirees
ret1 <- df2 %>% filter(str_detect(variable, "retire")) %>%
  select(year, variable, value) %>%
  spread(variable, value) %>%
  mutate(rstep1=retired.yy * retired[year==1] / retired.yy[year==1])
ret1 %>% gather(variable, value, -year) %>% qplot(year, value, data=., colour=variable, geom=c("point", "line"))


# step 2: adjust all values after year 1 so that they have the same pvsum as the target
# recompute pvs as that is easy
pvf <- function(vec, dr) {
  pv <- vec / ((1+dr)^((1:length(vec))-1))
  pv <- sum(pv, na.rm=TRUE)
  return(pv)
}

duration <- function(vec, dr) {
  pv <- pvf(vec, dr)
  pv.m1 <- pvf(vec, dr - .01)
  pv.p1 <- pvf(vec, dr + .01)
  duration <- (pv.m1 - pv.p1) / (2*pv) * 100
  return(duration)
}


ratio <- pvf(ret1$retired[ret1$year>1], dr) / pvf(ret1$rstep1[ret1$year>1], dr)
ret2 <- ret1 %>% mutate(rstep2=ifelse(year>1, rstep1 * ratio, rstep1))
ret2 %>% gather(variable, value, -year) %>% qplot(year, value, data=., colour=variable, geom=c("point", "line"))


pvf(ret2$retired, dr)
pvf(ret2$retired.yy, dr)
pvf(ret2$rstep1, dr)
pvf(ret2$rstep2, dr)

duration(ret2$retired, dr)
duration(ret2$retired.yy, dr)
duration(ret2$rstep1, dr)
duration(ret2$rstep2, dr)

# step 3: adjust step2 values so that they have the same duration for +/- .01

eval_f <- function(step3, step2, dr) {
  # objective function - evaluates to a single number
  # compute distortion from step 3 values - the further the movement from step 2, the greater the penalty
  # x (solution) must be in vector form
  # simple obj function: least squares of sum of differences weighted by an indicator of how far
  # the x is moving from the center of the area
  obji <- (step3 - step2)^2
  obj <- sum(obji)
  return(obj)
}

eval_grad_f<-function(step3, step2, dr){
  # gradient of objective function - a vector length x 
  # giving the partial derivatives of obj wrt each x[i]
  # must pass unnecessary parameters because ipoptr requires all functions to get the same parameters
  gradf <- 2 * (step3 - step2)
  return(gradf)
}

eval_g <- function(step3, step2, dr) {
  # constraints that must hold in the solution - just give the lhs of each constraint expression
  # return a vector where each element evaluates a constraint
  # must pass unnecessary parameters because ipoptr requires all functions to get the same parameters
  # we calculate 3 constraints: pv at dr, and at dr + and - .01
  c1.pv <- pvf(step3, dr)
  c2.pv.m1 <- pvf(step3, dr - .01)
  c3.pv.p1 <- pvf(step3, dr + .01)
  c4.year1 <- step3[1]
  constraints <- c(c1.pv, c2.pv.m1, c3.pv.p1, c4.year1)
  return(constraints)
}

eval_jac_g <- function(step3, step2, dr){
  # function to construct the Jacobian matrix of first partial derivatives of constraints at point x
  # one column per constraint
  # one row per x variable
  # easy here because constraints are linear
  # can be set up as dense (one value for each element in the matrix, whether or not the value is zero)
  # or sparse (one value for each nonzero element in the matrix)
  # this problem is small enough that dense will be fine
  pvchange <- function(vec, dr) return(1 / ((1+dr)^((1:length(vec))-1)))
  c1.pv.derivs <- pvchange(step3, dr)
  c2.pv.m1.derivs <- pvchange(step3, dr - .01)
  c3.pv.p1.derivs <- pvchange(step3, dr + .01)
  c4.year1.derivs <- c(1, rep(0, length(step3) - 1))

  derivs <- c(c1.pv.derivs, c2.pv.m1.derivs, c3.pv.p1.derivs, c4.year1.derivs)
  return(derivs) 
}

# define jacobian structure - a list with one item (a vector) for each consrtraint
jaclen <- nrow(ret2)
eval_jac_g_structure <- list(1:jaclen, 1:jaclen, 1:jaclen, 1:jaclen) # The Jacobian for this problem is dense


# now solve for the step3 values
# use these options to test the first derivatives against ipopt finite difference approach
optlist <- list("file_print_level"=5,
                "output_file"="test.out",
                "derivative_test"="first-order",
                "derivative_test_print_all"="yes",
                "max_iter"=20)

# once derivatives are correct, use these options
optlist <- list("file_print_level"=5,
                "output_file"="test.out")


targets <- c(pvf(ret2$retired, dr), pvf(ret2$retired, dr -.01), pvf(ret2$retired, dr + .01), ret2$retired[1])

result <- ipoptr(# x values - starting point and bounds (x is the new matrix to find, in vector shape)
  x0=ret2$rstep2,
  lb=rep(0, length(ret2$rstep2)),
  ub=rep(1e6, length(ret2$rstep2)), # make step3 upper bounds as large as we might plausibly need
  
  # objective function and its gradient
  eval_f=eval_f,
  eval_grad_f=eval_grad_f,
  
  # constraint function, constraint values, and constraint Jacobian
  eval_g=eval_g,
  # set constraint lower and upper bounds to the same values as we would like to hit them exactly
  constraint_lb=targets, 
  constraint_ub=targets,
  eval_jac_g=eval_jac_g,
  eval_jac_g_structure=eval_jac_g_structure,
  
  # Hessian - let's try running without the 
  # eval_h=eval_h,
  # eval_h_structure=eval_h_structure,
  
  # additional parameters passed to functions
  step2=ret2$rstep2,
  dr=dr,
  
  opts=optlist)


str(result)

ret3 <- ret2 %>% mutate(rstep3=result$solution)
ret3 %>% gather(variable, value, -year) %>% qplot(year, value, data=., colour=variable, geom=c("point", "line"))


pvf(ret3$retired, dr)
pvf(ret3$rstep3, dr)

duration(ret3$retired, dr)
duration(ret3$rstep2, dr)
duration(ret3$rstep3, dr)




