

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

library("btools") # library that I created (install from github) - maybe not needed here
library("pdata") # library that I created (install from github) - not needed here
library("ipoptr") # get from dropbox


#****************************************************************************************************
#                Globals ####
#****************************************************************************************************
dr <- 0.0725

#****************************************************************************************************
#                Functions ####
#****************************************************************************************************
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


#****************************************************************************************************
#                Get and adjust ucrp data ####
#****************************************************************************************************
df <- read_excel(paste0("ucrp_djb.xlsx"), sheet="2015UCRP_djb")
glimpse(df)
ht(df)
# check each variable for plausibility
df %>% ggplot(aes(year, actives)) + geom_line() # good
df %>% ggplot(aes(year, termv)) + geom_line() # first year termv not believable, lets apply year2 to year3 growth to get adjusted year1
df %>% ggplot(aes(year, retired)) + geom_line() # strange beginning, but we'll keep it

# adjust for oddball values
df1a <- df %>% arrange(year) %>%
  mutate(year=as.integer(year), ratio=termv / lead(termv)) %>%
  mutate(termv=ifelse(year==2015, ratio[year==2016] * termv[year==2016], termv)) %>%
  select(-ratio)
head(df1a)
df1a %>% ggplot(aes(year, actives)) + geom_line() # good
df1a %>% ggplot(aes(year, termv)) + geom_line() # first year termv not believable, lets apply year2 to year3 growth to get adjusted year1
df1a %>% ggplot(aes(year, retired)) + geom_line() # strange beginning, but we'll keep it

df2 <- df1a %>% mutate_each(funs(ifelse(is.na(.), 0, . / 1e6)), -year) %>%
  arrange(year) %>%
  mutate(year=row_number())


#****************************************************************************************************
#                Just playing around - fit ucrp data ####
#****************************************************************************************************
vsmean <- function(vec) {
  # function to index values to their mean of nonzeros
  xzmean <- mean(vec[vec>0])
  vsmean <- vec / xzmean
  return(vsmean)
}
# vsmean(c(0, 0, 0, 4, 5, 6))

# df3 <- df2 %>% mutate_each(funs(vsmean), -year)
df3 <- df2 %>% mutate_each(funs(. / .[1]), -year) # compute all values relative to first value
df3
glimpse(df3)

# pick one of the following to work with
df4 <- df3 %>% select(year, y=actives)
df4 <- df3 %>% select(year, y=termv)
df4 <- df3 %>% select(year, y=retired)
df4 %>% ggplot(aes(year, y)) + geom_line()

# Gaussian is probably best for this
# The model is y = b + a * exp(−.5 * ((x−m) / s)^2) + u
# where u is error. 
# m:  x value at which the curve peaks
# s: s > 0 the rate at which the curve tapers off
# a: a > 0 max value of y
# b: baseline
rhs <- function(x, m, s, a, b=0) {b + a * exp(-0.5 * ((x - m)/s)^2)}

# define starting values
# ucrp info:
#   actives m s a: ~35, 16, 2
#   termv m s a: ~23, 17, 2
#   retired m s a: ~8, 15, 4
# starting values of 20, 15, 3 seem good all around
m3 <- nls(y ~ rhs(year, m, s, a), data = df4, start = list(m=20, s=15, a=13), trace = TRUE)
summary(m3)

df4 %>% mutate(yfit=predict(m3)) %>% 
  gather(variable, value, -year) %>%
  filter(variable %in% c("y", "yfit")) %>%
  ggplot(aes(year, value, colour=variable)) +
  geom_line()
rhs(30:40, coef(m3)["m"], coef(m3)["s"], coef(m3)["a"])

# it might be better to use some sort of exponential decay function for retirees? but gaussian could be ok

# what if we perturb the parameters?
df5 <- df4 %>% mutate(yfit=predict(m3),
                      yfitadj=rhs(year, coef(m3)["m"]-3, coef(m3)["s"]+2, coef(m3)["a"]*1.1)) %>% 
  gather(variable, value, -year) %>%
  filter(variable %in% c("y", "yfit", "yfitadj"))
df5 %>% ggplot(aes(year, value, colour=variable)) + geom_line()
# now let's push yfit and yfitadj around and see how well we measure up against known y


#****************************************************************************************************
#                Create an initial guess based on model fit ####
#****************************************************************************************************
# The purpose here is to get an initial guess of annual values. In a later step we will adjust this with
# an NLP to hit constraints while minimizing a penalty function for large movements from the guess.
# In a real setting, we would not have "actual" (Segal estimates) data to base our guess on. We prob
# would base the guess on parameters from models fitted to simulated data. But here we form the intial
# guess based on curve fitting (to a gaussian curve, since that is how actives and termvesteds look)
# and to get a better sense of how the approach works, we also perturb the fitted values so we can
# see what happens when we optimize based on a guess that is perturbed significantly. 

# Gaussian is prob best for this
# http://stats.stackexchange.com/questions/83022/how-to-fit-data-that-looks-like-a-gaussian
# http://stats.stackexchange.com/questions/70153/linear-regression-best-polynomial-or-better-approach-to-use
# The model is y = b + a * exp(−.5 * ((x−m) / s)^2) + u
# where u is error. 
# m:  x value at which the curve peaks
# s: s > 0 the rate at which the curve tapers off; when s is bigger the curve is more spread out
# a: a > 0 max value of y
# b: baseline
rhs <- function(x, m, s, a, b=0) {b + a * exp(-0.5 * ((x - m)/s)^2)}
x <- 1:100; y <- rhs(x, 30, 20, 10); qplot(x, y, geom="line") # explore the gaussian shape to help id starting values

m3 <- nls(actives ~ rhs(year, m, s, a), data = df3, start = list(m=20, s=15, a=13), trace = TRUE)
summary(m3)

# perturb the parameters and create alternative guess(es)
df5 <- df3 %>% select(year, actives, actives.yy=actives_shift.yy) %>%
  mutate(actives.fit=predict(m3),
         actives.fit1=rhs(year, coef(m3)["m"]-3, coef(m3)["s"]+2, coef(m3)["a"]*1.1),
         actives.guess=rhs(year, 30, 17, 13)) %>% 
  gather(variable, value, -year)
df5 %>% ggplot(aes(year, value, colour=variable)) + geom_line()
# now let's push fit and adjusted fit around and see how well we measure up against known (Segal) actives


#****************************************************************************************************
#                Prepare the data for optimization ####
#****************************************************************************************************
glimpse(df2) # df2 is adjusted for oddballs, has no na values, and values are in $m
glimpse(df3) # df3 is the same, as ratio to mean
glimpse(df5) # df5 has actual, fit, and perturbed fit
df5 %>% filter(year==1) # note that they have different year1 values

# start by getting pv of the ratio-to-mean data
pv.years <- df5 %>% group_by(variable) %>%
  mutate(pv=value / (1 + dr)^(year - 1),
         pv.m1=value / (1 + dr - .01)^(year - 1),
         pv.p1=value / (1 + dr + .01)^(year - 1))
glimpse(pv.years)

pvsums <- pv.years %>% summarise_each(funs(sum), contains("pv")) %>%
  mutate(pch.m1=pv.m1 / pv * 100 - 100,
         pch.p1=pv.p1 / pv * 100 - 100,
         duration=(pv.m1 - pv.p1) / (2*pv) * 100)
pvsums

# for each alternative, get ratio of {pv for years 2-end} to same for actual actives
pvratios <- pv.years %>% filter(year==1) %>%
  select(variable, year1=value) %>%
  left_join(select(pvsums, variable, pv)) %>% 
  ungroup %>%
  mutate(year2plus=pv - year1, year2pratio=year2plus / year2plus[variable=="actives"])
pvratios

# steps 1 and 2: shift each fit curve so that it has the same year1 value as known actives and so that it has same pv
steps1and2 <- df5 %>% mutate(step1=ifelse(year==1, pvratios$year1[pvratios$variable=="actives"], value)) %>%
  mutate(step2=ifelse(year>1, step1 / pvratios$year2pratio[match(variable, pvratios$variable)], step1))
# verify that year1 values are the same
steps1and2 %>% filter(year<=3)
# verify that pv of each step2 group is the same
steps1and2 %>% group_by(variable) %>% summarise(pv=pvf(step2, dr)) # good
steps1and2 %>% group_by(variable) %>% summarise(dur=duration(step2, dr)) # perturbed fits have very different durations

# graphs
glimpse(steps1and2)
steps1and2 %>% ggplot(aes(year, step2, colour=variable)) + geom_line()

# verify that resulting data has the desired pv
steps1and2 %>% group_by(variable) %>%
  arrange(year) %>%
  summarise_each(funs(pv=pvf(., dr) %>% round(2)), value, step1, step2)


# step 3: adjust step2 values so that they have the same duration for +/- .01
#****************************************************************************************************
#                Objective functions, gradients, etc. needed by ipoptr ####
#****************************************************************************************************
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

eval_jac_g_structure_fn <- function(df) {
  # define jacobian structure - a list with one item (a vector) for each consrtraint
  jaclen <- nrow(df)
  eval_jac_g_structure <- list(1:jaclen, 1:jaclen, 1:jaclen, 1:jaclen) # The Jacobian for this problem is dense
  return(eval_jac_g_structure)
}


callipoptr <- function(optdf, targets, dr, optlist=NULL) {
    # set up a function to call ipoptr so that we can call multiple times
    
    if(is.null(optlist)) optlist <- list("file_print_level"=5, "output_file"="test.out")
    
    result <- ipoptr(# x values - starting point and bounds (x is the new matrix to find, in vector shape)
      x0=optdf$value,
      lb=rep(0, length(optdf$value)),
      ub=rep(1e6, length(optdf$value)), # make step3 upper bounds as large as we might plausibly need
      
      # objective function and its gradient
      eval_f=eval_f,
      eval_grad_f=eval_grad_f,
      
      # constraint function, constraint values, and constraint Jacobian
      eval_g=eval_g,
      # set constraint lower and upper bounds to the same values as we would like to hit them exactly
      constraint_lb=targets, 
      constraint_ub=targets,
      eval_jac_g=eval_jac_g,
      eval_jac_g_structure=eval_jac_g_structure_fn(optdf),
      
      # Hessian - let's try running without the 
      # eval_h=eval_h,
      # eval_h_structure=eval_h_structure,
      
      # additional parameters passed to functions
      step2=optdf$value,
      dr=dr,
      
      opts=optlist)
  }


# use these options to test the first derivatives against ipopt finite difference approach
# optlist <- list("file_print_level"=5,
#                 "output_file"="test.out",
#                 "derivative_test"="first-order",
#                 "derivative_test_print_all"="yes",
#                 "max_iter"=20)


#****************************************************************************************************
#                Actives: final setup and run ####
#****************************************************************************************************
glimpse(steps1and2)

step2 <- steps1and2 %>% select(year, variable, step2) %>% spread(variable, step2)
step2l <- step2 %>% gather(variable, value, -year)

# problem-specific information
# set up constraints
opttarget <- step2l %>% filter(variable=="actives")
targets <- c(pvf(opttarget$value, dr), 
             pvf(opttarget$value, dr -.01), 
             pvf(opttarget$value, dr + .01), 
             opttarget$value[1])

count(step2l, variable)

# now solve for the step3 values using different guesses
optdf <- step2l %>% filter(variable=="actives.fit")
step3.fit <- callipoptr(optdf, targets, dr)
res.fit <- step2 %>% select(year) %>% mutate(step3.fit=step3.fit$solution)

optdf <- step2l %>% filter(variable=="actives.fit1")
step3.fit1 <- callipoptr(optdf, targets, dr)
res.fit1 <- step2 %>% select(year) %>% mutate(step3.fit1=step3.fit$solution)

optdf <- step2l %>% filter(variable=="actives.guess")
step3.guess <- callipoptr(optdf, targets, dr)
res.guess <- step2 %>% select(year) %>% mutate(step3.guess=step3.guess$solution)

optdf <- step2l %>% filter(variable=="actives.yy")
step3.yy <- callipoptr(optdf, targets, dr)
res.yy <- step2 %>% select(year) %>% mutate(step3.yy=step3.yy$solution)

# combine results
resall <- step2 %>% left_join(res.fit) %>%
  left_join(res.fit1) %>%
  left_join(res.guess) %>%
  left_join(res.yy)
resl <- resall %>% gather(variable, value, -year)
glimpse(resl)
count(resl, variable)

# how well do we do when we base the guess on a curve fit to known values?
resl %>% filter(variable %in% c("actives", "actives.fit", "step3.fit")) %>%
  ggplot(aes(year, value, colour=variable)) +
    geom_line()

# how well do we do when we base the guess on a perturbed version of the curve fit to known values?
resl %>% filter(variable %in% c("actives", "actives.fit1", "step3.fit1")) %>%
  ggplot(aes(year, value, colour=variable)) +
  geom_line()

# how well do we do when we base the guess on an arbitrarily constructed but plausible gaussian curve?
resl %>% filter(variable %in% c("actives", "actives.guess", "step3.guess")) %>%
  ggplot(aes(year, value, colour=variable)) +
  geom_line()
# even a pretty wild guess is pulled fairly close to the actual

# how well do we do when we base the guess on Yimeng's estimates based on benefit calculator?
resl %>% filter(variable %in% c("actives", "actives.yy", "step3.yy")) %>%
  ggplot(aes(year, value, colour=variable)) +
  geom_line()
# note that the NLP pulled the early values down toward zero; presumably the penalty function
# did not penalize those moves by very much; we should think about whether there is any easy
# way to have higher penalties for moves that change the shape of the curve so much, especially in 
# the early years

# compare actual, and optimization using as guess: fit based on actual, 2 perturbations, and yy
resl %>% filter(variable %in% c("actives", "step3.fit", "step3.fit1", "step3.guess", "step3.yy")) %>%
  ggplot(aes(year, value, colour=variable)) +
  geom_line()

# compare inputs: actual, fit based on actual, and 2 perturbations
resl %>% filter(variable %in% c("actives", "actives.fit", "actives.fit1", "actives.guess")) %>%
  ggplot(aes(year, value, colour=variable)) +
  geom_line()


#****************************************************************************************************
#                Verify that we can walk back from indexed values to initial nominal values ####
#****************************************************************************************************
glimpse(df2) # original cash flows after minor adjustments
glimpse(resall)
resadj <- resall %>% select(year, starts_with("step3")) %>%
  mutate_each(funs(. * df2$actives[1]), -year) %>%
  left_join(select(df2, year, actives)) %>%
  gather(variable, value, -year)
resadj %>% ggplot(aes(year, value, colour=variable)) + geom_line()

# check the pvs
resadj %>% group_by(variable) %>%
  arrange(year) %>%
  summarise(value.pv=pvf(value, dr))



#****************************************************************************************************
#  Now estimate cash outflows by year by type (actives, term vesteds, retirees) assuming we only know pv's for the total and not by type ####
#****************************************************************************************************



