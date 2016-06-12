# Don Boyd
# 6/7/2016

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
library("gridExtra")

library("btools") # library that I created (install from github) - maybe not needed here
library("pdata") # library that I created (install from github) - not needed here

library("ipoptr") # get from dropbox
library("nloptr")


#****************************************************************************************************
#                Useful sites and documentation ####
#****************************************************************************************************
# https://www.mail-archive.com/nlopt-discuss@ab-initio.mit.edu/msg00616.html
# error in nloptr::auglag documentation on page 5 states that the argument 'hin' to the auglag function
# defines the inequalty constraints,  hin(x) >= 0. However, this is not true. It should say
# defines the inequality constraints, hin(x) <= 0.
# http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms#Local_derivative-free_optimization
# http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
# https://books.google.com/books?id=fxzcBQAAQBAJ&pg=PA176&lpg=PA176&dq=r+alabama+auglag+variable+upper+and+lower+bounds&source=bl&ots=v_ukv9vEiT&sig=sZovvSR42TiiBvrsSElop0IbIJc&hl=en&sa=X&ved=0ahUKEwih5rWbipLNAhUUKlIKHdLqAUwQ6AEIQTAF#v=onepage&q=r%20alabama%20auglag%20variable%20upper%20and%20lower%20bounds&f=false

# cobyla bug:
# https://github.com/scipy/scipy/issues/2891
# documented bug: successful termination when constraints violated

# ?nloptr.print.options()
# ?nl.opts


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
#                Get data ####
#****************************************************************************************************
# start from scratch
df <- read_excel(paste0("./data-raw/", "ucrp_djb.xlsx"), sheet="2015UCRP_djb")
glimpse(df)
ht(df)

dra <- .0725
# check each variable for plausibility
# df %>% ggplot(aes(year, actives)) + geom_line() # good
# df %>% ggplot(aes(year, termv)) + geom_line() # first year termv not believable, lets apply year2 to year3 growth to get adjusted year1
# df %>% ggplot(aes(year, retired)) + geom_line() # strange beginning, but we'll keep it

# adjust for oddball values
df1a <- df %>% arrange(year) %>%
  mutate(year=as.integer(year), ratio=termv / lead(termv)) %>%
  mutate(termv=ifelse(year==2015, ratio[year==2016] * termv[year==2016], termv),
         total=actives + termv + retired) %>%
  select(-ratio)
head(df1a)
# df1a %>% ggplot(aes(year, actives)) + geom_line() # good
# df1a %>% ggplot(aes(year, termv)) + geom_line() # first year termv not believable, lets apply year2 to year3 growth to get adjusted year1
# df1a %>% ggplot(aes(year, retired)) + geom_line() # strange beginning, but we'll keep it

ucrpadj <- df1a %>% mutate_each(funs(ifelse(is.na(.), 0, . / 1e6)), -year) %>%
  arrange(year) %>%
  mutate(year=row_number())

glimpse(ucrpadj) # df2a is adjusted for oddballs, has no na values, and values are in $m

# NO let's try without any further scaling (beyond putting in $ millions) and see how ipopt does
# scale all to actives year1
optbase <- ucrpadj %>% select(year, actives, termv, retired, total, actives_yy=actives_shift.yy) %>%
  mutate_each(funs(. / actives[1]), termv, retired, total, actives_yy, actives) # put actives last so that we divide by untransformed actives
glimpse(optbase)

optbase %>% gather(variable, value, -year) %>% ggplot(aes(year, value, colour=variable)) + geom_line()


#****************************************************************************************************
#                Determine "known" values ####
#****************************************************************************************************
# we "know" the pv and duration for total but not for the other flow types
# we also know year1 values

pvsums <- optbase %>% gather(variable, value, -year) %>%
  group_by(variable) %>%
  arrange(year) %>%
  summarise(pv=pvf(value, dra),
            pv.m1=pvf(value, dra - .01),
            pv.p1=pvf(value, dra + .01),
            duration=duration(value, dra),
            pch.m1=pv.m1 / pv - 1,
            pch.p1=pv.p1 / pv - 1)
pvsums # we "know" the pv and duration for total but not for the other flow types

# make a list of the knowns - these (or a subset) will be our targets (constraints)
# consider: actives pch abs value will be > total, retireds will be < total --> inequality constraints
# note that I include redundant constraints for now
knowns <- list()

knowns$actives$year1 <- optbase$actives[1]
knowns$actives$pv <- pvsums$pv[pvsums$variable=="actives"]
knowns$actives$pv.m1lb <- knowns$actives$pv * (1 + pvsums$pch.m1[pvsums$variable=="total"]) # actives pch must be at least total
knowns$actives$pv.p1ub <- knowns$actives$pv * (1 + pvsums$pch.p1[pvsums$variable=="total"])

knowns$retired$year1 <- optbase$retired[1]
knowns$retired$pv <- pvsums$pv[pvsums$variable=="retired"]
knowns$retired$pv.m1ub <- knowns$retired$pv * (1 + pvsums$pch.m1[pvsums$variable=="total"]) # actives pch must be at least total
knowns$retired$pv.p1lb <- knowns$retired$pv * (1 + pvsums$pch.p1[pvsums$variable=="total"])

knowns$termv$year1 <- optbase$termv[1]
knowns$termv$pv <- pvsums$pv[pvsums$variable=="termv"]

knowns$total$year1 <- with(knowns, actives$year1 + retired$year1 + termv$year1)
knowns$total$pv <- with(knowns, actives$pv + retired$pv + termv$pv)
knowns$total$pv.m1 <-  pvsums$pv.m1[pvsums$variable=="total"]
knowns$total$pv.p1 <-  pvsums$pv.p1[pvsums$variable=="total"]

knowns

optbase %>% gather(variable, value, -year)  %>%
  group_by(variable) %>%
  arrange(year) %>%
  summarise(maxyear=which.max(value), 
            maxyval=max(value), 
            year1val=first(value), 
            ratio=maxyval / first(value)) %>%
  kable(digits=1)


#****************************************************************************************************
#                Construct parameter list (plist) ####
#****************************************************************************************************

plist <- list()
plist <- within(plist, {
  years <- 1:105
  dr <- 0.0725
  targnames <- c("year1", "pv", "pv.m1", "pv.p1")
})

plist <- within(plist, {
  actives <- list()
  actives <- within(actives, {
    # year1, pv, pv.m1, pv.p1
    # set targets up as vectors - makes them easier to work with but harder to follow
    targ.lbtols <- c(0.001, 0.001, 0, 0.5) # m1 can't be lower than what we will estimate, p1 can be as low as needed
    targ.ubtols <- c(0.001, 0.001, 0.5, 0) # m1 can be as high as needed, p1 can't be more than what we estimate
    names(targ.lbtols) <- targnames
    names(targ.ubtols) <- targnames
    targs <- c(knowns$actives$year1, 
               knowns$actives$pv,
               # actives pv.m1 will be GREATER than with pch from total
               knowns$actives$pv * (1 + pvsums$pch.m1[pvsums$variable=="total"]),
               # actives pv.p1 will be LESS than with pch from total
               knowns$actives$pv * (1 + pvsums$pch.p1[pvsums$variable=="total"]))
    names(targs) <- targnames
    targs.lb <- targs * (1 - targ.lbtols)
    names(targs.lb) <- targnames
    targs.ub <- targs * (1 + targ.ubtols)
    names(targs.ub) <- targnames
  })
})

plist <- within(plist, {
  retired <- list()
  retired <- within(retired, {
    # year1, pv, pv.m1, pv.p1
    # set targets up as vectors - makes them easier to work with but harder to follow
    targ.lbtols <- c(0.001, .001, 0.5, 0) # m1, p1
    targ.ubtols <- c(0.001, .001, 0, 0.5) # m1, p1
    names(targ.lbtols) <- targnames
    names(targ.ubtols) <- targnames
    targs <- c(knowns$retired$year1, 
               knowns$retired$pv,
               # retired pv.m1 will be LESS than with pch from total
               knowns$retired$pv * (1 + pvsums$pch.m1[pvsums$variable=="total"]),
               # retired pv.p1 will be GREATER than with pch from total
               knowns$retired$pv * (1 + pvsums$pch.p1[pvsums$variable=="total"]))
    names(targs) <- targnames
    targs.lb <- targs * (1 - targ.lbtols)
    names(targs.lb) <- targnames
    targs.ub <- targs * (1 + targ.ubtols)
    names(targs.ub) <- targnames
  })
})

plist <- within(plist, {
  termv <- list()
  termv <- within(termv, {
    # year1, pv, pv.m1, pv.p1
    # set targets up as vectors - makes them easier to work with but harder to follow
    # we don't really know what the m1 p1 values should be like so allow a wide range
    targ.lbtols <- c(0.001, .001, 0.5, 0.5) # m1, p1
    targ.ubtols <- c(0.001, .001, 0.5, 0.5) # m1, p1
    names(targ.lbtols) <- targnames
    names(targ.ubtols) <- targnames
    targs <- c(knowns$termv$year1, 
               knowns$termv$pv,
               knowns$termv$pv * (1 + pvsums$pch.m1[pvsums$variable=="total"]),
               knowns$termv$pv * (1 + pvsums$pch.p1[pvsums$variable=="total"]))
    names(targs) <- targnames
    targs.lb <- targs * (1 - targ.lbtols)
    names(targs.lb) <- targnames
    targs.ub <- targs * (1 + targ.ubtols)
    names(targs.ub) <- targnames
  })
})
  
plist <- within(plist, {
  total <- list()
  total <- within(total, {
    # year1, pv, pv.m1, pv.p1
    # set targets up as vectors - makes them easier to work with but harder to follow
    targ.lbtols <- c(0.0, .001, .001, .001)
    targ.ubtols <- c(0.0, .001, .001, .001)
    names(targ.lbtols) <- targnames
    names(targ.ubtols) <- targnames
    targs <- c(plist$actives$targs["year1"] + plist$retired$targs["year1"] + plist$termv$targs["year1"],
               plist$actives$targs["pv"] + plist$retired$targs["pv"] + plist$termv$targs["pv"],
               pvsums$pv.m1[pvsums$variable=="total"], 
               pvsums$pv.p1[pvsums$variable=="total"])
    names(targs) <- targnames
    targs.lb <- targs * (1 - targ.lbtols)
    names(targs.lb) <- targnames
    targs.ub <- targs * (1 + targ.ubtols)
    names(targs.ub) <- targnames
  })
})
  
plist
pvsums
knowns


#****************************************************************************************************
#                Optimization functions ####
#****************************************************************************************************
gauss <- function(x, years) {
  # x[1:3]: mu, sigma, k
  x[3] * exp(-0.5 * ((years - x[1]) / x[2])^2)
}

eqfn_all <- function(x, plist=plist){
  iactives <- 1:3
  iretired <- 4:6
  itermv <- 7:9
  # x is concatenation of parameters for actives, retired, termv, so it has 9 elements
  # x.flow is a subset of x, for a given flow (actives, retired, termv), and so has 3 elements
  
  gauss <- function(x.flow, years) {x.flow[3] * exp(-0.5 * ((years - x.flow[1]) / x.flow[2])^2)}
  
  yvals.actives <- gauss(x[iactives], plist$years)
  pv.actives <- sum(yvals.actives / ((1 + plist$dr)^plist$years))
  
  yvals.retired <- gauss(x[iretired], plist$years)
  pv.retired <- sum(yvals.retired / ((1 + plist$dr)^plist$years))
  
  yvals.termv <- gauss(x[itermv], plist$years)
  pv.termv <- sum(yvals.termv / ((1 + plist$dr)^plist$years))
  
  pv <- pv.actives + pv.retired + pv.termv
  
  obj <- (pv - plist$total$targs["pv"])^2
  return(obj)
}


hineq_all <- function(x, ...){
  iactives <- 1:3
  iretired <- 4:6
  itermv <- 7:9
  # x is concatenation of parameters for actives, retired, termv, so it has 9 elements
  # x.flow is a subset of x, for a given flow (actives, retired, termv), and so has 3 elements
  gauss <- function(x.flow, years) {x.flow[3] * exp(-0.5 * ((years - x.flow[1]) / x.flow[2])^2)}
  
  pvf <- function(vec, dr) {
    pv <- vec / ((1+dr)^((1:length(vec))-1))
    pv <- sum(pv, na.rm=TRUE)
    return(pv)
  }
  
  evalcons <- function(x.flow, plist) {
    cval <- numeric(4)
    names(cval) <- c("year1", "pv", "pv.m1", "pv.p1")
    cval[1] <- gauss(x.flow, 1)
    cval[2] <- pvf(gauss(x.flow, plist$years), plist$dr)
    cval[3] <- pvf(gauss(x.flow, plist$years), plist$dr - .01)
    cval[4] <- pvf(gauss(x.flow, plist$years), plist$dr + .01)
    return(cval)
  }
  
  convals.actives <- evalcons(x[iactives], plist)
  convals.retired <- evalcons(x[iretired], plist)
  convals.termv <- evalcons(x[itermv], plist)
  
  h <- numeric(30)
  
  # actives constraints
  h[1:4] <- convals.actives - plist$actives$targs.lb[names(convals.actives)]
  h[5:8] <- plist$actives$targs.ub[names(convals.actives)] - convals.actives
  
  # retired constraints
  h[9:12] <- convals.retired - plist$retired$targs.lb[names(convals.retired)]
  h[13:16] <- plist$retired$targs.ub[names(convals.retired)] - convals.retired
  
  # termv constraints
  h[17:20] <- convals.termv - plist$termv$targs.lb[names(convals.termv)]
  h[21:24] <- plist$termv$targs.ub[names(convals.termv)] - convals.termv
  
  # total constraints
  tot.pv.year1 <- convals.actives["year1"] + convals.retired["year1"] + convals.termv["year1"]
  h[25] <- tot.pv.year1 - plist$total$targs.lb["year1"]
  h[26] <- plist$total$targs.ub["year1"] - tot.pv.year1
  
  # total constraints (duration)
  tot.pv.m1 <- convals.actives["pv.m1"] + convals.retired["pv.m1"] + convals.termv["pv.m1"]
  h[27] <- tot.pv.m1 - plist$total$targs.lb["pv.m1"]
  h[28] <- plist$total$targs.ub["pv.m1"] - tot.pv.m1
  
  tot.pv.p1 <- convals.actives["pv.p1"] + convals.retired["pv.p1"] + convals.termv["pv.p1"]
  h[29] <- tot.pv.p1 - plist$total$targs.lb["pv.p1"]
  h[30] <- plist$total$targs.ub["pv.p1"] - tot.pv.p1
  
  return(h)
}


#****************************************************************************************************
#                Call optimization routine ####
#****************************************************************************************************

# we stack all vectors as actives, retired, termv
# each row below has 3 variables: mu, sigma, k
# completely unconstrained parameters
x0_all <- c(30, 20, 10,
            30, 20, 10,
            30, 20, 10)
xlb_all <- c(0, 0, 0,
             0, 0, 0,
             0, 0, 0)
xub_all <- c(50, 50, 50,
             50, 50, 50,
             50, 50, 50)
# using a priori knowledge / belief:
# x0_all <- c(36, 15, 10,
#             30, 20, 10,
#             30, 20, 10)
# xlb_all <- c(34, 0, 0,
#              0, 0, 0,
#              0, 0, 0)
# xub_all <- c(50, 15.5, 50,
#              50, 50, 50,
#              50, 50, 50)

(x0_all >= xlb_all) & (x0_all <= xub_all)

opt <- slsqp(x0_all, eqfn_all, lower=xlb_all, upper=xub_all, hin=hineq_all, control=list(xtol_rel = 1e-8), plist=plist)
str(opt)


#****************************************************************************************************
#                Examine results ####
#****************************************************************************************************

fitdf <- data_frame(year=plist$years,
                   actives.fit=gauss(opt$par[1:3], year),
                   retired.fit=gauss(opt$par[4:6], year),
                   termv.fit=gauss(opt$par[7:9], year),
                   total.fit=actives.fit + retired.fit + termv.fit)
glimpse(fitdf)
knowns


# create table comparing results to targets and other knowns
knowndf <- pvsums %>% filter(variable!="actives_yy") %>%
  rename(flow=variable) %>%
  mutate(valtype="known") %>%
  left_join(filter(optbase, year==1) %>% 
              select(actives, termv, retired, total) %>%
              mutate(valtype="known") %>%
              gather(flow, year1, -valtype)) %>%
  select(flow, valtype, year1, everything())

estdf <- fitdf %>% gather(flow, value, -year) %>%
  group_by(flow) %>%
  arrange(year) %>%
  summarise(pv=pvf(value, dra),
            pv.m1=pvf(value, dra - .01),
            pv.p1=pvf(value, dra + .01),
            duration=duration(value, dra),
            pch.m1=pv.m1 / pv - 1,
            pch.p1=pv.p1 / pv - 1) %>%
  mutate(valtype="estimate") %>%
  left_join(filter(fitdf, year==1) %>% select(-year) %>% gather(flow, year1)) %>%
  select(flow, valtype, year1, everything()) %>%
  mutate(flow=str_replace(flow, ".fit", ""))

compdf <- bind_rows(knowndf, estdf) %>%
  gather(variable, value, -flow, -valtype) %>%
  spread(valtype, value) %>%
  mutate(est_minus_known=estimate - known,
         pctdiff=est_minus_known / known * 100,
         targeted=ifelse(variable %in% c("duration") & flow!="total", "no", ""),
         targeted=ifelse(flow=="actives" & variable %in% c("pch.m1", "pv.m1"), "lower_bound", targeted),
         targeted=ifelse(flow=="actives" & variable %in% c("pch.p1", "pv.p1"), "upper_bound", targeted),
         targeted=ifelse(flow=="retired" & variable %in% c("pch.m1", "pv.m1"), "upper_bound", targeted),
         targeted=ifelse(flow=="retired" & variable %in% c("pch.p1", "pv.p1"), "lower_bound", targeted),
         targeted=ifelse(flow=="termv" & !variable %in% c("year1", "pv"), "no", targeted))
compdf %>% kable(digits=1)


# graph results
tmp <- optbase %>% select(-actives_yy)
vnames <- paste0(names(tmp), ".dat")
names(tmp)[2:5] <- vnames[2:5]
p1 <- tmp %>% gather(variable, value, -year) %>% 
  ggplot(aes(year, value, colour=variable)) + 
    geom_line() + 
    scale_y_continuous(breaks=seq(0, 60, 2), limits=c(0, 20)) +
    scale_x_continuous(breaks=seq(0, 200, 10)) + 
    ggtitle("Original data")

p2 <- fitdf %>% gather(variable, value, -year) %>% 
  ggplot(aes(year, value, colour=variable)) + 
    geom_line() + 
    scale_y_continuous(breaks=seq(0, 60, 2), limits=c(0, 20)) +
    scale_x_continuous(breaks=seq(0, 200, 10)) + 
    ggtitle("Estimated data")

p3 <- left_join(select(tmp, year, total.dat),
                select(fitdf, year, total.fit)) %>%
  gather(variable, value, -year) %>%
  ggplot(aes(year, value, colour=variable)) + 
    geom_line() + 
    scale_y_continuous(breaks=seq(0, 60, 2), limits=c(0, 20)) +
    scale_x_continuous(breaks=seq(0, 200, 10)) + 
    ggtitle("Original vs. estimated total")

grid.arrange(p1, p2, p2, p3, ncol=2)


