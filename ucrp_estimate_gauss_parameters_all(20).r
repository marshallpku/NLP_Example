# Don Boyd
# 6/8/2016

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
#                Get ucrp data so that we can develop "known" values and get cash flows ####
#****************************************************************************************************
# start from scratch
df <- read_excel(paste0("./data-raw/", "ucrp_djb.xlsx"), sheet="2015UCRP_djb")
glimpse(df)
ht(df)


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

glimpse(ucrpadj) # ucrpadj is adjusted for oddballs, has no na values, and values are in $m

# scale all to actives year1
optbase <- ucrpadj %>% select(year, actives, termv, retired, total, actives_yy=actives_shift.yy) %>%
  mutate_each(funs(. / actives[1]), termv, retired, total, actives_yy, actives) # put actives last so that we divide by untransformed actives
glimpse(optbase)

optbase %>% gather(variable, value, -year) %>% ggplot(aes(year, value, colour=variable)) + geom_line()


#****************************************************************************************************
#                Determine "known" values for ucrp ####
#****************************************************************************************************
# we "know" the pv and duration for total but not for the other flow types
# we also know year1 values
dra <- .0725
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

optbase %>% gather(variable, value, -year)  %>%
  group_by(variable) %>%
  arrange(year) %>%
  summarise(maxyear=which.max(value), 
            maxyval=max(value), 
            year1val=first(value), 
            ratio=maxyval / first(value)) %>%
  kable(digits=1)

# **** Here is where we begin with data for real or constructed plans **** ########

#****************************************************************************************************
#                Read in raw data and make adjusted df ####
#****************************************************************************************************
plans.s <- 
"planid, planname, dr, pv_actives, pv_retired, pv_termv, pv_total, pv.m1_total, pv.p1_total, year1_actives, year1_retired, year1_termv, year1_total
1, ucrp, .0725, 88.030372, 49.740722, 7.840111, 145.611205, 171.043708, 125.793767, 1, 4.181388, 0.2862846, 5.467673"

plans.df <- read_csv(plans.s)
plans.df

as.list(plans.df[1, ])

theplanid <- 1 # the plan to examine
the_plan <- plans.df %>% filter(planid==theplanid)


#****************************************************************************************************
#                Set up targets (constraints) ####
#****************************************************************************************************

# create equality constraints for some targets and inequality bounds for others
# if we find it is too hard to hit targets exactly we can change each equality target
# to a set of 2 inequality targets - a lower bound and upper bound

targets.eq <- the_plan %>% select(-dr) %>%
  gather(variable, target, -planid, -planname) %>%
  separate(variable, c("targtype", "flowtype"), sep="_") %>%
  select(planid, planname, flowtype, targtype, everything()) %>%
  arrange(planid, planname, flowtype, targtype)
# keep the pv.m1 and pv.p1 values for total flows in here, just in case we want to target them as equalities
# however, as implemented, right now we target as inequalities with a small tolerance around them

totpvtol <- .005 # tolerance for establishing bounds on pv.m1_total and p1
targets.ineq <- the_plan %>% select(planid, planname, pv_actives, pv_retired, pv_total, pv.m1_total, pv.p1_total) %>%
  mutate(pv.m1_lower_actives=pv_actives * pv.m1_total / pv_total,
         pv.p1_upper_actives=pv_actives * pv.p1_total / pv_total,
         pv.m1_upper_retired=pv_retired * pv.m1_total / pv_total,
         pv.p1_lower_retired=pv_retired * pv.p1_total / pv_total,
         pv.m1_lower_total=pv.m1_total * (1 - totpvtol),
         pv.m1_upper_total=pv.m1_total * (1 + totpvtol),
         pv.p1_lower_total=pv.p1_total * (1 - totpvtol),
         pv.p1_upper_total=pv.p1_total * (1 + totpvtol)) %>%
  select(-starts_with("pv_"), -contains("1_total")) %>%
  gather(variable, bound, -planid, -planname) %>%
  separate(variable, c("targtype", "boundtype", "flowtype"), sep="_") %>%
  select(planid, planname, flowtype, targtype, boundtype, bound) %>%
  mutate(mult=ifelse(boundtype=="upper", -1, 1)) %>% # multiplier to use for: mult * (constraintval - bound) >= 0
  arrange(planid, planname, flowtype, targtype, boundtype)

targets.eq
targets.ineq


#****************************************************************************************************
#                Construct parameter list (plist) ####
#****************************************************************************************************

plist <- list()
plist$planid <- the_plan$planid
plist$planname <- the_plan$planname
plist$dr <- the_plan$dr
plist$years <- 1:105
plist$iactives <- 1:3 # indexes into the 9 element array of gauss function parameters
plist$iretired <- 4:6
plist$itermv <- 7:9
plist$pvf <- function(vec, dr) {pv <- vec / ((1+dr)^((1:length(vec))-1)); pv <- sum(pv, na.rm=TRUE); return(pv)}
plist$gauss <- function(pars, years) {pars[3] * exp(-0.5 * ((years - pars[1]) / pars[2])^2)} # pars[1:3]: mu, sigma, k

# convert the equality targets to a matrix with names, for ease of use
targs <- targets.eq %>% select(flowtype, targtype, target) %>% spread(targtype, target)
targets.equal.mat <- targs[, -1] %>% as.matrix()
rownames(targets.equal.mat) <- targs$flowtype
plist$targets.equal <- targets.equal.mat # we keep all totals in here even though we don't target p1, m1 as equalities

# get lower and upper bound matrices for inequality targets
targslower <- targets.ineq %>% filter(boundtype=="lower") %>%
  select(flowtype, targtype, bound) %>% 
  spread(targtype, bound)
lower.mat <- targslower[, -1] %>% as.matrix()
rownames(lower.mat) <- targslower$flowtype
plist$targets.lower <- lower.mat 

targsupper <- targets.ineq %>% filter(boundtype=="upper") %>%
  select(flowtype, targtype, bound) %>% 
  spread(targtype, bound)
upper.mat <- targsupper[, -1] %>% as.matrix()
rownames(upper.mat) <- targsupper$flowtype
plist$targets.upper <- upper.mat 

plist


#****************************************************************************************************
#                Optimization functions ####
#****************************************************************************************************

objfn <- function(x, plist=plist){
  # x is concatenation of parameters for actives, retired, termv, so it has 9 elements
  yvals.actives <- plist$gauss(x[plist$iactives], plist$years)
  yvals.retired <- plist$gauss(x[plist$iretired], plist$years)
  yvals.termv <- plist$gauss(x[plist$itermv], plist$years)
  yvals.total <- yvals.actives + yvals.retired + yvals.termv
  
  pv_total <- plist$pvf(yvals.total, plist$dr)
  obj <- (pv_total - plist$targets.eq["total", "pv"])^2
  
  # adding p1 and m1 to the obj function, even downweighted, does not work very well
  # pv.p1_total <- plist$pvf(yvals.total, plist$dr + .01)
  # pv.m1_total <- plist$pvf(yvals.total, plist$dr - .01)
  # 
  # obj <- (pv_total - plist$targets.eq["total", "pv"])^2 +
  #   (pv.p1_total - plist$targets.eq["total", "pv.p1"])^2 / 1e3 +
  #   (pv.m1_total - plist$targets.eq["total", "pv.m1"])^2 / 1e3
    
  return(obj)
}


heq <- function(x, ...) {
  # define equality constraints
  getyear1 <- function(x, flowtype, plist) {
    flowparam.indexes <- plist[[paste0("i", flowtype)]]
    year1.cons <- plist$gauss(x[flowparam.indexes], 1)
    return(year1.cons)
  }
  
  getpv <- function(x, flowtype, plist) {
    flowparam.indexes <- plist[[paste0("i", flowtype)]]
    yvals <- plist$gauss(x[flowparam.indexes], plist$years)
    pv.cons <- pvf(yvals, plist$dr)
    return(pv.cons)
  }
  
  cons.eq <- numeric(6)
  cons.eq[1] <- getyear1(x, "actives", plist) - plist$targets.equal["actives", "year1"]
  cons.eq[2] <- getyear1(x, "retired", plist) - plist$targets.equal["retired", "year1"]
  cons.eq[3] <- getyear1(x, "termv", plist) - plist$targets.equal["termv", "year1"]
  
  cons.eq[4] <- getpv(x, "actives", plist) - plist$targets.equal["actives", "pv"]
  cons.eq[5] <- getpv(x, "retired", plist) - plist$targets.equal["retired", "pv"]
  cons.eq[6] <- getpv(x, "termv", plist) - plist$targets.equal["termv", "pv"]
   
  return(cons.eq)
}


hin <- function(x, ...) {
  # define inequality constraints
  actives <- plist$gauss(x[plist$iactives], plist$years)
  actives.m1 <- pvf(actives, plist$dr - .01)
  actives.p1 <- pvf(actives, plist$dr + .01)
  
  retired <- plist$gauss(x[plist$iretired], plist$years)
  retired.m1 <- pvf(retired, plist$dr - .01)
  retired.p1 <- pvf(retired, plist$dr + .01)
    
  termv <- plist$gauss(x[plist$itermv], plist$years)
  termv.m1 <- pvf(termv, plist$dr - .01)
  termv.p1 <- pvf(termv, plist$dr + .01)
    
  total.pv.m1 <- actives.m1 + retired.m1 + termv.m1
  total.pv.p1 <- actives.p1 + retired.p1 + termv.p1

  cons.ineq <- numeric(8)
  cons.ineq[1] <- total.pv.m1  - plist$targets.lower["total", "pv.m1"]
  cons.ineq[2] <- plist$targets.upper["total", "pv.m1"] - total.pv.m1
  
  cons.ineq[3] <- total.pv.p1  - plist$targets.lower["total", "pv.p1"]
  cons.ineq[4] <- plist$targets.upper["total", "pv.p1"] - total.pv.p1
  
  # now put in one-sided bounds - these don't seem to have much impact
  cons.ineq[5] <- actives.m1  - plist$targets.lower["actives", "pv.m1"]
  cons.ineq[6] <- retired.p1  - plist$targets.lower["retired", "pv.p1"]
  
  cons.ineq[7] <- plist$targets.upper["actives", "pv.p1"] - actives.p1
  cons.ineq[8] <- plist$targets.upper["retired", "pv.m1"] - retired.m1

  return(cons.ineq)
}


#****************************************************************************************************
#                Call optimization routine ####
#****************************************************************************************************
# parameter starting values, lower bounds, and upper bounds
# mu, sigma, and k
x0 <- c(30, 20, 10, # actives
        30, 20, 10, # retired
        30, 20, 10) # termv
xlb <- c(0, 0, 0,
         0, 0, 0,
         0, 0, 0)
xub <- c(50, 50, 50,
         50, 50, 50,
         50, 50, 50)
(x0 >= xlb) & (x0 <= xub)

objfn(x0, plist)
heq(x0, plist)
hin(x0, plist)


opt <- slsqp(x0, objfn, lower=xlb, upper=xub, heq=heq, hin=hin, control=list(xtol_rel = 1e-6), plist=plist)
str(opt)
opt$par %>% round(2)


#****************************************************************************************************
#                Examine results ####
#****************************************************************************************************

fitdf <- data_frame(year=plist$years,
                    actives.fit=plist$gauss(opt$par[1:3], year),
                    retired.fit=plist$gauss(opt$par[4:6], year),
                    termv.fit=plist$gauss(opt$par[7:9], year),
                    total.fit=actives.fit + retired.fit + termv.fit)
glimpse(fitdf)

fitdf %>% gather(variable, value, -year) %>% 
  ggplot(aes(year, value, colour=variable)) + 
  geom_line() + 
  scale_y_continuous(breaks=seq(0, 60, 2), limits=c(0, 20)) +
  scale_x_continuous(breaks=seq(0, 200, 10)) + 
  ggtitle("Estimated data")


# create table comparing results to targets and other knowns
knowndf <- pvsums %>% filter(variable!="actives_yy") %>%
  rename(flow=variable) %>%
  mutate(valtype="known", pch.m1=pch.m1 * 100, pch.p1=pch.p1 * 100) %>%
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
            pch.m1=pv.m1 / pv * 100 - 100,
            pch.p1=pv.p1 / pv * 100 - 100) %>%
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
         targeted=ifelse(flow=="termv" & !variable %in% c("year1", "pv"), "no", targeted),
         targeted=ifelse(flow=="total" & !variable %in% c("pv.m1", "pv.p1"), "bounds", targeted),
         targeted=ifelse(flow=="total" & !variable %in% c("pch.m1", "pch.p1"), "indirect", targeted))
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


