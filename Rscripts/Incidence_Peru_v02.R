rm(list = ls()) #only delete stuff in the visible environment; the best approach is to restart the session

# load packages
library(readxl)
library(data.table)
library(ggplot2)
library(splines)

## 1.2 Functions for prevalence data ----
logit  <- function(x){log(x/(1-x)) }
expit  <- function(x){1/(1+exp(-x))}

# List of years
years <- 2016:2023

# Initialize an empty data.table
dt.endes_all <- data.table()

# Loop through each year and bind to endes_all
for (year in years) {
  # Construct the file path
  file_path <- c("data/Estimated.prev.2016-2023.xlsx")

  # Read and process the data
  dt.prev <- data.table(read_xlsx(file_path, sheet = as.character(year)))
  setnames(dt.prev, "median prev", "prev")
  dt.prev[, `:=`(
    age = 2.5 + as.numeric(substr(BAND, 1, 2)),
    year = year+0.5,
    prev = 0.01 * prev,
    SEX = fcase(
      SEX == 1, "Men",
      SEX == 2, "Women"
    )
    )]

  dt.prev[, `:=`(
    prev_logit = logit(prev),
    low_logit = logit(low/100))]
  dt.prev[ , SE_logit := (prev_logit - low_logit) / 1.96]

  # Append to the full dataset
  dt.endes_all <- rbind(dt.endes_all, dt.prev, fill = TRUE)
}


### function returns time, age and sex-specific prevalence
#pmod <- lm(prev_logit ~ ns(age, df = 3) * year * SEX, data = dt.endes_all)

# Best fit so far:
#pmod <- lm(prev_logit ~ ns(age, df = 3) + year + SEX, data = dt.endes_all)
# do interactions improve fit?
pmod1 <- lm(prev_logit ~ ns(age, df = 3) + year + SEX, data = dt.endes_all)
pmod2 <- lm(prev_logit ~ ns(age, df = 3) + year * SEX, data = dt.endes_all)
pmod3 <- lm(prev_logit ~ ns(age, df = 3) * year + SEX, data = dt.endes_all)
pmod4 <- lm(prev_logit ~ ns(age, df = 3) * year * SEX, data = dt.endes_all)
pmod5 <- lm(prev_logit ~ ns(age, df = 3) + ns(year, df = 2) + SEX, data = dt.endes_all)
pmod6 <- lm(prev_logit ~ ns(age, df = 3) + ns(year, df = 2) * SEX, data = dt.endes_all)
pmod7 <- lm(prev_logit ~ ns(age, df = 3) * ns(year, df = 2) + SEX, data = dt.endes_all)
pmod8 <- lm(prev_logit ~ ns(age, df = 3) * ns(year, df = 2) * SEX, data = dt.endes_all)
# maybe slightly for year*sex when considering R and AIC
BIC(pmod1, pmod2, pmod3, pmod4, pmod5, pmod6, pmod7, pmod8)
pmod <- lm(prev_logit ~ ns(age, df = 3) * year + SEX, data = dt.endes_all)
summary(pmod)

dt.endes_all$prev_logit.predicted <- predict(pmod)

plt.input_prev <- ggplot(data = dt.endes_all, aes(x = age, y = expit(prev_logit.predicted) * 100)) +
  geom_line() +
  geom_point(aes(y = prev*100)) +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  ylab("Prevalence") +
  xlab("Age in years") +
  theme_bw() +
  theme(strip.text = element_text(size = 13), axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_y_continuous(breaks= seq(0, 35, 5)) +
  facet_grid(rows = vars(SEX), cols = vars(year))

plt.input_prev2 <- ggplot(data = dt.endes_all, aes(x = year, y = expit(prev_logit.predicted) * 100)) +
  geom_line() +
  geom_point(aes(y = prev*100)) +
  geom_errorbar(aes(ymin = low, ymax = high)) +
  ylab("Prevalence") +
  xlab("year") +
  theme_bw() +
  theme(strip.text = element_text(size = 13), axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_y_continuous(breaks= seq(0, 35, 5)) +
  facet_grid(rows = vars(SEX), cols = vars(BAND))


get_prev <- function(t, a, p.Sex, rnd = 0){
  ln.p  <- predict(pmod, newdata = data.frame(year = t, age = a, SEX = p.Sex))

  p.SE_ <- dt.endes_all[year %in% t & SEX == p.Sex, SE_logit]
  a_    <- dt.endes_all[year %in% t & SEX == p.Sex, age]
  p.SE  <- approx(a_, p.SE_, xout = a, rule = 2)$y

  return(pmax(0, expit(ln.p + rnd * p.SE)))
}

plot(20:80, get_prev(2023.5, 20:80, p.Sex = "Men", rnd = 0), ylim = c(0, 0.25))
lines(20:80, get_prev(2023.5, 20:80, p.Sex = "Men", rnd = 1.96), ylim = c(0, 0.25))
lines(20:80, get_prev(2023.5, 20:80, p.Sex = "Men", rnd = -1.96), ylim = c(0, 0.25))

### function returns change prevalence between survey years
get_dp <- function(t1, t2, dp.a, dp.Sex, rnd.p.1 = 0, rnd.p.2 = 0){
  dT  <- t2 - t1
  p1  <- get_prev(t1, dp.a - dT / 2, p.Sex = dp.Sex, rnd = rnd.p.1)
  p2  <- get_prev(t2, dp.a + dT / 2, p.Sex = dp.Sex, rnd = rnd.p.2)
  dp  <- (p2 - p1) / dT
  return(dp)
}

### function returns  prevalence midway between survey years
get_p.mid <- function(t1, t2, p.a, dp.Sex, rnd.p.1 = 0, rnd.p.2 = 0){
  dT  <- t2 - t1
  p1  <- get_prev(t1, p.a - dT / 2, p.Sex = dp.Sex, rnd = rnd.p.1)
  p2  <- get_prev(t2, p.a + dT / 2, p.Sex = dp.Sex, rnd = rnd.p.2)
  pt  <- 0.5 * (p2 + p1)
  return(pt)
}

## 1.3 Data and function for mortality rate ratio (diab vs. no diab) ----
# MRR function as in India paper
get_R <- function(a_){
  R30 <- runif(1, min = 7.5, max = 8.5)
  R90 <- runif(1, min = 2, max = 2.5)
  logR_ <- approx(c(30, 90), log(c(R30, R90)), xout = a_, rule = 2)$y
  return(exp(logR_))
}
#r_vec <- rep(NA, 1000000)
#for(i in 1:1000000){
#  r_i_ <- get_R(80)
#  r_vec[i] <- r_i_
#}
#range(r_vec)
# range at age 80: 2.493110 3.065325

### mortality rate ratio (MRR) from  Alegre-Díaz et al. 2016
age.MRR <- c(40, 50, 60, 70)
MRR.input <- c(7.5, 4.8, 3.3, 2.7)

### model relationship between age and MRR
lm.MRR <- lm(log(MRR.input) ~ ns(age.MRR, df = 2))

### predicted age-specific MRR based on the model
cbind(c(20:100), exp(predict(lm.MRR, newdata = data.frame(age.MRR = 20:100))))

### compare with input data
plot(20:100, predict(lm.MRR, newdata = data.frame(age.MRR = 20:100)), type = "l")
points(age.MRR, log(MRR.input), pch = 3) # points for input MRR

### functions returns age-specific MRR from the model above
fct.MRR.Mex <- function(a){
  MRR.mod_ <- exp(predict(lm.MRR, newdata = data.frame(age.MRR = a)))
  return(pmax(MRR.mod_, 1.1))
}
plot(20:100, fct.MRR.Mex(20:100))
fct.MRR.Mex(20)


## 1.4 Data and function for general mortality ----
# general mortality from here: https://population.un.org/dataportal/data/indicators/80/locations/604/start/2010/end/2023/table/pivotbysex?df=99779a07-5b67-43bb-93aa-b0b03b82a778
gen_mort <- fread("data/unpopulation_dataportal_20250428135409.csv")
gen_mort[ , `:=`(
  Age = as.numeric(substr(Age, 1, 3)),
  Sex = fcase(
    Sex == "Male", "Men",
    Sex == "Female", "Women")
  )]
setorder(gen_mort, Time, Sex, Age)
gen_mort <- gen_mort[, .(Time, Age, Sex, Value)]

fct.gen_mort <- function(t, a, sex_mort) {
  mort_dat <- gen_mort[Time == t & Sex == sex_mort]
  mort_    <- approx(mort_dat$Age, mort_dat$Value, xout = a, rule = 2)$y
  return(mort_)
}

plt.mortality <- ggplot(data = gen_mort, aes(x = Age, y = log(Value), colour = Sex)) +
  geom_line(linewidth = 1) +
  #geom_point() +
  facet_wrap(vars(Time), ncol = 2)

plt.heat <- ggplot(data = gen_mort, aes(x = Time, y = Age, fill = log(Value))) +
  geom_tile() +
  scale_fill_viridis_c()


# 2. Estimate incidence ----

### create matrix to store the results
out.ages <- 15:90 # age range in years
out.years <- seq(2017, 2023, 1) # year of incidence estimation
svy.years <- seq(2016.5, 2023.5, 1) # ENDES survey years
nboots <- 1000 # no. of repetitions

i.res <- array(dim = c(length(out.ages), length(svy.years) - 1, 2, nboots),
               dimnames = list(c(out.ages), c(out.years), c("Men", "Women"), 1:nboots))

### apply partial differential equation (PDE) to input data
set.seed(68763)
for(bb in 1:nboots){ # nboots loops for random error estimation
  for(ii in c("Men", "Women")){ # repeat analysis for men and women separately
    for(tt in 1:(length(svy.years) - 1)){ # and midpoints between surveys
      #rnd.vec <- c(0,0)
      rnd.vec <- rnorm(2) # draw random error
      MRR_  <- fct.MRR.Mex(out.ages) # get MRR
      #MRR_  <- get_R(out.ages) # get MRR
      t1_   <- svy.years[tt] # define t1
      t2_   <- svy.years[tt + 1] # define t2
      t_    <- 0.5 * (t1_ + t2_) # midpoint between t1 and t2

      # get general mortality at midpoint between t1 and t2
      m_    <- fct.gen_mort(t = round(t_), a = out.ages, sex_mort = ii)

      # get prevalence at midpoint between t1 and t2
      p_    <- get_p.mid(t1 = t1_, t2 = t2_, p.a = out.ages, dp.Sex = ii,
                         rnd.p.1 = rnd.vec[1], rnd.p.2 = rnd.vec[2])

      # get change in prevalence at midpoint between t1 and t2
      dp_   <- get_dp(t1 = t1_, t2 = t2_, dp.a = out.ages, dp.Sex = ii,
                      rnd.p.1 = rnd.vec[1], rnd.p.2 = rnd.vec[2])

      # calculate incidence at midpoint between t1 and t2 using the PDE
      i_  <- dp_/(1-p_) + m_ * p_ * (MRR_ - 1) / (1 + p_ * (MRR_ - 1))

      # put incidence into results matrix
      i.res[ , as.character(t_), ii, bb] <- pmax(i_, 0)
    }
  }
}

### matrix as data.table
dt.incidence <- as.data.table(i.res, value.name = "incidence")
setnames(dt.incidence, 1:4, c("age", "Year", "Sex", "boot"))

### change type of age and year to numeric
dt.incidence[ , c("age", "Year") := lapply(.SD, as.numeric),
              .SDcols = c("age", "Year")]

### calculate median and 95-% confidence intervals
dt.incidence.CI <- dt.incidence[ , .(inc = median(incidence),
                                     inc.lcl = quantile(incidence, probs = 0.025),
                                     inc.ucl = quantile(incidence, probs = 0.975)),
                                 by = .(age, Year, Sex)]

### change type of age and year to numeric
dt.incidence[ , c("age", "Year") := lapply(.SD, as.numeric),
              .SDcols = c("age", "Year")]
dt.incidence.CI[ , c("age", "Year") := lapply(.SD, as.numeric),
                 .SDcols = c("age", "Year")]

### plot results
plt.inc <- ggplot(data = dt.incidence.CI[age < 81],
                  aes(x = age, y = inc * 1000,
                      color = factor(Sex, levels = c("Women", "Men")),
                      fill = factor(Sex, levels = c("Women", "Men")))) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = inc.lcl * 1000, ymax = inc.ucl * 1000), alpha = 0.3,
              show.legend = FALSE) +
  facet_wrap(vars(Year)) +
  xlab("Age in years") +
  ylab("Incidence per 1,000 person years") +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 10))

ggsave("plt.inc_withCI.tiff", plt.inc, compression = "lzw", height = 10, width = 20, units = "cm")
ggsave("plt_inc_withCI.svg", plt.inc, height = 10, width = 20, units = "cm")
ggsave("plt_inc_withCI.pdf", plt.inc, height = 10, width = 20, units = "cm")
ggsave("plt_inc_withCI.png", plt.inc, height = 10, width = 20, units = "cm")

#the best graphic representation of data is svg--the resolution is better
