# load packages
library(readxl)
library(data.table)
library(ggplot2)

# load prevalence data
prev <- data.table(read_xlsx("data/Final_Results_2023.xlsx"))
prev
setorder(prev, SEX, AGE_GROUP)

#plot(prev[1:17, DIABETES])
#plot(prev[18:34, DIABETES])
prev[ , age := 2.5 + as.numeric(substr(AGE_GROUP, 1, 2))]

## 1.2 Functions for prevalence data ----
logit  <- function(x){log(x/(1-x)) }
expit  <- function(x){1/(1+exp(-x))}

### function returns time, age and sex-specific prevalence
pmod <- lm(DIABETES ~ ns(age, df = 4) * Year * Sex, data = prev)
summary(pmod)
prev$link_lm <- predict(pmod)

get_prev <- function(t, a, p.Sex, rnd){
  a_ <- prev.predicted[Year %in% t & Sex == p.Sex, Age]
  p_ <- prev.predicted[Year %in% t & Sex == p.Sex, link_lm]
  #p.SE_ <- prev.predicted[Year %in% t & Sex == p.Sex, SE]

  ln.p <- approx(a_, p_, xout = a, rule = 2)$y
  #ln.p.SE <- approx(a_, p.SE_, xout = a, rule = 2)$y

  return(expit(ln.p + rnd * ln.p.SE))
}

plt.input_prev <- ggplot(data = prev[Age > 19 & Age < 80],
                         aes(x = Age, y = expit(link_lm) * 100)) +
  geom_line() +
  geom_ribbon(aes(ymin = expit(link_lm-1.96*SE) * 100,
                  ymax = expit(link_lm+1.96*SE) * 100),
              alpha = 0.4) +
  ylab("Prevalence") +
  xlab("Age in years") +
  theme_bw() +
  theme(strip.text = element_text(size = 13), axis.text = element_text(size = 13),
        axis.title = element_text(size = 13),
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_y_continuous(breaks= seq(0, 35, 5)) +
  facet_grid(rows = vars(Sex), cols = vars(Year))

### function returns change prevalence between survey years
get_dp <- function(t1, t2, dp.a, dp.Sex, rnd.p.1, rnd.p.2){
  dT  <- t2 - t1
  p1  <- get_prev(t1, dp.a - dT / 2, p.Sex = dp.Sex, rnd = rnd.p.1)
  p2  <- get_prev(t2, dp.a + dT / 2, p.Sex = dp.Sex, rnd = rnd.p.2)
  dp  <- (p2 - p1) / dT
  return(dp)
}

### function returns  prevalence midway between survey years
get_p.mid <- function(t1, t2, p.a, dp.Sex, rnd.p.1, rnd.p.2){
  dT  <- t2 - t1
  p1  <- get_prev(t1, p.a - dT / 2, p.Sex = dp.Sex, rnd = rnd.p.1)
  p2  <- get_prev(t2, p.a + dT / 2, p.Sex = dp.Sex, rnd = rnd.p.2)
  pt   <- 0.5 * (p2 + p1)
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
r_vec <- rep(NA, 1000000)
for(i in 1:1000000){
  r_i_ <- get_R(80)
  r_vec[i] <- r_i_
}
range(r_vec)
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
gen_mort[ , Age := as.numeric(substr(Age, 1, 3))]
setorder(gen_mort, Time, Sex, Age)
gen_mort <- gen_mort[, .(Time, Age, Sex, Value)]

fct.gen_mort <- function(t, a, sex_mort) {
  mort_dat <- gen_mort[Time == t & Sex == sex_mort]
  mort_    <- approx(mort_dat$Age, mort_dat$Value, xout = a, rule = 2)$y
  return(mort_)
}

fct.gen_mort(2010, 96:100, sex_mort = "Female")


plt.mortality <- ggplot(data = gen_mort, aes(x = Age, y = log(Value), colour = Sex)) +
  geom_line(linewidth = 1) +
  #geom_point() +
  facet_wrap(vars(Time), ncol = 2)

plt.heat <- ggplot(data = gen_mort, aes(x = Time, y = Age, fill = log(Value))) +
  geom_tile() +
  scale_fill_viridis_c()
