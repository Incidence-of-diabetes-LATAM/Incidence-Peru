rm(list = ls())

# load packages
library(readxl)
library(data.table)
library(ggplot2)
library(splines)

## 1.2 Functions for prevalence data ----
logit  <- function(x){log(x/(1-x)) }
expit  <- function(x){1/(1+exp(-x))}

# List of years
years <- 2015:2023

# Initialize an empty data.table
dt.endes_all <- data.table()

# Loop through each year and bind to endes_all
for (year in years) {
  # Construct the file path
  file_path <- paste0("data/Results_ENDES", year, ".xlsx")

  # Read and process the data
  dt.prev <- data.table(read_xlsx(file_path))
  dt.prev[, `:=`(
    age = 2.5 + as.numeric(substr(AGE_GROUP, 1, 2)),
    year = year+0.5,
    prev = 0.01 * DIABETES,
    SEX = fcase(
      SEX == 1, "Men",
      SEX == 2, "Women"
    )
    )]

  dt.prev[, prev_logit := logit(prev)]

  # Append to the full dataset
  dt.endes_all <- rbind(dt.endes_all, dt.prev, fill = TRUE)
}

dt.endes_all[ , SE := (prev - ci_l/100) / 1.96]
dt.endes_all[ , SE := (prev - ci_l/100) / 1.96]


### function returns time, age and sex-specific prevalence
pmod <- lm(prev_logit ~ poly(age, 3) + SEX + year, data = dt.endes_all)

summary(pmod)
AIC(pmod)
BIC(pmod)
dt.endes_all$prev_logit.predicted <- predict(pmod)

plt.input_prev <- ggplot(data = dt.endes_all, aes(x = age, y = expit(prev_logit.predicted) * 100)) +
  geom_line() +
  geom_point(aes(y = DIABETES)) +
  #geom_ribbon(aes(ymin = ci_l, ymax = ci_u), alpha = 0.4) +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u)) +
  ylab("Prevalence") +
  xlab("Age in years") +
  theme_bw() +
  theme(strip.text = element_text(size = 10), axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_y_continuous(breaks= seq(0, 35, 5)) +
  facet_grid(rows = vars(SEX), cols = vars(year))

ggsave("Results/plt.input_prev.tiff", plt.input_prev, compression = "lzw", height = 10, width = 20, units = "cm")

get_prev <- function(t, a, p.Sex){
  ln.p  <- expit(predict(pmod, newdata = data.frame(year = t, age = a, SEX = p.Sex)))
  return(ln.p)
}

plot(20:80, get_prev(2023.5, 20:80, p.Sex = "Men"), ylim = c(0, 0.25))

### function returns change prevalence between survey years
get_dp <- function(t1, t2, dp.a, dp.Sex){
  dT  <- t2 - t1
  p1  <- get_prev(t1, dp.a - dT / 2, p.Sex = dp.Sex)
  p2  <- get_prev(t2, dp.a + dT / 2, p.Sex = dp.Sex)
  dp  <- (p2 - p1) / dT
  return(dp)
}

### function returns  prevalence midway between survey years
get_p.mid <- function(t1, t2, p.a, dp.Sex){
  dT  <- t2 - t1
  p1  <- get_prev(t1, p.a - dT / 2, p.Sex = dp.Sex)
  p2  <- get_prev(t2, p.a + dT / 2, p.Sex = dp.Sex)
  pt  <- 0.5 * (p2 + p1)
  return(pt)
}

## 1.3 Data and function for mortality rate ratio (diab vs. no diab) ----
get_R <- function(a_){
  R30 <- runif(1, min = 7.5, max = 8.5)
  R90 <- runif(1, min = 2, max = 2.5)
  logR_ <- approx(c(30, 90), log(c(R30, R90)), xout = a_, rule = 2)$y
  return(exp(logR_))
}

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
out.ages <- 18:84 # age range in years
out.years <- seq(2016, 2023, 1) # year of incidence estimation
svy.years <- seq(2015.5, 2023.5, 1) # ENDES survey years
nboots <- 1000 # no. of repetitions

i.res <- array(dim = c(length(out.ages), length(svy.years) - 1, 2, nboots),
               dimnames = list(c(out.ages), c(out.years), c("Men", "Women"), 1:nboots))

### apply partial differential equation (PDE) to input data
set.seed(68763)
for(bb in 1:nboots){ # nboots loops for random error estimation
  for(ii in c("Men", "Women")){ # repeat analysis for men and women separately
    for(tt in 1:(length(svy.years) - 1)){ # and midpoints between surveys
      pdat_ <- dt.endes_all
      rnd.vec <- rnorm(nrow(dt.endes_all))
      pdat_[ , prev_err := pmax(1e-10, prev + rnd.vec * SE)]
      pmod <- lm(logit(prev_err) ~ poly(age, 3) + year + SEX, data = pdat_)


      MRR_  <- get_R(out.ages) # get MRR

      t1_   <- svy.years[tt] # define t1
      t2_   <- svy.years[tt + 1] # define t2
      t_    <- 0.5 * (t1_ + t2_) # midpoint between t1 and t2

      # get general mortality at midpoint between t1 and t2
      m_    <- fct.gen_mort(t = round(t_), a = out.ages, sex_mort = ii)

      # get prevalence at midpoint between t1 and t2
      p_    <- get_p.mid(t1 = t1_, t2 = t2_, p.a = out.ages, dp.Sex = ii)

      # get change in prevalence at midpoint between t1 and t2
      dp_   <- get_dp(t1 = t1_, t2 = t2_, dp.a = out.ages, dp.Sex = ii)

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

### plot age-specific incidence
plt.inc_age <- ggplot(data = dt.incidence.CI,
                  aes(x = age, y = inc * 1000,
                      color = factor(Sex, levels = c("Women", "Men")),
                      fill = factor(Sex, levels = c("Women", "Men")))) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = inc.lcl * 1000, ymax = inc.ucl * 1000), alpha = 0.2,
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

ggsave("Results/plt.inc_age.tiff", plt.inc_age, compression = "lzw", height = 12, width = 20, units = "cm")


## 2.1 Estimate age-standardized incidence rate ----

# population data for age standardization from here: https://population.un.org/dataportal/data/indicators/47/locations/604/start/2010/end/2023/table/pivotbyvariant?df=99779a07-5b67-43bb-93aa-b0b03b82a778
dt.population <- fread("Data/unpopulation_dataportal_population_Peru.csv")
dt.population <- dt.population[Age %in% out.ages & Time == 2023, .(Sex, Age, Value)]
setnames(dt.population, "Value", "population")

# calculate weights: proportion of age group in relation to total population
dt.population[ , totalpop := sum(population), by = Sex]
dt.population[ , weight := population / totalpop]
dt.population[ , sum(weight), by = Sex]
dt.population[ , Sex := fcase(
  Sex == "Male", "Men",
  Sex == "Female", "Women"
)]
dt.population[ , Age := as.numeric(Age)]

# merge weights with incidence
dt.incidence.weights <- merge(dt.incidence,
                              dt.population[ , .(Sex, Age, weight)],
                              by.x = c("Sex", "age"),
                              by.y = c("Sex", "Age"))

# check if sum(weight) = 1
dt.incidence.weights[boot == 500 , sum(weight), by = .(Sex, Year)]

# calculate age-standardized incidence
dt.incidence.standardized <-
  dt.incidence.weights[ , .(incidence = sum(incidence * weight)),
                        by = .(Sex, Year, boot)]


dt.incidence.standardized <-
  merge(dt.incidence.standardized, dt.incidence.standardized[Year == 2016],
        by = c("Sex", "boot"))

names(dt.incidence.standardized) <-  c("Sex", "boot", "Year", "incidence",
                                       "Year.ref", "incidence.ref")

dt.incidence.standardized[ , IRR := incidence / incidence.ref]

dt.incidence.standardized.CI <- dt.incidence.standardized[ ,
                                                           .(inc = median(incidence),
                                                             inc.lcl = quantile(incidence, probs = 0.025),
                                                             inc.ucl = quantile(incidence, probs = 0.975),

                                                             IRR = median(IRR),
                                                             IRR.lcl = quantile(IRR, probs = 0.025),
                                                             IRR.ucl = quantile(IRR, probs = 0.975)),
                                                           by = .(Sex, Year)]
# table includes age-standardized incidence rate (inc) with 95%-CI (inc.lcl, inc.ucl) and
# incidence rate ratio (IRR) with 95%-CI (IRR.lcl, IRR.ucl). The IRRs indicate the temporal trend in
# in incidence with 2016 as the reference year.
write.table(dt.incidence.standardized.CI, file = "Results/inc.standardized.csv")

# plot for temporal trend in age-standardized incidence
plt.standardized <- ggplot(data = dt.incidence.standardized.CI,
                      aes(x = Year, y = inc * 1000,
                          color = factor(Sex, levels = c("Women", "Men")),
                          fill = factor(Sex, levels = c("Women", "Men")))) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(ymin = inc.lcl * 1000, ymax = inc.ucl * 1000), alpha = 0.2,
              show.legend = FALSE) +
  xlab("Age in years") +
  ylab("Incidence per 1,000 person years") +
  ylim(c(0,17)) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.title = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        legend.text = element_text(size = 11))

ggsave("Results/plt.standardized.tiff", plt.standardized, compression = "lzw", height = 10, width = 20, units = "cm")


## 2.2 Incidence in 10 year age steps and incidence rate ratios for temporal trend ----
dt.inc.10yearage <- dt.incidence.weights[age %in% seq(20, 80, 10) & Year %in% c(2016, 2023)]
dt.inc.10yearage <- merge(dt.inc.10yearage,
                          dt.inc.10yearage[Year == 2016,
                                           .(Sex, age, Year, boot, incidence)],
                          by = c("Sex", "age", "boot"))

names(dt.inc.10yearage) <-  c("Sex", "age", "boot", "Year","incidence",
                              "weight", "Year.ref", "incidence.ref")
dt.inc.10yearage[ , IRR := incidence / incidence.ref]


# remove NA
nrow(dt.inc.10yearage)
dt.inc.10yearage <- dt.inc.10yearage[!(is.na(IRR))]
nrow(dt.inc.10yearage)

# remove Inf
dt.inc.10yearage[(IRR == "Inf")]
dt.inc.10yearage <- dt.inc.10yearage[!(IRR == "Inf")]
nrow(dt.inc.10yearage)

# calculate 95%-CI for incidence in age groups
dt.inc.10yearage.CI <- dt.inc.10yearage[ ,
                                         .(inc = median(incidence),
                                           inc.lcl = quantile(incidence, probs = 0.025),
                                           inc.ucl = quantile(incidence, probs = 0.975),

                                           IRR = median(IRR),
                                           IRR.lcl = quantile(IRR, probs = 0.025),
                                           IRR.ucl = quantile(IRR, probs = 0.975)
                                         ),
                                         by = .(Sex, age, Year)]

# table includes incidence rate (inc) with 95%-CI (inc.lcl, inc.ucl) and
# incidence rate ratio (IRR) with 95%-CI (IRR.lcl, IRR.ucl). The IRR indicates the temporal trend in
# in incidence comparin inc in 2023 vs 2016.
write.table(dt.inc.10yearage.CI, file = "Results/inc.agespecific.csv")

