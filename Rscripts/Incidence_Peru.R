# load packages
library(readxl)
library(data.table)

# load prevalence data
prev <- data.table(read_xlsx("data/Final_Results_2023.xlsx"))
prev
setorder(prev, SEX, AGE_GROUP)

plot(prev[1:17, DIABETES])
plot(prev[18:34, DIABETES])

