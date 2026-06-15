
# OBS:
# From the folder "Data", the files "Results_ENDES2016-Results_ENDES2023",
# have one of two workbooks that contain generic information onn the total
# number of people, the breakdown for men/women, and diabetes/no diabetes
# I have used that information to manually create "table1_clean.xlsx" document"

library(readxl)

table1_clean <- read_excel("Data/table1_clean.xlsx")

saveRDS(table1_clean, file = "Data/table1_clean.rds")
