install.packages("tidyverse")
install.packages(c("haven", "survey"))

library(tidyverse)
library(haven)
library(survey)

##################################

setwd("C:/Users/Usuario/Documents/OBEMET Centro de Obesidad y Salud Metabólica/Artículos/Incidencia Diabetes Peru Colombia/ENDES/2023")

getwd()


# Load the ENDES 2023 datasets using the read_sav() function

SALUD <- read_sav("CSALUD01_2023.sav")
MUJER_OBS <- read_sav("RE223132_2023.SAV")
MUJER_LAC <- read_sav("REC42_2023.SAV")
MUJER_ANT <- read_sav("RECH5_2023.SAV")
PERSONA <- read_sav("RECH1_2023.SAV")
VIVIENDA <- read_sav("RECH23_2023.SAV")
HOGAR <- read_sav("RECH0_2023.SAV")

# Process the linking variables (EXTRACTING HHID AND QSNUMERO FROM CASEID):
# 1. Remove the last 3 characters of CASEID to create HHID
# 2. Extract the last 2 characters of CASEID to create QSNUMERO
# 3. Convert the text variable to numeric format
MUJER_OBS$HHID <- str_sub(MUJER_OBS$CASEID,1,(str_length(MUJER_OBS$CASEID)-3))
MUJER_OBS$QSNUMERO <- str_sub(MUJER_OBS$CASEID,-2,-1)
MUJER_OBS$QSNUMERO <- as.numeric(MUJER_OBS$QSNUMERO)
MUJER_LAC$HHID <- str_sub(MUJER_LAC$CASEID,1,(str_length(MUJER_LAC$CASEID)-3))
MUJER_LAC$QSNUMERO <- str_sub(MUJER_LAC$CASEID,-2,-1)
MUJER_LAC$QSNUMERO <- as.numeric(MUJER_LAC$QSNUMERO)
PERSONA$QSNUMERO <- PERSONA$HVIDX
MUJER_ANT$QSNUMERO <- MUJER_ANT$HA0

# Merge datasets using the linking identifiers HHID and QSNUMERO
BASE1 <- left_join(PERSONA, SALUD, by = c("HHID","QSNUMERO"))
BASE2 <- left_join(MUJER_OBS, MUJER_LAC, by = c("HHID", "QSNUMERO"))
BASE2 <- left_join(BASE2, MUJER_ANT, by = c("HHID", "QSNUMERO"))
BASE3 <- left_join(HOGAR, VIVIENDA, by = "HHID")
BASE4 <- left_join(BASE1,BASE2,by = c("HHID","QSNUMERO"))
BASE5 <- left_join(BASE3,BASE4,by = 'HHID')


# Filter dataset by removing rows with NA values in QS109, QS23, or QSSEXO
BASE_FINAL <- BASE5 %>%
  filter(!is.na(QS109), !is.na(QS23), !is.na(QSSEXO), !is.na(Peso15_amas))

nrow(BASE_FINAL)
unique(BASE_FINAL$QS109)
table(BASE_FINAL$QS109, useNA = "ifany")


## Create the Diabetes variable
BASE_FINAL <- BASE_FINAL %>%
  mutate(DIABETES = ifelse(QS109 == 1, 1, 0))  # 1 = Has diabetes, 0 = Does not have diabetes

table(BASE_FINAL$DIABETES, useNA = "ifany")


###Create a simplified dataset with only the necessary variables
ENDES <- BASE_FINAL %>%
  select(HHID, QSNUMERO, CONGLOMERADO = HV001, ESTRATO = HV022,
         PONDERACION = Peso15_amas, SEXO = HV104, EDAD = QS23, DIABETES,AREA_RESIDENCIA = HV025,
         QUINTIL_BIENESTAR = HV270, REGION_NATURAL = SHREGION)


# Specify the survey design
diseño <- svydesign(id =~ CONGLOMERADO, strata =~ ESTRATO, weights=~ PONDERACION, data=ENDES)

##Create age group variable
ENDES <- ENDES %>%
  mutate(
    GRUPO_EDAD = cut(
      EDAD,
      breaks = c(15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
                 65, 70, 75, 80, 85, 90, 95, 100),
      right = FALSE,
      labels = c("15-19","20-24","25-29","30-34","35-39","40-44",
                 "45-49","50-54","55-59","60-64","65-69","70-74",
                 "75-79","80-84","85-89","90-94","95-99")
    )
  )

## Verifying
table(ENDES$GRUPO_EDAD, useNA = "ifany")

# Specifying the survey design
diseño <- svydesign(id =~ CONGLOMERADO, strata =~ ESTRATO, weights=~ PONDERACION, data=ENDES)

# Calculate prevalence and 95% CI by sex and age group (%)
resultados <- svyby(
  ~ DIABETES,                  # Outcome variable (0/1)
  ~ SEXO + GRUPO_EDAD,         # Stratification variables
  design = diseño,             # Survey design
  FUN = svymean,
  na.rm = TRUE,
  vartype = c("ci")            # Confidence intervals
)

# Convert to %
resultados <- resultados %>%
  mutate(
    DIABETES = DIABETES * 100,   # Convert proportion to percentage
    ci_l = ci_l * 100,          # Lower bound of CI in percentage
    ci_u = ci_u * 100           # Upper bound of CI in percentage
  )

# Check results in percentage format
print(resultados)

# Load library
library(writexl)

# Save the table in Excel
write_xlsx(resultados, "resultados_diabetes.xlsx")



# Calculate prevalence and 95% CI using logit-transformed intervals
resultados <- svyby(
  ~ DIABETES,                  # Outcome variable (0/1)
  ~ SEXO + GRUPO_EDAD,         # Stratification variables
  design = diseño,             # Survey design
  FUN = svymean,
  na.rm = TRUE,
  vartype = c("ci"),           # Intervalos de confianza
  method = "logit"             # Use logit method for CI
)

# Convert to %
resultados <- resultados %>%
  mutate(
    DIABETES = DIABETES * 100,  # Convert to %
    ci_l = pmax(ci_l * 100, 0), # Ensure lower bounds are not negative
    ci_u = pmin(ci_u * 100, 100) # Ensure upper bounds do not exceed 100%
  )

# Verifyin results
print(resultados)


# Load library
library(writexl)

# Save the table in Excel
write_xlsx(resultados, "resultados22_diabetes.xlsx")


######################
#######################
#########################
###########################


###CROSS-CHECK:
## OVERALL DIABETES PREVALENCE BY SEX:
# Prevalence by sex
resultados_sexo <- svyby(
  ~ DIABETES,          # Outcome variable
  ~ SEXO,              # Stratifying by sex
  design = diseño,     # Survey design
  FUN = svymean,
  na.rm = TRUE,
  vartype = c("ci"),   # Confidence intervals
  method = "logit"     # Logit method for 95%CI
)

# Convert to %
resultados_sexo <- resultados_sexo %>%
  mutate(
    DIABETES = DIABETES * 100,
    ci_l = pmax(ci_l * 100, 0),
    ci_u = pmin(ci_u * 100, 100)
  )

# Check results
print(resultados_sexo)



###CROSS-CHECK:
## DIABETES PREVALENCE BY AGE GROUP:
## 15-19, 20-29, 30-39, 40-49, 50-59, 60 Y MÁS

# Define new age group categories
ENDES <- ENDES %>%
  mutate(
    GRUPO_EDAD_ESPECIAL = case_when(
      EDAD >= 15 & EDAD <= 19 ~ "15-19",
      EDAD >= 20 & EDAD <= 29 ~ "20-29",
      EDAD >= 30 & EDAD <= 39 ~ "30-39",
      EDAD >= 40 & EDAD <= 49 ~ "40-49",
      EDAD >= 50 & EDAD <= 59 ~ "50-59",
      EDAD >= 60 ~ "60 y más",
      TRUE ~ NA_character_
    ),
    GRUPO_20MAS = ifelse(EDAD >= 20, "20 y más", "Menos de 20")
  )

# Specify the survey design
diseño <- svydesign(id =~ CONGLOMERADO, strata =~ ESTRATO, weights=~ PONDERACION, data=ENDES)

# Prevalence by GRUPO_EDAD_ESPECIAL
resultados_edad <- svyby(
  ~ DIABETES,
  ~ GRUPO_EDAD_ESPECIAL,
  design = diseño,
  FUN = svymean,
  na.rm = TRUE,
  vartype = c("ci"),
  method = "logit"
)

# Convert to %
resultados_edad <- resultados_edad %>%
  mutate(
    DIABETES = DIABETES * 100,
    ci_l = pmax(ci_l * 100, 0),
    ci_u = pmin(ci_u * 100, 100)
  )

# Check results
print(resultados_edad)

