install.packages("tidyverse")
install.packages(c("haven", "survey"))

library(tidyverse)
library(haven)
library(survey)

##################################

setwd("C:/Users/Usuario/Documents/OBEMET Centro de Obesidad y Salud Metabólica/Artículos/Incidencia Diabetes Peru Colombia/ENDES/2023")

getwd()


# Cargamos las bases de datos de la ENDES 2023 con la función read_sav()

SALUD <- read_sav("CSALUD01_2023.sav")
MUJER_OBS <- read_sav("RE223132_2023.SAV")
MUJER_LAC <- read_sav("REC42_2023.SAV")
MUJER_ANT <- read_sav("RECH5_2023.SAV")
PERSONA <- read_sav("RECH1_2023.SAV")
VIVIENDA <- read_sav("RECH23_2023.SAV")
HOGAR <- read_sav("RECH0_2023.SAV")

# Transformamos las variables de unión (DESCOMPOSICIÓN DEL CASEID EN HHID Y QSNUMERO):
# 1. Descartamos últimos 3 caracteres de CASEID y lo convertimos a HHID
# 2. Extraemos los últimos 2 caracteres de CASEID y lo convertimos a QSNUMERO
# 3. Convertimos la variable de texto a numérica
MUJER_OBS$HHID <- str_sub(MUJER_OBS$CASEID,1,(str_length(MUJER_OBS$CASEID)-3))
MUJER_OBS$QSNUMERO <- str_sub(MUJER_OBS$CASEID,-2,-1)
MUJER_OBS$QSNUMERO <- as.numeric(MUJER_OBS$QSNUMERO)
MUJER_LAC$HHID <- str_sub(MUJER_LAC$CASEID,1,(str_length(MUJER_LAC$CASEID)-3))
MUJER_LAC$QSNUMERO <- str_sub(MUJER_LAC$CASEID,-2,-1)
MUJER_LAC$QSNUMERO <- as.numeric(MUJER_LAC$QSNUMERO)
PERSONA$QSNUMERO <- PERSONA$HVIDX
MUJER_ANT$QSNUMERO <- MUJER_ANT$HA0

# Realizamos la unión mediante los identificadores o llaves HHID y QSNUMERO
BASE1 <- left_join(PERSONA, SALUD, by = c("HHID","QSNUMERO"))
BASE2 <- left_join(MUJER_OBS, MUJER_LAC, by = c("HHID", "QSNUMERO"))
BASE2 <- left_join(BASE2, MUJER_ANT, by = c("HHID", "QSNUMERO"))
BASE3 <- left_join(HOGAR, VIVIENDA, by = "HHID")
BASE4 <- left_join(BASE1,BASE2,by = c("HHID","QSNUMERO"))
BASE5 <- left_join(BASE3,BASE4,by = 'HHID')


# Filtrar eliminando filas con NA en QS109, QS23 o QSSEXO
BASE_FINAL <- BASE5 %>%
  filter(!is.na(QS109), !is.na(QS23), !is.na(QSSEXO), !is.na(Peso15_amas))

nrow(BASE_FINAL)
unique(BASE_FINAL$QS109)
table(BASE_FINAL$QS109, useNA = "ifany")


##Calcular variable Diabetes
BASE_FINAL <- BASE_FINAL %>%
  mutate(DIABETES = ifelse(QS109 == 1, 1, 0))  # 1 = Sí tiene diabetes, 0 = No

table(BASE_FINAL$DIABETES, useNA = "ifany")


###Crear una base de datos simplificada con lo necesario
ENDES <- BASE_FINAL %>%
  select(HHID, QSNUMERO, CONGLOMERADO = HV001, ESTRATO = HV022,
         PONDERACION = Peso15_amas, SEXO = HV104, EDAD = QS23, DIABETES,AREA_RESIDENCIA = HV025,
         QUINTIL_BIENESTAR = HV270, REGION_NATURAL = SHREGION)


# Especificamos el diseño muestral de la encuesta
diseño <- svydesign(id =~ CONGLOMERADO, strata =~ ESTRATO, weights=~ PONDERACION, data=ENDES) 

##Crear la variable de grupos de edad
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

## Verificando
table(ENDES$GRUPO_EDAD, useNA = "ifany")

# Especificamos el diseño muestral de la encuesta
diseño <- svydesign(id =~ CONGLOMERADO, strata =~ ESTRATO, weights=~ PONDERACION, data=ENDES) 

# Calcular prevalencias e IC95% por sexo y grupo de edad en porcentaje
resultados <- svyby(
  ~ DIABETES,                  # La variable de interés (0/1)
  ~ SEXO + GRUPO_EDAD,         # Variables de estratificación
  design = diseño,             # Diseño actualizado
  FUN = svymean, 
  na.rm = TRUE, 
  vartype = c("ci")            # Intervalos de confianza
)

# Convertir a porcentaje
resultados <- resultados %>%
  mutate(
    DIABETES = DIABETES * 100,   # Convertir la proporción a porcentaje
    ci_l = ci_l * 100,          # Límite inferior del IC en porcentaje
    ci_u = ci_u * 100           # Límite superior del IC en porcentaje
  )

# Verificar resultados en porcentaje
print(resultados)

# Cargar librería
library(writexl)

# Guardar la tabla en Excel
write_xlsx(resultados, "resultados_diabetes.xlsx")



# Calcular prevalencias e IC95% usando logit-transformed intervals
resultados <- svyby(
  ~ DIABETES,                  # La variable de interés (0/1)
  ~ SEXO + GRUPO_EDAD,         # Variables de estratificación
  design = diseño,             # Diseño actualizado
  FUN = svymean, 
  na.rm = TRUE, 
  vartype = c("ci"),           # Intervalos de confianza
  method = "logit"             # Usar método logit para los IC
)

# Convertir las proporciones e IC a porcentajes
resultados <- resultados %>%
  mutate(
    DIABETES = DIABETES * 100,  # Convertir proporción a porcentaje
    ci_l = pmax(ci_l * 100, 0), # Asegurar que los límites inferiores no sean negativos
    ci_u = pmin(ci_u * 100, 100) # Asegurar que los límites superiores no excedan 100%
  )

# Verificar resultados
print(resultados)


# Cargar librería
library(writexl)

# Guardar la tabla en Excel
write_xlsx(resultados, "resultados22_diabetes.xlsx")


######################
#######################
#########################
###########################


###CROSS-CHECK:
##PREVALENCIA DE DIABETES GENERAL PARA HOMBRES Y MUJERES: 
# Prevalencias por sexo
resultados_sexo <- svyby(
  ~ DIABETES,          # Variable de interés
  ~ SEXO,              # Estratificación por sexo
  design = diseño,     # Diseño muestral
  FUN = svymean, 
  na.rm = TRUE, 
  vartype = c("ci"),   # Intervalos de confianza
  method = "logit"     # Método logit para IC
)

# Convertir a porcentaje
resultados_sexo <- resultados_sexo %>%
  mutate(
    DIABETES = DIABETES * 100,
    ci_l = pmax(ci_l * 100, 0),
    ci_u = pmin(ci_u * 100, 100)
  )

# Ver resultados
print(resultados_sexo)




###CROSS-CHECK:
##PREVALENCIA DE DIABETES PARA GRUPOS DE EDAD: 
## 15-19, 20-29, 30-39, 40-49, 50-59, 60 Y MÁS

# Agregar nuevas categorías de edad
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

# Especificamos el diseño muestral de la encuesta
diseño <- svydesign(id =~ CONGLOMERADO, strata =~ ESTRATO, weights=~ PONDERACION, data=ENDES) 

# Prevalencias por GRUPO_EDAD_ESPECIAL
resultados_edad <- svyby(
  ~ DIABETES,
  ~ GRUPO_EDAD_ESPECIAL,
  design = diseño,
  FUN = svymean,
  na.rm = TRUE,
  vartype = c("ci"),
  method = "logit"
)

# Convertir a porcentaje
resultados_edad <- resultados_edad %>%
  mutate(
    DIABETES = DIABETES * 100,
    ci_l = pmax(ci_l * 100, 0),
    ci_u = pmin(ci_u * 100, 100)
  )

# Ver resultados
print(resultados_edad)

