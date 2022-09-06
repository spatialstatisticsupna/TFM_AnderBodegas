rm(list = ls())

library(INLA)
library(sf)
library(tmap)
library(tmaptools)

#------------------------------------------------------------------------------#
# LECTURA
#------------------------------------------------------------------------------#

# Mis rutas
ruta_datos <- "Datos/"                  # "C:/Users/Ander/Documents/TFM/Datos/"
ruta_resultados <- "Resultados/"        # "C:/Users/Ander/Documents/TFM/Resultados/"

# Leemos los datos de cáncer de mujeres y hombres
final_gb_f <- read.csv(paste(ruta_datos, "Final/final_gb_f.csv", sep = ""), header = T)
final_gb_m <- read.csv(paste(ruta_datos, "Final/final_gb_m.csv", sep = ""), header = T)

# Leemos la cartografía de Gran Bretaña
carto_gb <- st_read(paste(ruta_datos, "Cartografía/carto_gb.shp", sep = ""))

# Leemos la matriz de adyacencia (tiene añadidas las conexiones de las islas etc.)
adj_gb <- as.matrix(read.table(paste(ruta_datos, "Cartografía/adj_gb.txt", sep = "")))

# 

#------------------------------------------------------------------------------#
# DATOS
#------------------------------------------------------------------------------#

# Ordenamos por Región y Año
final_gb_f <- final_gb_f[order(final_gb_f$Code), ]
final_gb_f <- final_gb_f[order(final_gb_f$Year), ]
final_gb_m <- final_gb_m[order(final_gb_m$Code), ]
final_gb_m <- final_gb_m[order(final_gb_m$Year), ]
carto_gb <- carto_gb[order(carto_gb$Code), ]

# Nos quedamos con los datos para el ejemplo espacial (Nuevos casos de cáncer de pulmón femenino 2019)
inc_lung_2019_f <- final_gb_f[which(final_gb_f$Year == 2019), c(1, 3, 22)]

# Nos quedamos con los datos para el ejemplo espacio-temporal (Nuevos casos de cáncer de pulmón femenino 2002-2019)
inc_lung_f <- final_gb_f[, c(1:3, 22)]

# Nos quedamos con el resto de datos para el análisis final (Nuevos casos de cáncer de pulmón masculino y femenino 2002-2019)
mor_lung_f <- final_gb_f[, c(1:3, 23)]
inc_lung_m <- final_gb_m[, c(1:3, 22)]
mor_lung_m <- final_gb_m[, c(1:3, 23)]

# Cambiamos el nombre de la variable para que cuadre con las fórmulas
names(inc_lung_2019_f)[3] <- "Obs"
names(inc_lung_f)[4] <- "Obs"
names(mor_lung_f)[4] <- "Obs"
names(inc_lung_m)[4] <- "Obs"
names(mor_lung_m)[4] <- "Obs"

# Añadimos las variales ID.Code e ID.Year para poder usar las interacciones
aux_code <- data.frame(unique(inc_lung_f$Code), 1:142)
aux_year <- data.frame(unique(inc_lung_f$Year), 1:18)
names(aux_code) <- c("Code", "ID.Code")
names(aux_year) <- c("Year", "ID.Year")
inc_lung_f <- merge(inc_lung_f, aux_code)
inc_lung_f <- merge(inc_lung_f, aux_year)
inc_lung_m <- merge(inc_lung_m, aux_code)
inc_lung_m <- merge(inc_lung_m, aux_year)
mor_lung_f <- merge(mor_lung_f, aux_code)
mor_lung_f <- merge(mor_lung_f, aux_year)
mor_lung_m <- merge(mor_lung_m, aux_code)
mor_lung_m <- merge(mor_lung_m, aux_year)

#------------------------------------------------------------------------------#
# EJEMPLO MEDIDAS CLÁSICAS (2019)
#------------------------------------------------------------------------------#

# Añadimos a la cartografía los casos osbservados, la población y la tasa cruda
# (usamos merge para que no haya errores asignando los valores a las regiones)
aux <- data.frame(inc_lung_2019_f$Code,
                  inc_lung_2019_f$Obs,
                  inc_lung_2019_f$Population,
                  100000 * inc_lung_2019_f$Obs / inc_lung_2019_f$Population)
names(aux) <- c("Code", "Casos", "Población", "TC")
carto <- merge(carto_gb, aux)

# Creamos los tres mapas y los guardamos en la carpeta de resultados
paleta <- get_brewer_pal("Blues", n = 7, contrast = 0.8)
archivo <- paste(ruta_resultados, "Mapas/Ejemplo_Clasicas.jpeg", sep = "")
jpeg(file = archivo, width = 2000, height = 1200)
print(tm_shape(carto) + 
        tm_polygons(c("Casos", "Población", "TC"), style = "equal", n = 7, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 2, legend.text.size = 1.5, legend.position = c("left", "top")))
dev.off()

# Habría que indicarle los breaks manualmente

#------------------------------------------------------------------------------#
# EJEMPLO MODELOS ESPACIALES (2019)
#------------------------------------------------------------------------------#

# Ajustamos el modelo iCAR
codes <- unique(inc_lung_2019_f$Code)
form1 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb)
model1 <- inla(formula = form1, family = "poisson", E = Population, data = inc_lung_2019_f)

# Representamos las distribuciones marginales posteriores del valor base, del parámetro espacial correspondiente
# a la primera región (Barnsley) y del hiperparámetro. Lo guardamos en la carpeta de resultados
archivo <- paste(ruta_resultados, "Gráficos/Posteriores_SP.pdf", sep = "")
pdf(file = archivo, width = 10, height = 5)
par(mfrow = c(1, 3))
plot(model1$marginals.fixed$`(Intercept)`, type = "l", lwd = 1.5, main = "Valor base", xlab = expression(eta), ylab = expression(tilde(p)(paste(eta, "|", y))))
abline(v = model1$summary.fixed$`0.5quant`, col = "blue")
plot(model1$marginals.random$Code$index.1, type = "l", lwd = 1.5, main = "Parámetro espacial", xlab = expression(u[1]), ylab = expression(tilde(p)(paste(u[1], "|", y))))
abline(v = model1$summary.random$Code$`0.5quant`[1], col = "blue")
plot(model1$marginals.hyperpar$`Precision for Code`, type = "l", lwd = 1.5, main = "Hiperparámetro", xlab = expression(tau), ylab = expression(tilde(p)(paste(tau, "|", y))))
abline(v = model1$summary.hyperpar$`0.5quant`, col = "blue")
par(mfrow = c(1, 1))
dev.off()

# Valor base
model1$summary.fixed$`0.5quant`
exp(model1$summary.fixed$`0.5quant`)
100000 * exp(model1$summary.fixed$`0.5quant`)

# Parámetro espacial de Barnsley
model1$summary.random$Code$`0.5quant`[1]
exp(model1$summary.random$Code$`0.5quant`[1])

# Barnsley (manualmente y con fitted values)
100000 * exp(model1$summary.fixed$`0.5quant`) * exp(model1$summary.random$Code$`0.5quant`[1])
100000 * model1$summary.fitted.values$`0.5quant`[1]

# Hiperparámetro
model1$summary.hyperpar$`0.5quant`

# Ajustamos todos los modelos con el WAIC y el DIC

# i.i.d.
form0 <- Obs ~ f(Code, model = "iid", values = codes, graph = adj_gb)
model0 <- inla(formula = form0, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_2019_f)

# iCAR
form1 <- Obs ~ f(Code, model = "besag", values = codes, graph = adj_gb)
model1 <- inla(formula = form1, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_2019_f)

# BYM
form2 <- Obs ~ f(Code, model = "bym", values = codes, graph = adj_gb)
model2 <- inla(formula = form2, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_2019_f)

# Leroux
Q <- diag(abs(colSums(adj_gb))) - adj_gb
C <- diag(rep(1, nrow(Q))) - Q
form3 <- Obs ~ f(Code, model = "generic1", values = codes, Cmatrix = C)
model3 <- inla(formula = form3, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_2019_f)

# WAIC
model0$waic$waic
model1$waic$waic
model2$waic$waic
model3$waic$waic # menor

# DIC
model0$dic$dic
model1$dic$dic
model2$dic$dic
model3$dic$dic # menor

# Añadimos a la cartografía los resultados del modelo de Leroux
model <- model3
aux <- data.frame(model$summary.random$Code$ID,
                  100000 * model$summary.fitted.values)
names(aux) <- c("Code", "TA")
carto <- merge(carto_gb, aux)

# Creamos el mapa y lo guardamos en la carpeta de resultados
paleta <- get_brewer_pal("-RdYlGn", n = 7, contrast = 0.8)
archivo <- paste(ruta_resultados, "Mapas/Ejemplo_SP.jpeg", sep = "")
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("TA", style = "equal", n = 7, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()

# Habría que indicarle los breaks manualmente
# Habría que añadir el mapa de exeedance probabilities

#------------------------------------------------------------------------------#
# EJEMPLO MODELOS ESPACIO-TEMPORALES (2002-2019)
#------------------------------------------------------------------------------#

# Ajustamos el modelo iCAR + PA1 + TipoIV
years <- unique(inc_lung_f$Year)
form114 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
model114 <- inla(formula = form114, family = "poisson", E = Population, data = inc_lung_f)

# Representamos las distribuciones marginales posteriores del parámetro espacial correspondiente a la primera región
# (Barnsley), del parámetro temporal correspondiente al primer año (2002) y del parámetro correspondiente a la interacción
# Barnsley-2002
archivo <- paste(ruta_resultados, "Gráficos/Posteriores_SP_TP_param.pdf", sep = "")
pdf(file = archivo, width = 10, height = 5)
par(mfrow = c(1, 3))
plot(model114$marginals.random$Code$index.1, type = "l", lwd = 1.5, main = "Parámetro espacial", xlab = expression(u[1]), ylab = expression(tilde(p)(paste(u[1], "|", y))))
abline(v = model114$summary.random$Code$`0.5quant`[1], col = "blue")
plot(model114$marginals.random$Year$index.1, type = "l", lwd = 1.5, main = "Parámetro temporal", xlab = expression(phi[1]), ylab = expression(tilde(p)(paste(phi[1], "|", y))))
abline(v = model114$summary.random$Year$`0.5quant`[1], col = "blue")
plot(model114$marginals.random$ID.Code$index.1, type = "l", lwd = 1.5, main = "Parámetro de interacción", xlab = expression(delta[1][1]), ylab = expression(tilde(p)(paste(delta[1][1], "|", y))))
abline(v = model114$summary.random$ID.Code$`0.5quant`[1], col = "blue")
par(mfrow = c(1, 1))
dev.off()

# Representamos las distribuciones marginales posteriores del hiperparámetro espacial, temporal y de interacción
archivo <- paste(ruta_resultados, "Gráficos/Posteriores_SP_TP_hiperparam.pdf", sep = "")
pdf(file = archivo, width = 10, height = 5)
par(mfrow = c(1, 3))
plot(model114$marginals.hyperpar$`Precision for Code`[1:65, ], type = "l", lwd = 1.5, main = "Hiperparámetro espacial", xlab = expression(tau[u]), ylab = expression(tilde(p)(paste(tau[u], "|", y))))
abline(v = model114$summary.hyperpar$`0.5quant`[1], col = "blue")
plot(model114$marginals.hyperpar$`Precision for Year`[1:70, ], type = "l", lwd = 1.5, main = "Hiperparámetro temporal", xlab = expression(tau[phi]), ylab = expression(tilde(p)(paste(tau[phi], "|", y))))
abline(v = model114$summary.hyperpar$`0.5quant`[2], col = "blue")
plot(model114$marginals.hyperpar$`Precision for ID.Code`, type = "l", lwd = 1.5, main = "Hiperparámetro de interacción", xlab = expression(tau[delta]), ylab = expression(tilde(p)(paste(tau[delta], "|", y))))
abline(v = model114$summary.hyperpar$`0.5quant`[3], col = "blue")
par(mfrow = c(1, 1))
dev.off()

# Valor base
model114$summary.fixed$`0.5quant`
exp(model114$summary.fixed$`0.5quant`)
100000 * exp(model114$summary.fixed$`0.5quant`)

# Parámetro espacial (Barnsley)
model114$summary.random$Code$`0.5quant`[1]
exp(model114$summary.random$Code$`0.5quant`[1])

# Parámetro temporal (2002)
model114$summary.random$Year$`0.5quant`[1]
exp(model114$summary.random$Year$`0.5quant`[1])

# Parámetro de interacción (Barnsley-2002)
model114$summary.random$ID.Code$`0.5quant`[1]
exp(model114$summary.random$ID.Code$`0.5quant`[1])

# Barnsley en 2002 (manualmente y con fitted values)
100000 * exp(model114$summary.fixed$`0.5quant`) * exp(model114$summary.random$Code$`0.5quant`[1]) * exp(model114$summary.random$Year$`0.5quant`[1]) * exp(model114$summary.random$ID.Code$`0.5quant`[1])
100000 * model114$summary.fitted.values$`0.5quant`[1]

# Hiperparámetro espacial
model114$summary.hyperpar$`0.5quant`[1] # muy alto

# Hiperparámetro temporal
model114$summary.hyperpar$`0.5quant`[2]

# Hiperparámetro de interacción
model114$summary.hyperpar$`0.5quant`[3]

# Ajustamos todos los modelos con el WAIC y el DIC (iCAR, BYM, Leroux + PA1, PA2 + TipoI, TipoII + TipoIII + TipoIV)
# Tarda un rato, con onlybest = T sólo calcula el mejor (Leroux + PA2 + TipoIV)
onlybest <- T

if (onlybest) {
  
  # Leroux + PA2 + TipoIV
  form324 <- Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model324 <- inla(formula = form324, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # WAIC
  print(model324$waic$waic)
  
  # DIC
  print(model324$dic$dic)
  
  
} else {
  
  print("iCAR")
  
  # iCAR + PA1 + TipoI
  form111 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model111 <- inla(formula = form111, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # iCAR + PA1 + TipoII
  form112 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model112 <- inla(formula = form112, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # iCAR + PA1 + TipoIII
  form113 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model113 <- inla(formula = form113, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # iCAR + PA1 + TipoIV
  form114 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model114 <- inla(formula = form114, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # iCAR + PA2 + TipoI
  form121 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model121 <- inla(formula = form121, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # iCAR + PA2 + TipoII
  form122 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model122 <- inla(formula = form122, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # iCAR + PA2 + TipoIII
  form123 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model123 <- inla(formula = form123, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # iCAR + PA2 + TipoIV
  form124 <- Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model124 <- inla(formula = form124, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  print("BYM")
  
  # BYM + PA1 + TipoI
  form211 <- Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model211 <- inla(formula = form211, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # BYM + PA1 + TipoII
  form212 <- Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model212 <- inla(formula = form212, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # BYM + PA1 + TipoIII
  form213 <- Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model213 <- inla(formula = form213, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # BYM + PA1 + TipoIV
  form214 <- Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model214 <- inla(formula = form214, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # BYM + PA2 + TipoI
  form221 <- Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model221 <- inla(formula = form221, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # BYM + PA2 + TipoII
  form222 <- Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model222 <- inla(formula = form222, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # BYM + PA2 + TipoIII
  form223 <- Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model223 <- inla(formula = form223, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # BYM + PA2 + TipoIV
  form224 <- Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model224 <- inla(formula = form224, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  print("Leroux")
  
  # Leroux + PA1 + TipoI
  form311 <- Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model311 <- inla(formula = form311, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # Leroux + PA1 + TipoII
  form312 <- Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model312 <- inla(formula = form312, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # Leroux + PA1 + TipoIII
  form313 <- Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model313 <- inla(formula = form313, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # Leroux + PA1 + TipoIV
  form314 <- Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model314 <- inla(formula = form314, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # Leroux + PA2 + TipoI
  form321 <- Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model321 <- inla(formula = form321, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # Leroux + PA2 + TipoII
  form322 <- Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model322 <- inla(formula = form322, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # Leroux + PA2 + TipoIII
  form323 <- Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid"))
  model323 <- inla(formula = form323, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # Leroux + PA2 + TipoIV
  form324 <- Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1"))
  model324 <- inla(formula = form324, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = inc_lung_f)
  
  # WAIC
  model111$waic$waic
  model112$waic$waic
  model113$waic$waic
  model114$waic$waic
  model121$waic$waic
  model122$waic$waic
  model123$waic$waic
  model124$waic$waic
  model211$waic$waic
  model212$waic$waic
  model213$waic$waic
  model214$waic$waic
  model221$waic$waic
  model222$waic$waic
  model223$waic$waic
  model224$waic$waic
  model311$waic$waic
  model312$waic$waic
  model313$waic$waic
  model314$waic$waic
  model321$waic$waic
  model322$waic$waic
  model323$waic$waic
  model324$waic$waic # menor
  
  # DIC
  model111$dic$dic
  model112$dic$dic
  model113$dic$dic
  model114$dic$dic
  model121$dic$dic
  model122$dic$dic
  model123$dic$dic
  model124$dic$dic
  model211$dic$dic
  model212$dic$dic
  model213$dic$dic
  model214$dic$dic
  model221$dic$dic
  model222$dic$dic
  model223$dic$dic
  model224$dic$dic
  model311$dic$dic
  model312$dic$dic
  model313$dic$dic
  model314$dic$dic
  model321$dic$dic
  model322$dic$dic
  model323$dic$dic
  model324$dic$dic # menor
  
}

# Añadimos a la cartografía la TA de los 18 años (no me aclaro con el orden de los fitted values, lo he hecho sumando parámetros)
carto <- carto_gb
model <- model324
carto$Año <- 2002
aux <- data.frame(aux_code$Code, 100000 * exp(model$summary.fixed$`0.5quant` + model$summary.random$Code$`0.5quant`[1:142] + model$summary.random$Year$`0.5quant`[1] + model$summary.random$ID.Code$`0.5quant`[1:142]))
names(aux) <- c("Code", "TA")
carto <- merge(carto, aux)

for (i in 2:18) {
  aux1 <- carto_gb
  aux1$Año <- 2002 + i - 1
  aux <- data.frame(aux_code$Code, 100000 * exp(model$summary.fixed$`0.5quant` + model$summary.random$Code$`0.5quant`[1:142] + model$summary.random$Year$`0.5quant`[i] + model$summary.random$ID.Code$`0.5quant`[(1 + 142 * (i - 1)):(142 + 142 * (i - 1))]))
  names(aux) <- c("Code", "TA")
  aux1 <- merge(aux1, aux)
  carto <- rbind(carto, aux1)
}

# Representamos los 18 mapas y los guardamos en la carpeta de resultados
archivo <- paste(ruta_resultados, "Mapas/Ejemplo_SP_TP_All.jpeg", sep = "")
jpeg(file = archivo, width = 7000, height = 3000)
print(tm_shape(carto) + 
        tm_polygons("TA", style = "equal", n = 7, palette = paleta, border.col = "black", legend.reverse = T) + 
        tm_layout(scale = 5, title = "", main.title = "", legend.outside = T, panel.labels = as.character(2002:2019)) +
        tm_facets("Año", nrow = 3, ncol = 6))
dev.off()

# Habría que indicarle los breaks manualmente
# Se podrían hacer de 3 en 3 o de 6 en 6 años

#------------------------------------------------------------------------------#
# INCIDENCIA Y MORTALIDAD MASCULINA Y FEMENINA PULMÓN (2002-2019)
#------------------------------------------------------------------------------#

# Ya hemos ajustado incidencia femenina, falta ajustar mortalidad masculina y femenina e incidencia masculina

# Vamos a crear una tabla para guardar los indicadores de todos los modelos
models <- c("iCAR + PA1 + TipoI", "iCAR + PA1 + TipoII", "iCAR + PA1 + TipoIII", "iCAR + PA1 + TipoIV",
            "iCAR + PA2 + TipoI", "iCAR + PA2 + TipoII", "iCAR + PA2 + TipoIII", "iCAR + PA2 + TipoIV",
            "BYM + PA1 + TipoI", "BYM + PA1 + TipoII", "BYM + PA1 + TipoIII", "BYM + PA1 + TipoIV",
            "BYM + PA2 + TipoI", "BYM + PA2 + TipoII", "BYM + PA2 + TipoIII", "BYM + PA2 + TipoIV",
            "Leroux + PA1 + TipoI", "Leroux + PA1 + TipoII", "Leroux + PA1 + TipoIII", "Leroux + PA1 + TipoIV",
            "Leroux + PA2 + TipoI", "Leroux + PA2 + TipoII", "Leroux + PA2 + TipoIII", "Leroux + PA2 + TipoIV")
n_models <- length(models)
indic <- data.frame(rep(models, 4),
                    rep(rep(c("Incidencia", "Mortalidad"), each = n_models), 2),
                    rep(c("Mujeres", "Hombres"), each = 2 * n_models))
names(indic) <- c("Modelo", "Datos", "Sexo")
indic$WAIC <- 0
indic$DIC <- 0

# Creamos una lista con todas las fórmulas
form_list <- list(Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "besag", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "bym", values = codes, graph = adj_gb) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw1", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "iid", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")),
                  Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "iid")),
                  Obs ~ 1 + f(Code, model = "generic1", values = codes, Cmatrix = C) + f(Year, model = "rw2", values = years) + f(ID.Code, model = "besag", values = 1:142, graph = adj_gb, group = ID.Year, control.group = list(model = "rw1")))

# Tarda mucho, con skip = T carga directamente la tabla de indicadores ya calculada
skip <- T

if (skip) {
  
  # Cargamos la tabla ya calculada
  indic <- read.csv(paste(ruta_resultados, "indic_final.csv", sep = ""), header = T)
  
  # Seleccionamos los mejores modelos puntuando aquellos que tengan un menor WAIC y DIC. Hay contradicciones,
  # así que seleccionamos el que tenga una suma menor del WAIC y el DIC
  indic$Suma <- indic$WAIC + indic$DIC
  indic$Puntos <- 0
  for (i in 1:(dim(indic)[1] / n_models)) {
    indic$Puntos[n_models * (i - 1) + which(indic$Suma[(1 + n_models * (i - 1)):(n_models + n_models * (i - 1))] == min(indic$Suma[(1 + n_models * (i - 1)):(n_models + n_models * (i - 1))]))] <- 1
  }
  
  indic_best <- indic[which(indic$Puntos == 1), -c(6, 7)]
  indic_best
  
  # Ajustamos finalmente los modelos seleccionados
  
  # Incidencia femenina
  form_inc_f <- form_list[[which(models == indic_best$Modelo[1])]]
  model_inc_f <- inla(formula = form_inc_f, family = "poisson", E = Population, data = inc_lung_f)
  
  # Mortalidad femenina
  form_mor_f <- form_list[[which(models == indic_best$Modelo[2])]]
  model_mor_f <- inla(formula = form_mor_f, family = "poisson", E = Population, data = mor_lung_f)
  
  # Incidencia masculina
  form_inc_m <- form_list[[which(models == indic_best$Modelo[3])]]
  model_inc_m <- inla(formula = form_inc_m, family = "poisson", E = Population, data = inc_lung_m)
  
  # Mortalidad masculina
  form_mor_m <- form_list[[which(models == indic_best$Modelo[4])]]
  model_mor_m <- inla(formula = form_mor_m, family = "poisson", E = Population, data = mor_lung_m)
  
} else {
  
  # Añadimos a la tabla los datos que ya tenemos de la sección anterior
  lista <- list(model111, model112, model113, model114, model121, model122, model123, model124,
                model211, model212, model213, model214, model221, model222, model223, model224,
                model311, model312, model313, model314, model321, model322, model323, model324)
  for (i in 1:n_models) {
    indic$WAIC[i] <- lista[[i]]$waic$waic
    indic$DIC[i] <- lista[[i]]$dic$dic
  }
  
  # Calculamos todos los modelos restantes, lo hacemos dentro de un for para simplificarlo
  data_list <- list(mor_lung_f, inc_lung_m, mor_lung_m)
  
  cont <- n_models + 1
  for (i in 1:3) {
    datos <- data_list[[i]]
    for (j in 1:n_models) {
      form <- form_list[[j]]
      model <- inla(formula = form, family = "poisson", E = Population, control.compute = list(waic = T, dic = T), data = datos)
      indic$WAIC[cont] <- model$waic$waic
      indic$DIC[cont] <- model$dic$dic
      cont <- cont + 1
    }
  }
  
  # Guardamos la tabla con los indicadores
  write.csv(indic, "C:/Users/Ander/Documents/TFM/Resultados/indic_final.csv", row.names = F)
  
  # Seleccionamos los mejores modelos puntuando aquellos que tengan un menor WAIC y DIC. Hay contradicciones,
  # así que seleccionamos el que tenga una suma menor del WAIC y el DIC
  indic$Suma <- indic$WAIC + indic$DIC
  indic$Puntos <- 0
  for (i in 1:(dim(indic)[1] / n_models)) {
    indic$Puntos[n_models * (i - 1) + which(indic$Suma[(1 + n_models * (i - 1)):(n_models + n_models * (i - 1))] == min(indic$Suma[(1 + n_models * (i - 1)):(n_models + n_models * (i - 1))]))] <- 1
  }
  
  indic_best <- indic[which(indic$Puntos == 1), -c(6, 7)]
  
  # Ajustamos finalmente los modelos seleccionados
  
  # Incidencia femenina
  form_inc_f <- form_list[[which(models == indic_best$Modelo[1])]]
  model_inc_f <- inla(formula = form_inc_f, family = "poisson", E = Population, data = inc_lung_f)
  
  # Mortalidad femenina
  form_mor_f <- form_list[[which(models == indic_best$Modelo[2])]]
  model_mor_f <- inla(formula = form_mor_f, family = "poisson", E = Population, data = mor_lung_f)
  
  # Incidencia masculina
  form_inc_m <- form_list[[which(models == indic_best$Modelo[3])]]
  model_inc_m <- inla(formula = form_inc_m, family = "poisson", E = Population, data = inc_lung_m)
  
  # Mortalidad masculina
  form_mor_m <- form_list[[which(models == indic_best$Modelo[4])]]
  model_mor_m <- inla(formula = form_mor_m, family = "poisson", E = Population, data = mor_lung_m)
  
}


