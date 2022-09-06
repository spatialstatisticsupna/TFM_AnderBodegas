rm(list = ls())

library(INLA)
library(sf)
library(spdep)
library(tmap)
library(tmaptools)
library(bigDM)

#------------------------------------------------------------------------------#
# 0) FUNCIONES
#------------------------------------------------------------------------------#

select_data <- function(cancer, sexo, tipo) {
  
  if (cancer == "Cérvix") {sexo <- "Mujeres"}
  if (cancer == "Próstata") {sexo <- "Hombres"}
  
  if (sexo == "Mujeres") {
    
    Datos_gb_f <- read.csv("../Datos/Datos_gb_f.csv", header = T)
    
    cancer_list <- c("Todos", "Leucemia", "Mama", "Cérvix", "Melanoma", "Vejiga", "Colorrectal", "Páncreas", "Estómago", "Pulmón", "Hígado", "Esófago")
    tipo_list <- c("Incidencia", "Mortalidad")
    
    index <- 3 + which(tipo_list == tipo) + 2 * (which(cancer_list == cancer) - 1)
    
    datos <- Datos_gb_f[, c(1:3, index)]
    names(datos)[4] <- "Obs"
    datos <- datos[order(datos$Year, datos$Code), ]
    
  }
  
  if (sexo == "Hombres") {
    
    Datos_gb_m <- read.csv("../Datos/Datos_gb_m.csv", header = T)
    
    cancer_list <- c("Todos", "Leucemia", "Mama", "Próstata", "Melanoma", "Vejiga", "Colorrectal", "Páncreas", "Estómago", "Pulmón", "Hígado", "Esófago")
    tipo_list <- c("Incidencia", "Mortalidad")
    
    index <- 3 + which(tipo_list == tipo) + 2 * (which(cancer_list == cancer) - 1)
    
    datos <- Datos_gb_m[, c(1:3, index)]
    names(datos)[4] <- "Obs"
    datos <- datos[order(datos$Year, datos$Code), ]
    
  }
  
  return(datos)
  
}

#------------------------------------------------------------------------------#
# 1) OPCIONES
#------------------------------------------------------------------------------#

# Seleccionar tipo de cáncer:
# "Todos", "Leucemia", "Mama", "Melanoma", "Vejiga", "Colorrectal",
# "Páncreas", "Estómago", "Pulmón", "Hígado", "Esófago", "Cérvix", "Próstata"
cancer <- "Melanoma"

# Seleccionar sexo (para Cérvix y Próstata no es necesario):
# "Hombres", "Mujeres"
sexo <- "Hombres"

# Seleccionar tipo de dato:
# "Incidencia", "Mortalidad"
tipo <- "Incidencia"

# Seleccionar año para el análsis espacial:
# uno entre 2002 y 2019
año <- 2019

# Seleccionar estrategia para los modelos:
# "gaussian", "simplified.laplace", "laplace"
estrategia <- "gaussian"

# Seleccionar Tasa Cruda para exeedance probabilities
# TC 2020 Europa Cáncer de pulmón masculino: 87.1	(https://gco.iarc.fr/)
# TC 2020 Europa Cáncer de pulmón femenino:  42.0	(https://gco.iarc.fr/)
# TC 2020 Europa Melanoma masculino:         21.1	(https://gco.iarc.fr/)
# TC 2020 Europa Melanoma femenino:          19.2	(https://gco.iarc.fr/)
umbral <- 21.1

#------------------------------------------------------------------------------#
# 2) LECTURA DE DATOS
#------------------------------------------------------------------------------#

# Leemos los datos seleccionados
datos <- select_data(cancer, sexo, tipo)
datosAño <- datos[which(datos$Year == año), ]

# Leemos la cartografía de Gran Bretaña
carto_gb <- st_read("../Datos/Carto/carto_gb.shp")
head(carto_gb)

# Leemos la matriz de adyacencia (grafo conexo)
adj_gb <- as.matrix(read.table("../Datos/Carto/adj_gb.txt"))

# Representamos el grafo de vecindad
carto.nb <- mat2listw(adj_gb)$neighbours
if(n.comp.nb(carto.nb)$nc == 1) {print("El grafo es conexo")} else {"El grafo no es conexo"}

plot(carto_gb$geometry, main = "Grafo de vecindad modificado (conexo)")
plot(carto.nb, st_centroid(st_geometry(carto_gb), of_largest_polygon = TRUE),
     pch = 19, cex = 0.5, col = "blue", add = TRUE)

#------------------------------------------------------------------------------#
# 3) TASA CRUDA: Casos/Población x 100.000
#------------------------------------------------------------------------------#
datosAño$TC <- 100000 * datosAño$Obs / datosAño$Population
carto <- merge(carto_gb, datosAño)

# Creamos los tres mapas y los guardamos en la carpeta de figuras
archivo <- paste(paste(paste(paste("Figuras/Espacial/Mapa_Cancer", cancer, sep = "_"), sexo, sep = "_"), año, sep = "_"), ".jpeg", sep = "")
n.color <- 7
paleta <-  get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)
breaks.obs <- round(seq(min(carto$Obs), max(carto$Obs), length.out = n.color + 1))
breaks.pop <- round(seq(min(carto$Population), max(carto$Population), length.out = n.color + 1))
breaks.TC <- round(seq(min(carto$TC), max(carto$TC), length.out = n.color + 1))
breaks.TC[1] <- floor(min(carto$TC))
breaks.TC[n.color + 1] <- ceiling(max(carto$TC))

graphics.off()
jpeg(file = archivo, width = 2000, height = 1200)
print(tm_shape(carto) + 
        tm_polygons(c("Obs","Population","TC"), style = "fixed", n = n.color,
                    breaks = list(breaks.obs, breaks.pop,breaks.TC),
                    lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 2, legend.text.size = 1.5, legend.position = c("left", "top")))
dev.off()

#------------------------------------------------------------------------------#
# 4) MODELOS ESPACIALES 
#------------------------------------------------------------------------------#

# Ajustamos el modelo iCAR
model1 <- CAR_INLA(carto = carto, ID.area = "Code", O = "Obs", E = "Population", W = adj_gb,
                   prior = "intrinsic", model = "global", strategy = estrategia)
summary(model1)

# Representamos las distribuciones marginales a posteriori del valor base, del parámetro espacial correspondiente
# a la primera región (Barnsley) y del hiperparámetro. Lo guardamos en la carpeta de figuras
archivo <- paste(paste(paste(paste("Figuras/Espacial/Posteriores", cancer, sep = "_"), sexo, sep = "_"), año, sep = "_"), "S.pdf", sep = "_")

graphics.off()
pdf(file = archivo, width = 10, height = 5)
par(mfrow = c(1, 3))
plot(inla.smarginal(model1$marginals.fixed$`(Intercept)`), type = "l", lwd = 1.5, main = "Valor base", xlab = expression(eta), ylab = expression(tilde(p)(paste(eta, "|", y))))
abline(v = model1$summary.fixed$`0.5quant`, col = "blue")
plot(inla.smarginal(model1$marginals.random$ID.area$index.1), type = "l", lwd = 1.5, main = "Parámetro espacial", xlab = expression(u[1]), ylab = expression(tilde(p)(paste(u[1], "|", y))))
abline(v = model1$summary.random$ID.area$`0.5quant`[1], col = "blue")
plot(inla.smarginal(model1$marginals.hyperpar$`Precision for ID.area`), type = "l", lwd = 1.5, main = "Hiperparámetro", xlab = expression(tau), ylab = expression(tilde(p)(paste(tau, "|", y))))
abline(v = model1$summary.hyperpar$`0.5quant`, col = "blue")
par(mfrow = c(1, 1))
dev.off()

# Valor base
model1$summary.fixed$`0.5quant`
inla.emarginal(function(x) 100000 * exp(x), model1$marginals.fixed$`(Intercept)`)

# Parámetro espacial de Barnsley
model1$summary.random$Code$`0.5quant`[1]
inla.emarginal(function(x) exp(x), model1$marginals.random$ID.area[[1]])

# Barnsley (manualmente y con fitted values)
100000 * exp(model1$summary.fixed$`0.5quant`) * exp(model1$summary.random$ID.area$`0.5quant`[1])
100000 * model1$summary.fitted.values$`0.5quant`[1]

# Distribuciones a posteriori de los hiperparámetros
model1$summary.hyperpar

# Ajustamos todos los modelos espaciales

# i.i.d.
form0 <- Obs ~ f(Code, model = "iid", values = unique(datosAño$Code), graph = adj_gb)
model0 <- inla(formula = form0, family = "poisson", data = datosAño, E = Population,
               control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE),
               control.inla = list(strategy = "laplace"))

# iCAR
model1 <- CAR_INLA(carto = carto, ID.area = "Code", O = "Obs", E = "Population", W = adj_gb,
                   prior = "intrinsic", model = "global", strategy = estrategia)

# BYM
model2 <- CAR_INLA(carto = carto, ID.area = "Code", O = "Obs", E = "Population", W = adj_gb,
                   prior = "BYM", model = "global", strategy = estrategia)

# Leroux
model3 <- CAR_INLA(carto = carto, ID.area = "Code", O = "Obs", E = "Population", W = adj_gb,
                   prior = "Leroux", model = "global", strategy = estrategia)

# Comparamos los modelos
model_list <- list(iid = model0, iCAR = model1, BYM = model2, Leroux = model3)

Tabla <- data.frame(mean.deviance = unlist(lapply(model_list, function(x) x$dic$mean.deviance)),
                    p.eff = unlist(lapply(model_list, function(x) x$dic$p.eff)),
                    DIC = unlist(lapply(model_list, function(x) x$dic$dic)),
                    WAIC = unlist(lapply(model_list, function(x) x$waic$waic)),
                    LS = unlist(lapply(model_list, function(x) -sum(log(x$cpo$cpo)))))
print(Tabla)

#------------------------------------------------------------------------------#
# 5) REPRESENTACIONES GRÁFICAS
#------------------------------------------------------------------------------#

# Escogemos el de Leroux
model <- model3

# Representamos en un mapa las tasas ajustadas y la probabilidad de que supere a la tasa media de GB
tasa.global <- exp(model$summary.fixed$`0.5quant`)

carto$TA <- 100000 * model$summary.fitted.values$`0.5quant`
carto$Prob <- unlist(lapply(model$marginals.fitted.values, function(x) 1 - inla.pmarginal(umbral / 100000, x)))

# Creamos los mapas y los guardamos en la carpeta figuras
archivo <- paste(paste(paste(paste("Figuras/Espacial/Mapa_Cancer", cancer, sep = "_"), sexo, sep = "_"), año, sep = "_"), "TA.jpeg", sep = "_")

n.color <- 7
paleta <- get_brewer_pal("-RdYlGn", n = n.color, contrast = 0.8, plot = F)
values <- 100000 * seq(min(model$summary.fitted.values$mean), max(model$summary.fitted.values$mean),length.out = n.color + 1)
values[1] <- floor(values[1])
values[n.color + 1] <- ceiling(values[n.color + 1])
values <- round(values)

graphics.off()
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("TA", style = "fixed", breaks = values, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()

archivo <- paste(paste(paste(paste("Figuras/Espacial/Mapa_Cancer", cancer, sep = "_"), sexo, sep = "_"), año, sep = "_"), "Prob.jpeg", sep = "_")
n.color <- 5
paleta <- get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)
                            
graphics.off()
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("Prob", style = "fixed", interval.closure = "left", legend.reverse = T,
                    breaks = c(0, 0.1, 0.2, 0.8, 0.9, 1), labels = c("[0-0.1)", "[0.1-0.2)", "[0.2-0.8)", "[0.8-0.9)", "[0.9-1]"),
                    lwd = 2, border.col = "black", palette = paleta) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()

