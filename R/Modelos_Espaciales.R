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
cancer <- "Pulmón"

# Seleccionar sexo (para Cérvix y Próstata está predefinido):
# "Hombres", "Mujeres"
sexo <- "Mujeres"

# Seleccionar tipo de dato:
# "Incidencia", "Mortalidad"
tipo <- "Incidencia"

# Seleccionar año para el análisis espacial:
# uno entre 2002 y 2019
año <- 2019

# Seleccionar estrategia para los modelos:
# "gaussian", "simplified.laplace", "laplace"
estrategia <- "laplace"

# Seleccionar Tasa Cruda para las probabilidades de exceso
# (en https://gco.iarc.fr/ se pueden encontrar tasas a nivel europeo,
# si se deja en NULL se calcula automáticamente una TC media general)
umbral <- NULL

#------------------------------------------------------------------------------#
# 2) DATOS
#------------------------------------------------------------------------------#

# Leemos y arreglamos los datos seleccionados
datos <- select_data(cancer, sexo, tipo)
datosAnyo <- datos[which(datos$Year == año), -2]
datosAnyo$TC <- 100000 * datosAnyo$Obs / datosAnyo$Population
if (is.null(umbral)) {umbral <- mean(datosAnyo$TC)}

# Leemos y ordenamos la cartografía de Gran Bretaña
carto_gb <- st_read("../Datos/Carto/carto_gb.shp")
carto_gb <- carto_gb[order(carto_gb$Code), ]

# Leemos la matriz de adyacencia (grafo conexo)
adj_gb <- as.matrix(read.table("../Datos/Carto/adj_gb.txt"))

# Representamos y guardamos el grafo de vecindad en la carpeta de Figuras
carto.nb <- mat2listw(adj_gb)$neighbours
if(n.comp.nb(carto.nb)$nc == 1) {print("El grafo es conexo")} else {"El grafo no es conexo"}

archivo <- "Figuras/Grafo.png"
graphics.off()
png(file = archivo, width = 1200, height = 2000)
plot(carto_gb$geometry, col = "grey90")
plot(carto.nb, st_centroid(st_geometry(carto_gb), of_largest_polygon = TRUE),
     pch = 19, cex = 0.5, col = "blue", add = TRUE)
dev.off()

#------------------------------------------------------------------------------#
# 3) ANÁLISIS DESCRIPTIVO
#------------------------------------------------------------------------------#
carto <- merge(carto_gb, datosAnyo)

# Creamos los tres mapas y los guardamos en la carpeta de Figuras Espaciales
n.color <- 7
paleta <-  get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)

breaks.obs <- round(seq(min(carto$Obs), max(carto$Obs), length.out = n.color + 1))
breaks.pop <- round(seq(min(carto$Population), max(carto$Population), length.out = n.color + 1))
breaks.TC <- round(seq(min(carto$TC), max(carto$TC), length.out = n.color + 1))

breaks.TC[1] <- floor(min(carto$TC))
breaks.TC[n.color + 1] <- ceiling(max(carto$TC))

names(carto)[3:4] <- c("Población", "Casos")

archivo <- paste(paste(paste(paste(paste("Figuras/Espacial/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), año, sep = "_"), ".jpeg", sep = "")
graphics.off()
jpeg(file = archivo, width = 2000, height = 1200)
print(tm_shape(carto) + 
        tm_polygons(c("Casos","Población","TC"), style = "fixed", n = n.color,
                    breaks = list(breaks.obs, breaks.pop, breaks.TC),
                    lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 2, legend.text.size = 1.5, legend.position = c("left", "top")))
dev.off()

names(carto)[3:4] <- c("Population", "Obs")

#------------------------------------------------------------------------------#
# 4) MODELOS ESPACIALES 
#------------------------------------------------------------------------------#

# Ajustamos el modelo iCAR
model1 <- CAR_INLA(carto = carto, ID.area = "Code", O = "Obs", E = "Population", W = adj_gb,
                   prior = "intrinsic", model = "global", strategy = estrategia)

# Representamos las distribuciones marginales a posteriori del valor base, del parámetro espacial correspondiente
# a la primera región (Barnsley) y del hiperparámetro. Lo guardamos en la carpeta de Figuras Espaciales
archivo <- paste(paste(paste(paste(paste("Figuras/Espacial/Posteriores", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), año, sep = "_"), ".pdf", sep = "")
graphics.off()
pdf(file = archivo, width = 10, height = 5)
par(mfrow = c(1, 3))
plot(inla.smarginal(model1$marginals.fixed$`(Intercept)`), type = "l", lwd = 1.5, main = "Valor base", xlab = expression(eta), ylab = expression(tilde(p)(paste(eta, "|", y))))
abline(v = model1$summary.fixed$`0.5quant`, col = "blue")
plot(inla.smarginal(model1$marginals.random$ID.area$index.1), type = "l", lwd = 1.5, main = "Parámetro espacial", xlab = expression(u[1]), ylab = expression(tilde(p)(paste(u[1], "|", y))))
abline(v = model1$summary.random$ID.area$`0.5quant`[1], col = "blue")
plot(inla.smarginal(model1$marginals.hyperpar$`Precision for ID.area`), type = "l", lwd = 1.5, main = "Hiperparámetro", xlab = expression(tau[u]), ylab = expression(tilde(p)(paste(tau[u], "|", y))))
abline(v = model1$summary.hyperpar$`0.5quant`, col = "blue")
par(mfrow = c(1, 1))
dev.off()

# Ajustamos todos los modelos espaciales

# i.i.d.
form0 <- Obs ~ f(Code, model = "iid", values = unique(datosAnyo$Code), graph = adj_gb)
model0 <- inla(formula = form0, family = "poisson", data = datosAnyo, E = Population,
               control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE, return.marginals.predictor = TRUE),
               control.inla = list(strategy = estrategia))

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

Tabla <- data.frame(DIC = unlist(lapply(model_list, function(x) x$dic$dic)),
                    WAIC = unlist(lapply(model_list, function(x) x$waic$waic)),
                    LS = unlist(lapply(model_list, function(x) -sum(log(x$cpo$cpo)))))
print(Tabla)

#------------------------------------------------------------------------------#
# 5) REPRESENTACIONES GRÁFICAS
#------------------------------------------------------------------------------#

# Escogemos el de Leroux
model <- model3

# Representamos la distribución de lambda y el área bajo la curva de los fitted values
den <- as.data.frame(100000 * model$marginals.fitted.values$fitted.Predictor.001)

archivo <- paste(paste(paste(paste(paste("Figuras/Espacial/Probabilities", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), año, sep = "_"), ".pdf", sep = "")
graphics.off()
pdf(file = archivo, width = 10, height = 5)
par(mfrow = c(1, 2))
plot(model$marginals.hyperpar$`Beta for ID.area`, type = "l", lwd = 1.5, main = "Hiperparámetro de Leroux", xlab = expression(lambda), ylab = expression(tilde(p)(paste(lambda, "|", y))))
abline(v = model$summary.hyperpar$`0.5quant`[2], col = "blue")
plot(den, type = "l", lwd = 1.5, main = "Predictor lineal", xlab = expression(eta[1]), ylab = expression(tilde(p)(paste(eta[1], "|", y))))
polygon(c(den$x[den$x >= umbral ], umbral),
        c(den$y[den$x >= umbral ], 0),
        col = "skyblue",
        border = 1)
par(mfrow = c(1, 1))
dev.off()

# Representamos en un mapa las tasas ajustadas y la probabilidad de que supere el umbral establecido. Lo guardamos
# en la carpeta de Figuras Espaciales
carto$TA <- 100000 * model$summary.fitted.values$`0.5quant`
carto$Prob <- unlist(lapply(model$marginals.fitted.values, function(x) 1 - inla.pmarginal(umbral / 100000, x)))

# Creamos los mapas y los guardamos en la carpeta de Figuras Espaciales
n.color <- 7
paleta <- get_brewer_pal("-RdYlGn", n = n.color, contrast = 0.8, plot = F)
values <- 100000 * seq(min(model$summary.fitted.values$mean), max(model$summary.fitted.values$mean), length.out = n.color + 1)
values[1] <- floor(values[1])
values[n.color + 1] <- ceiling(values[n.color + 1])
values <- round(values)

archivo <- paste(paste(paste(paste(paste("Figuras/Espacial/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), año, sep = "_"), "TA.jpeg", sep = "_")
graphics.off()
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("TA", style = "fixed", breaks = values, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()

n.color <- 7
paleta <- get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)

archivo <- paste(paste(paste(paste(paste("Figuras/Espacial/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), año, sep = "_"), "Prob.jpeg", sep = "_")
graphics.off()
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("Prob", style = "fixed", interval.closure = "left", legend.reverse = T,
                    breaks = c(0, 0.05, 0.1, 0.2, 0.8, 0.9, 0.95, 1), labels = c("[0-0.05)", "[0.05-0.1)", "[0.1-0.2)", "[0.2-0.8)", "[0.8-0.9)", "[0.9-0.95)", "[0.95-1]"),
                    lwd = 2, border.col = "black", palette = paleta) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()
