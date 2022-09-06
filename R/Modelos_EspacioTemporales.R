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

# Seleccionar sexo (para Cérvix y Próstata no es necesario):
# "Hombres", "Mujeres"
sexo <- "Hombres"

# Seleccionar tipo de dato:
# "Incidencia", "Mortalidad"
tipo <- "Incidencia"

# Seleccionar estrategia para los modelos:
# "gaussian", "simplified.laplace", "laplace"
estrategia <- "gaussian"

# Seleccionar Tasa Cruda para exeedance probabilities
# TC Incidencia 2020 Europa Cáncer de pulmón masculino: 87.1	(https://gco.iarc.fr/)
# TC Incidencia 2020 Europa Cáncer de pulmón femenino:  42.0	(https://gco.iarc.fr/)
# TC Incidencia 2020 Europa Melanoma masculino:         21.1	(https://gco.iarc.fr/)
# TC Incidencia 2020 Europa Melanoma femenino:          19.2	(https://gco.iarc.fr/)
# TC Mortalidad 2020 Europa Cáncer de pulmón masculino: 71.9	(https://gco.iarc.fr/)
# TC Mortalidad 2020 Europa Cáncer de pulmón femenino:  32.1	(https://gco.iarc.fr/)
# TC Mortalidad 2020 Europa Melanoma masculino:         4.1	(https://gco.iarc.fr/)
# TC Mortalidad 2020 Europa Melanoma femenino:          3	(https://gco.iarc.fr/)
umbral <- 87.1

#------------------------------------------------------------------------------#
# 2) LECTURA DE DATOS
#------------------------------------------------------------------------------#

# Leemos los datos seleccionados
datos <- select_data(cancer, sexo, tipo)

# Leemos la cartografía de Gran Bretaña
carto_gb <- st_read("../Datos/Carto/carto_gb.shp")
head(carto_gb)

# Leemos la matriz de adyacencia (grafo conexo)
adj_gb <- as.matrix(read.table("../Datos/Carto/adj_gb.txt"))

#------------------------------------------------------------------------------#
# 3) TASA CRUDA: Casos/Población x 100.000
#------------------------------------------------------------------------------#

# Calculamos la tasa cruda por área y año
S <- nrow(carto_gb)
T <- length(unique(datos$Year))
T.from <- min(datos$Year)
T.to <- max(datos$Year)

datos$TC <- 100000 * datos$Obs / datos$Population

carto <- do.call("rbind", replicate(18, carto_gb, simplify = F))
carto <- cbind(carto, datos[which(is.element(datos$Year, seq(T.from, T.to, 2))), c(2, 5)])
head(carto)

# Evolución geográfica de las TC por año
archivo <- paste(paste(paste("Figuras/EspacioTemporal/Mapa_Cancer", cancer, sep = "_"), sexo, sep = "_"), "2002_2019.jpeg", sep = "_")

n.color <- 7
paleta <-  get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)
breaks <- round(seq(min(datos$TC), max(datos$TC), length.out = n.color + 1))
breaks[1] <- floor(min(datos$TC))
breaks[n.color + 1] <- ceiling(max(datos$TC))

graphics.off()
jpeg(file = archivo, width = 3000, height = 3000)
print(tm_shape(carto) + 
        tm_polygons("TC", style = "fixed", n = n.color, breaks = breaks, 
                    border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(scale = 5, title = "", legend.outside = T, panel.show = T, panel.labels = seq(T.from, T.to, 2)) +
        tm_facets("Year", ncol = 3, nrow = 3))
dev.off()


# Evolución temporal de las TC por región
aux <- carto[,as.character(T.from:T.to)]
aux$geometry <- NULL

matplot(t(aux), type = "l", xlab = "Año", ylab = "Tasa cruda", xaxt = "n")
axis(1, at = seq(1, T, 2), labels = seq(T.from, T.to, 2))

#------------------------------------------------------------------------------#
# 4) MODELOS ESPACIO-TEMPORALES (Leroux + RW1)
#------------------------------------------------------------------------------#

## Leroux + RW1 (modelo aditivo) ##
RW1.Additive <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                           W = adj_gb, spatial = "Leroux", temporal = "rw1", interaction = "none",
                           model = "global", strategy = estrategia)

## Leroux + RW1 + TypeI ##
RW1.TypeI <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                        W = adj_gb, spatial = "Leroux", temporal = "rw1", interaction = "TypeI",
                        model = "global", strategy = estrategia)

## Leroux + RW1 + TypeII ##
RW1.TypeII <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                         W = adj_gb, spatial = "Leroux", temporal = "rw1", interaction = "TypeII",
                         model = "global", strategy = estrategia)

## Leroux + RW1 + TypeIII ##
RW1.TypeIII <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                          W = adj_gb, spatial = "Leroux", temporal = "rw1", interaction = "TypeIII",
                          model = "global", strategy = estrategia)

## Leroux + RW1 + TypeIV ##
RW1.TypeIV <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                         W = adj_gb, spatial = "Leroux", temporal = "rw1", interaction = "TypeIV",
                         model = "global", strategy = estrategia)

## Leroux + RW2 (modelo aditivo) ##
RW2.Additive <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                           W = adj_gb, spatial = "Leroux", temporal = "rw2", interaction = "none",
                           model = "global", strategy = estrategia)

## Leroux + RW2 + TypeI ##
RW2.TypeI <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                        W = adj_gb, spatial = "Leroux", temporal = "rw2", interaction = "TypeI",
                        model = "global", strategy = estrategia)

## Leroux + RW2 + TypeII ##
RW2.TypeII <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                         W = adj_gb, spatial = "Leroux", temporal = "rw2", interaction = "TypeII",
                         model = "global", strategy = estrategia)

## Leroux + RW2 + TypeIII ##
RW2.TypeIII <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                          W = adj_gb, spatial = "Leroux", temporal = "rw2", interaction = "TypeIII",
                          model = "global", strategy = estrategia)

## Leroux + RW2 + TypeIV ##
RW2.TypeIV <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                         W = adj_gb, spatial = "Leroux", temporal = "rw2", interaction = "TypeIV",
                         model = "global", strategy = estrategia)

# Comparamos los modelos
model_list <- list(RW1.Aditivo = RW1.Additive, RW1.TipoI = RW1.TypeI, RW1.TipoII = RW1.TypeII, RW1.TipoIIII = RW1.TypeIII, RW1.TipoIV = RW1.TypeIV,
                   RW2.Aditivo = RW2.Additive, RW2.TipoI = RW2.TypeI, RW2.TipoII = RW2.TypeII, RW2.TipoIIII = RW2.TypeIII, RW2.TipoIV = RW1.TypeIV)

Tabla <- data.frame(mean.deviance = unlist(lapply(model_list, function(x) x$dic$mean.deviance)),
                    p.eff = unlist(lapply(model_list, function(x) x$dic$p.eff)),
                    DIC = unlist(lapply(model_list, function(x) x$dic$dic)),
                    WAIC = unlist(lapply(model_list, function(x) x$waic$waic)),
                    LS = unlist(lapply(model_list, function(x) -sum(log(x$cpo$cpo)))))
print(Tabla)

#------------------------------------------------------------------------------#
# 5) REPRESENTACIONES GRÁFICAS
#------------------------------------------------------------------------------#
model <- RW1.TypeII

# Distribuciones a posteriori de los hiperparámetros y del valor base
archivo <- paste(paste(paste("Figuras/EspacioTemporal/Posteriores", cancer, sep = "_"), sexo, sep = "_"), "ST.pdf", sep = "_")

model$summary.fixed
model$summary.hyperpar

graphics.off()
pdf(file = archivo, width = 10, height = 5)
par(mfrow = c(2, 3))
plot(inla.smarginal(model$marginals.fixed$`(Intercept)`), type = "l", lwd = 1.5, main = "Valor base", xlab = expression(eta), ylab = expression(tilde(p)(paste(eta, "|", y))))
abline(v = model$summary.fixed$`0.5quant`, col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Precision for ID.area`), type = "l", lwd = 1.5, main = "Hiperparámetro espacial", xlab = expression(tau[u]), ylab = expression(tilde(p)(paste(tau[u], "|", y))))
abline(v = model$summary.hyperpar["Precision for ID.area","0.5quant"], col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Precision for ID.year`), type = "l", lwd = 1.5, main = "Hiperparámetro temporal", xlab = expression(tau[phi]), ylab = expression(tilde(p)(paste(tau[phi], "|", y))))
abline(v = model$summary.hyperpar["Precision for ID.year","0.5quant"], col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Precision for ID.area.year`), type = "l", lwd = 1.5, main = "Hiperparámetro de interacción", xlab = expression(tau[delta]), ylab = expression(tilde(p)(paste(tau[delta], "|", y))))
abline(v = model$summary.hyperpar["Precision for ID.area.year","0.5quant"], col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Beta for ID.area`), type = "l", lwd = 1.5, main = "Hiperparámetro de Leroux", xlab = expression(lambda), ylab = expression(tilde(p)(paste(lambda, "|", y))))
abline(v = model$summary.hyperpar["Beta for ID.area","0.5quant"], col = "blue")
par(mfrow = c(1, 1))
dev.off()

# Patrón espacial (común para todo el periodo)
tasa.global <- exp(model$summary.fixed$`0.5quant`)
aux <- lapply(model$marginals.random$ID.area, function(x) inla.tmarginal(function(y) 100000 * tasa.global * exp(y), x))
pattern.S <- unlist(lapply(aux, function(x) inla.qmarginal(0.5, x)))
carto$TA <- pattern.S

archivo <- paste(paste(paste("Figuras/EspacioTemporal/Mapa_Cancer", cancer, sep = "_"), sexo, sep = "_"), "PatrónEspacial.jpeg", sep = "_")
n.color <- 7
paleta <- get_brewer_pal("-RdYlGn", n = n.color, contrast = 0.8, plot = F)

graphics.off()
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("TA", style = "equal", n = n.color, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()

# Patrón temporal (común para todas las áreas)
archivo <- paste(paste(paste("Figuras/EspacioTemporal/Grafica_Cancer", cancer, sep = "_"), sexo, sep = "_"), "PatrónTemporal.pdf", sep = "_")

aux <- lapply(model$marginals.random$ID.year, function(x) inla.tmarginal(function(y) 100000 * tasa.global * exp(y), x))
temporal.risk <- unlist(lapply(aux, function(x) inla.emarginal(function(y) y, x)))
q1 <- unlist(lapply(aux, function(x) inla.qmarginal(0.025, x)))
q2 <- unlist(lapply(aux, function(x) inla.qmarginal(0.5, x)))
q3 <- unlist(lapply(aux, function(x) inla.qmarginal(0.975, x)))

graphics.off()
pdf(file = archivo)
x <- 1:T
plot(temporal.risk, type = "l", ylim = range(temporal.risk) * c(1 / 1.1, 1.1),
     xaxt = "n", xlab = "Año", ylab = "Casos x 100,000 habitantes", main = "Patrón temporal")
axis(1, at = seq(1, T, 2), labels = seq(T.from, T.to, 2))
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1, tail(q3, 1), rev(q3), q1[1])
polygon(X.Vec, Y.Vec, col = "gray", border = NA)
lines(q2)
abline(h = 100000 * tasa.global, lty = 2, col = "blue")
dev.off()

# Mapas con las tasas ajustadas por año
TA <- data.frame(rep(2002:2019, each = 142), 100000 * model$summary.fitted.values$`0.5quant`)
names(TA) <- c("Year", "TA")
carto <- do.call("rbind", replicate(18, carto_gb, simplify = F))
carto <- cbind(carto, TA[which(is.element(TA$Year, seq(T.from, T.to, 2))), ])
head(carto)

archivo <- paste(paste(paste("Figuras/EspacioTemporal/Mapa_Cancer", cancer, sep = "_"), sexo, sep = "_"), "TA.jpeg", sep = "_")

n.color <- 7
paleta <- get_brewer_pal("-RdYlGn", n = n.color, contrast = 0.8, plot = F)
breaks <- 100000 * seq(min(model$summary.fitted.values$`0.5quant`), max(model$summary.fitted.values$`0.5quant`), length.out = n.color + 1)
breaks[1] <- floor(breaks[1])
breaks[n.color + 1] <- floor(breaks[n.color + 1])
breaks <- round(breaks)

graphics.off()
jpeg(file = archivo, width = 3000, height = 3000)
print(tm_shape(carto) + 
        tm_polygons("TA", style = "fixed", n = n.color, breaks = breaks, 
                    border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(scale = 5, title = "", legend.outside = T, panel.show = T, panel.labels = seq(T.from, T.to, 2)) +
        tm_facets("Year", ncol = 3, nrow = 3))
dev.off()


# Mapas con las probabilidades a posteriori
Prob <- unlist(lapply(model$marginals.fitted.values, function(x) 1 - inla.pmarginal(umbral / 100000, x)))
Prob <- data.frame(rep(2002:2019, each = 142), Prob)
names(Prob) <- c("Year", "Prob")
carto <- do.call("rbind", replicate(18, carto_gb, simplify = F))
carto <- cbind(carto, Prob[which(is.element(TA$Year, seq(T.from, T.to, 2))), ])
head(carto)

archivo <- paste(paste(paste("Figuras/EspacioTemporal/Mapa_Cancer", cancer, sep = "_"), sexo, sep = "_"), "Prob.jpeg", sep = "_")

n.color <- 5
paleta <- get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)

graphics.off()
jpeg(file = archivo, width = 3000, height = 3000)
print(tm_shape(carto) + 
        tm_polygons("Prob", style = "fixed", n = n.color, interval.closure = "left", 
                    breaks = c(0, 0.1, 0.2, 0.8, 0.9, 1), labels = c("[0-0.1)", "[0.1-0.2)", "[0.2-0.8)", "[0.8-0.9)", "[0.9-1]"),
                    border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(scale = 5, title = "", legend.outside = T, panel.show = T, panel.labels = seq(T.from, T.to, 2)) +
        tm_facets("Year", ncol = 3, nrow = 3))
dev.off()
