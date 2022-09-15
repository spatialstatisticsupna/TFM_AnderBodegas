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

# Seleccionar estrategia para los modelos:
# "gaussian", "simplified.laplace", "laplace"
estrategia <- "gaussian"

# Seleccionar Tasa Cruda para las probabilidades de exceso
# (en https://gco.iarc.fr/ se pueden encontrar tasas a nivel europeo,
# si se deja en NULL se calcula automáticamente una TC media general)
umbral <- NULL

# ¿Saltar el ajuste de los 10 modelos y trabajar con Leroux + PA1 + TipoII?
# (TRUE = Saltar, FALSE = No saltar)
saltar <- TRUE

# ¿Representar los mapas todos juntos o por separado?
# (TRUE = juntos, FALSE = por separado)
juntos <- TRUE

#------------------------------------------------------------------------------#
# 2) LECTURA DE DATOS
#------------------------------------------------------------------------------#

# Leemos y arreglamos los datos seleccionados
datos <- select_data(cancer, sexo, tipo)
datos$TC <- 100000 * datos$Obs / datos$Population
if (is.null(umbral)) {
  umbral <- aggregate(x = datos$TC, by = list(datos$Year), FUN = mean)
  umbral <- mean(umbral$x)
}

# Leemos y ordenamos la cartografía de Gran Bretaña
carto_gb <- st_read("../Datos/Carto/carto_gb.shp")
carto_gb <- carto_gb[order(carto_gb$Code), ]

# Leemos la matriz de adyacencia (grafo conexo)
adj_gb <- as.matrix(read.table("../Datos/Carto/adj_gb.txt"))

#------------------------------------------------------------------------------#
# 3) ANÁLISIS DESCRIPTIVO
#------------------------------------------------------------------------------#

# Representamos las tasas crudas por año y guardamos en la carpeta de Figuras Espaciotemporales
S <- length(unique(datos$Code))
T <- length(unique(datos$Year))
T.from <- min(datos$Year)
T.to <- max(datos$Year)
TC <- matrix(datos$TC, S , T , byrow = F)

carto <- carto_gb
for(i in seq(T)){
  carto$var <- TC[, i]
  names(carto)[ncol(carto)] <- T.from + i - 1
}

n.color <- 7
paleta <-  get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)
breaks <- round(seq(min(datos$TC), max(datos$TC), length.out = n.color + 1))
breaks[1] <- floor(min(datos$TC))
breaks[n.color + 1] <- ceiling(max(datos$TC))
breaks <- round(breaks)

if (juntos) {
  
  nombres <- c("TC", as.character(2003:2019))
  names(carto)[4:21] <- nombres
  
  archivo <- paste(paste(paste(paste("Figuras/EspacioTemporal/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), "TC_2002_2019.jpeg", sep = "_")
  graphics.off()
  jpeg(file = archivo, width = 2000, height = 3000)
  print(tm_shape(carto) + 
          tm_polygons(nombres[seq(1, 18, 2)], style = "fixed", n = n.color, breaks = breaks, 
                      border.col = "black", palette = paleta, legend.reverse = T) + 
          tm_layout(scale = 5, legend.outside = T, title = "", panel.show = T, panel.labels = seq(T.from, T.to, 2)) + 
          tm_facets(ncol = 3))
  dev.off()
  
  names(carto)[4] <- "2002"
  
} else {
  
  for (i in 1:18) {
    
    carto_aux <- carto[, c(1:3, i + 3)]
    names(carto_aux)[4] <- "TC"
    
    archivo <- paste(paste(paste(paste(paste(paste("Figuras/EspacioTemporal/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), "TC", sep = "_"), as.character(2001 + i), sep = "_"), ".jpeg", sep = "")
    graphics.off()
    jpeg(file = archivo, width = 1500, height = 2500)
    print(tm_shape(carto_aux) + 
            tm_polygons("TC", style = "fixed", n = n.color, breaks = breaks, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
            tm_layout(main.title = as.character(2001 + i), main.title.size = 5, legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
    dev.off()
    
  }
  
}

# Representamos la evolución de las TC por región y guardamos en la carpeta de
# Figuras Espaciotemporales

aux <- carto[,as.character(T.from:T.to)]
aux$geometry <- NULL

archivo <- paste(paste(paste(paste("Figuras/EspacioTemporal/Evolución_TC", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), ".pdf", sep = "")
graphics.off()
pdf(file = archivo)
matplot(t(aux), type = "l", xlab = "Año", ylab = "Tasa cruda", xaxt = "n")
axis(1, at = seq(1, T, 2), labels = seq(T.from, T.to, 2))
dev.off()

#------------------------------------------------------------------------------#
# 4) MODELOS ESPACIO-TEMPORALES (Leroux + RW1, RW2)
#------------------------------------------------------------------------------#

if (saltar) {
  ## Leroux + RW1 + TypeII ##
  RW1.TypeII <- STCAR_INLA(carto = carto, data = datos, ID.area = "Code", ID.year = "Year", O = "Obs", E = "Population",
                           W = adj_gb, spatial = "Leroux", temporal = "rw1", interaction = "TypeII",
                           model = "global", strategy = estrategia)
} else {
  
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

  Tabla <- data.frame(DIC = unlist(lapply(model_list, function(x) x$dic$dic)),
                      WAIC = unlist(lapply(model_list, function(x) x$waic$waic)),
                      LS = unlist(lapply(model_list, function(x) -sum(log(x$cpo$cpo)))))
  print(Tabla)
  
}

#------------------------------------------------------------------------------#
# 5) REPRESENTACIONES GRÁFICAS
#------------------------------------------------------------------------------#
model <- RW1.TypeII

# Representamos las distribuciones a posteriori de los hiperparámetros y del valor base
# y guardamos en la carpeta de Figuras Espaciotemporales
archivo <- paste(paste(paste(paste("Figuras/EspacioTemporal/Posteriores", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), "ST.pdf", sep = "_")
graphics.off()
pdf(file = archivo, width = 10, height = 7)
par(mfrow = c(2, 2))
plot(inla.smarginal(model$marginals.hyperpar$`Beta for ID.area`), type = "l", lwd = 1.5, main = "Hiperparámetro de Leroux", xlab = expression(lambda), ylab = expression(tilde(p)(paste(lambda, "|", y))))
abline(v = model$summary.hyperpar["Beta for ID.area","0.5quant"], col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Precision for ID.area`), type = "l", lwd = 1.5, main = "Hiperparámetro espacial", xlab = expression(tau[u]), ylab = expression(tilde(p)(paste(tau[u], "|", y))))
abline(v = model$summary.hyperpar["Precision for ID.area","0.5quant"], col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Precision for ID.year`), type = "l", lwd = 1.5, main = "Hiperparámetro temporal", xlab = expression(tau[phi]), ylab = expression(tilde(p)(paste(tau[phi], "|", y))))
abline(v = model$summary.hyperpar["Precision for ID.year","0.5quant"], col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Precision for ID.area.year`), type = "l", lwd = 1.5, main = "Hiperparámetro de interacción", xlab = expression(tau[delta]), ylab = expression(tilde(p)(paste(tau[delta], "|", y))))
abline(v = model$summary.hyperpar["Precision for ID.area.year","0.5quant"], col = "blue")
par(mfrow = c(1, 1))
dev.off()

# Representamos el Patrón Espacial (común para todo el periodo) y guardamos en la carpeta de Figuras Espaciotemporales
aux <- lapply(model$marginals.random$ID.area, function(x) inla.tmarginal(function(y) 100000 * exp(model$summary.fixed$`0.5quant`) * exp(y), x))
pattern.S <- unlist(lapply(aux, function(x) inla.qmarginal(0.5, x)))
carto$TA <- pattern.S

n.color <- 7
paleta <- get_brewer_pal("-RdYlGn", n = n.color, contrast = 0.8, plot = F)
breaks <- 100000 * seq(min(pattern.S), max(pattern.S), length.out = n.color + 1)
breaks[1] <- floor(breaks[1])
breaks[n.color + 1] <- ceiling(breaks[n.color + 1])
breaks <- round(breaks)
titulo <- paste(tipo, sexo, sep = " ")

archivo <- paste(paste(paste(paste("Figuras/EspacioTemporal/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), "Patrón_Espacial.jpeg", sep = "_")
graphics.off()
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("TA", style = "equal", n = n.color, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(main.title = titulo, main.title.size = 5, legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()

# Representamos el Patrón Temporal (común para todas las áreas) y guardamos en la carpeta de Figuras Espaciotemporales
archivo <- paste(paste(paste(paste("Figuras/EspacioTemporal/Grafica", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), "Patrón_Temporal.pdf", sep = "_")

aux <- lapply(model$marginals.random$ID.year, function(x) inla.tmarginal(function(y) 100000 * exp(model$summary.fixed$`0.5quant`) * exp(y), x))
temporal.risk <- unlist(lapply(aux, function(x) inla.emarginal(function(y) y, x)))
q1 <- unlist(lapply(aux, function(x) inla.qmarginal(0.025, x)))
q2 <- unlist(lapply(aux, function(x) inla.qmarginal(0.5, x)))
q3 <- unlist(lapply(aux, function(x) inla.qmarginal(0.975, x)))

graphics.off()
pdf(file = archivo)
x <- 1:T
plot(temporal.risk, type = "l", ylim = range(temporal.risk) * c(1 / 1.1, 1.1),
     xaxt = "n", xlab = "Año", ylab = "TA (casos x 100,000 habitantes)", main = titulo)
axis(1, at = seq(1, T, 2), labels = seq(T.from, T.to, 2))
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1, tail(q3, 1), rev(q3), q1[1])
polygon(X.Vec, Y.Vec, col = "gray", border = NA)
lines(q2)
abline(h = 100000 * exp(model$summary.fixed$`0.5quant`), lty = 2, col = "blue")
dev.off()

# Representamos las tasas ajustadas por año y guardamos en la carpeta de Figuras Espaciotemporales
TA <- matrix(model$summary.fitted.values$`0.5quant` * 100000, S , T , byrow = F)

carto <- carto_gb
for(i in seq(T)){
  carto$var <- TA[, i]
  names(carto)[ncol(carto)] <- T.from + i - 1
}

n.color <- 7
paleta <- get_brewer_pal("-RdYlGn", n = n.color, contrast = 0.8, plot = F)
breaks <- 100000 * seq(min(model$summary.fitted.values$`0.5quant`), max(model$summary.fitted.values$`0.5quant`), length.out = n.color + 1)
breaks[1] <- floor(breaks[1])
breaks[n.color + 1] <- ceiling(breaks[n.color + 1])
breaks <- round(breaks)

if (juntos) {
  
  nombres <- c("TA", as.character(2003:2019))
  names(carto)[4:21] <- nombres
  
  archivo <- paste(paste(paste(paste("Figuras/EspacioTemporal/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), "TA_2002_2019.jpeg", sep = "_")
  graphics.off()
  jpeg(file = archivo, width = 2000, height = 3000)
  print(tm_shape(carto) + 
          tm_polygons(nombres[seq(1, 18, 2)], style = "fixed", n = n.color, breaks = breaks, 
                      border.col = "black", palette = paleta, legend.reverse = T) + 
          tm_layout(scale = 5, legend.outside = T, title = "", panel.show = T, panel.labels = seq(T.from, T.to, 2)) + 
          tm_facets(ncol = 3))
  dev.off()
  
  names(carto)[4] <- "2002"
  
} else {
  
  for (i in 1:18) {
    
    carto_aux <- carto[, c(1:3, i + 3)]
    names(carto_aux)[4] <- "TA"
    
    archivo <- paste(paste(paste(paste(paste(paste("Figuras/EspacioTemporal/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), "TA", sep = "_"), as.character(2001 + i), sep = "_"), ".jpeg", sep = "")
    graphics.off()
    jpeg(file = archivo, width = 1500, height = 2500)
    print(tm_shape(carto_aux) + 
            tm_polygons("TA", style = "fixed", n = n.color, breaks = breaks, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
            tm_layout(main.title = as.character(2001 + i), main.title.size = 5, legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
    dev.off()
    
  }
  
}

# Representamos las probabilidades de exceso y guardamos en la carpeta de Figuras Espaciotemporales
Prob <- unlist(lapply(model$marginals.fitted.values, function(x) 1 - inla.pmarginal(umbral / 100000, x)))

Prob <- matrix(Prob, S, T, byrow = F)

carto <- carto_gb
for(i in seq(T)){
  carto$var <- Prob[,i]
  names(carto)[ncol(carto)] <- T.from + i - 1
}

n.color <- 7
paleta <- get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)

if (juntos) {
  
  nombres <- c("Prob", as.character(2003:2019))
  names(carto)[4:21] <- nombres
  
  archivo <- paste(paste(paste(paste("Figuras/EspacioTemporal/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), "Prob_2002_2019.jpeg", sep = "_")
  graphics.off()
  jpeg(file = archivo, width = 2000, height = 3000)
  print(tm_shape(carto) + 
          tm_polygons(nombres[seq(1, 18, 2)], style = "fixed", n = n.color, interval.closure = "left", 
                      breaks = c(0, 0.05, 0.1, 0.2, 0.8, 0.9, 0.95, 1), labels = c("[0-0.05)", "[0-05-0.1)", "[0.1-0.2)", "[0.2-0.8)", "[0.8-0.9)", "[0.9-0.95)", "[0.95-1]"),
                      border.col = "black", palette = paleta, legend.reverse = T) + 
          tm_layout(scale = 5, title = "", legend.outside = T, panel.show = T, panel.labels = seq(T.from, T.to, 2)) +
          tm_facets(ncol = 3))
  dev.off()
  
  names(carto)[4] <- "2002"
  
} else {
  
  for (i in 1:18) {
    
    carto_aux <- carto[, c(1:3, i + 3)]
    names(carto_aux)[4] <- "Prob"
    
    archivo <- paste(paste(paste(paste(paste(paste("Figuras/EspacioTemporal/Mapa", tipo, sep = "_"), cancer, sep = "_"), sexo, sep = "_"), "Prob", sep = "_"), as.character(2001 + i), sep = "_"), ".jpeg", sep = "")
    graphics.off()
    jpeg(file = archivo, width = 1500, height = 2500)
    print(tm_shape(carto_aux) + 
            tm_polygons("Prob", style = "fixed", n = n.color, interval.closure = "left",
                        breaks = c(0, 0.05, 0.1, 0.2, 0.8, 0.9, 0.95, 1), labels = c("[0-0.05)", "[0-05-0.1)", "[0.1-0.2)", "[0.2-0.8)", "[0.8-0.9)", "[0.9-0.95)", "[0.95-1]"),
                        border.col = "black", palette = paleta, legend.reverse = T) + 
            tm_layout(main.title = as.character(2001 + i), main.title.size = 5, legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
    dev.off()
    
  }
  
}
