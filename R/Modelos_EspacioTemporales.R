rm(list = ls())

library(bigDM)
library(INLA)
library(sf)
library(spdep)
library(tmap)
library(tmaptools)

#------------------------------------------------------------------------------#
# 1) LECTURA DE DATOS
#------------------------------------------------------------------------------#

# Leemos los datos de cáncer de mujeres y hombres
datos_gb_f <- read.csv("../Datos/Datos_gb_f.csv", header=T)
datos_gb_m <- read.csv("../Datos/Datos_gb_m.csv", header=T)

# Leemos la cartografía de Gran Bretaña
carto_gb <- st_read("../Datos/Carto/carto_gb.shp")
head(carto_gb)

# Leemos la matriz de adyacencia (grafo conexo)
adj_gb <- as.matrix(read.table("../Datos/Carto/adj_gb.txt"))


#------------------------------------------------------------------------------------------------------#
# 2) Análisis espacio-temporal: Nuevos casos de cáncer de pulmón femenino/masculino (periodo 2002-2019)
#------------------------------------------------------------------------------------------------------#

# Ordenamos por Región y Año
datos_gb_f <- datos_gb_f[order(datos_gb_f$Year,datos_gb_f$Code), ]
datos_gb_m <- datos_gb_m[order(datos_gb_m$Year,datos_gb_m$Code), ]
carto_gb <- carto_gb[order(carto_gb$Code), ]

# Seleccionamos las variables de interés
inc_lung_f <- datos_gb_f[, c("Code","Year","Population","Count_Lung")]
inc_lung_m <- datos_gb_m[, c("Code","Year","Population","Count_Lung")]

# Cambiamos el nombre de la variable para que cuadre con las fórmulas
names(inc_lung_f)[4] <- "Obs"
names(inc_lung_m)[4] <- "Obs"


#------------------------------------------------------------------------------#
# 3) TASA CRUDA: Casos/Población x 100.000
#------------------------------------------------------------------------------#

# Calculamos la tasa cruda por área y año
S <- nrow(carto_gb)
T <- length(unique(inc_lung_f$Year))
T.from <- min(inc_lung_f$Year)
T.to <- max(inc_lung_f$Year)

inc_lung_f$TC <- inc_lung_f$Obs/inc_lung_f$Population*100000

TC.year <- as.data.frame(matrix(inc_lung_f$TC,S,T,byrow=F))
colnames(TC.year) <- unique(inc_lung_f$Year)
TC.year$Code <- unique(inc_lung_f$Code)

carto <- merge(carto_gb,TC.year)
head(carto)


# Evolución geográfica de las TC por año
archivo <- "Figuras/Mapa_CancerPulmon_f.jpeg"

n.color <- 7
paleta <-  get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)
breaks <- round(seq(min(inc_lung_f$TC),max(inc_lung_f$TC),length.out=n.color+1))
breaks[1] <- floor(min(inc_lung_f$TC))
breaks[n.color+1] <- ceiling(max(inc_lung_f$TC))

graphics.off()
jpeg(file = archivo, width = 2000, height = 1200)
print(tm_shape(carto) + 
        tm_polygons(as.character(seq(T.from,T.to,2)), style = "fixed", n = n.color, breaks=breaks, 
                    lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 2, legend.text.size = 1.5, legend.position = c("left", "top")))
dev.off()

# Evolución temporal de las TC por región
aux <- carto[,as.character(T.from:T.to)]
aux$geometry <- NULL

matplot(t(aux), type="l", xlab="Año", ylab="Tasa cruda", xaxt="n")
axis(1, at=seq(1,T,2), labels=seq(T.from,T.to,2))


#------------------------------------------------------------------------------#
# 4) MODELOS ESPACIO-TEMPORALES (Leroux + RW1)
#------------------------------------------------------------------------------#

## Leroux + RW1 (modelo aditivo) ##
RW1.Additive <- STCAR_INLA(carto=carto, data=inc_lung_f, ID.area="Code", ID.year="Year", O="Obs", E="Population",
                           W=adj_gb, spatial="Leroux", temporal="rw1", interaction="none",
                           model="global", strategy="simplified.laplace")

## Leroux + RW1 + TypeI ##
RW1.TypeI <- STCAR_INLA(carto=carto, data=inc_lung_f, ID.area="Code", ID.year="Year", O="Obs", E="Population",
                        W=adj_gb, spatial="Leroux", temporal="rw1", interaction="TypeI",
                        model="global", strategy="simplified.laplace")

## Leroux + RW1 + TypeII ##
RW1.TypeII <- STCAR_INLA(carto=carto, data=inc_lung_f, ID.area="Code", ID.year="Year", O="Obs", E="Population",
                         W=adj_gb, spatial="Leroux", temporal="rw1", interaction="TypeII",
                         model="global", strategy="simplified.laplace")

## Leroux + RW1 + TypeIII ##
RW1.TypeIII <- STCAR_INLA(carto=carto, data=inc_lung_f, ID.area="Code", ID.year="Year", O="Obs", E="Population",
                          W=adj_gb, spatial="Leroux", temporal="rw1", interaction="TypeIII",
                          model="global", strategy="simplified.laplace")

## Leroux + RW1 + TypeIV ##
RW1.TypeIV <- STCAR_INLA(carto=carto, data=inc_lung_f, ID.area="Code", ID.year="Year", O="Obs", E="Population",
                         W=adj_gb, spatial="Leroux", temporal="rw1", interaction="TypeIV",
                         model="global", strategy="simplified.laplace")

# Comparamos los modelos
MODELOS <- list(Aditivo=RW1.Additive,TypeI=RW1.TypeI,TypeII=RW1.TypeII,TypeIII=RW1.TypeIII,TypeIV=RW1.TypeIV)

Tabla <- data.frame(mean.deviance=unlist(lapply(MODELOS, function(x) x$dic$mean.deviance)),
                    p.eff=unlist(lapply(MODELOS, function(x) x$dic$p.eff)),
                    DIC=unlist(lapply(MODELOS, function(x) x$dic$dic)),
                    WAIC=unlist(lapply(MODELOS, function(x) x$waic$waic)),
                    LS=unlist(lapply(MODELOS, function(x) -sum(log(x$cpo$cpo)))))
print(Tabla)


#------------------------------------------------------------------------------#
# 5) REPRESENTACIONES GRÁFICAS
#------------------------------------------------------------------------------#
model <- RW1.TypeII

# Distribuciones a posteriori de los hiperparámetros
archivo <- "Figuras/Posteriores_ST.pdf"

model$summary.fixed
model$summary.hyperpar

graphics.off()
pdf(file = archivo, width = 10, height = 5)
par(mfrow = c(2, 3))
plot(inla.smarginal(model$marginals.fixed$`(Intercept)`), type = "l", lwd = 1.5, main = "Valor base", xlab = expression(eta), ylab = expression(tilde(p)(paste(eta, "|", y))))
abline(v = model$summary.fixed$`0.5quant`, col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Precision for ID.area`), type = "l", lwd = 1.5, main = "Hiperparámetro", xlab = expression(tau[s]), ylab = expression(tilde(p)(paste(tau[s], "|", y))))
abline(v = model$summary.hyperpar["Precision for ID.area","0.5quant"], col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Beta for ID.area`), type = "l", lwd = 1.5, main = "Hiperparámetro", xlab = expression(lambda[s]), ylab = expression(tilde(p)(paste(lambda[s], "|", y))))
abline(v = model$summary.hyperpar["Beta for ID.area","0.5quant"], col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Precision for ID.year`), type = "l", lwd = 1.5, main = "Hiperparámetro", xlab = expression(tau[t]), ylab = expression(tilde(p)(paste(tau[t], "|", y))))
abline(v = model$summary.hyperpar["Precision for ID.year","0.5quant"], col = "blue")
plot(inla.smarginal(model$marginals.hyperpar$`Precision for ID.area.year`), type = "l", lwd = 1.5, main = "Hiperparámetro", xlab = expression(tau[st]), ylab = expression(tilde(p)(paste(tau[st], "|", y))))
abline(v = model$summary.hyperpar["Precision for ID.area.year","0.5quant"], col = "blue")
par(mfrow = c(1, 1))
dev.off()


# Patrón espacial (común para todo el periodo)
pattern.S <- unlist(lapply(model$marginals.random$ID.area, function(x) inla.emarginal(function(y) 100000*tasa.global*exp(y),x)))
carto$pattern.S <- pattern.S


archivo <- "Figuras/Mapa_CancerPulmon_f_SpatialPattern.jpeg"
n.color <- 7
paleta <- get_brewer_pal("-RdYlGn", n = n.color, contrast = 0.8, plot = F)

graphics.off()
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("pattern.S", style = "equal", n=n.color, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()


# Patrón temporal (común para todas las áreas)
archivo <- "Figuras/Mapa_CancerPulmon_f_TemporalPattern.jpeg"

aux <- lapply(model$marginals.random$ID.year, function(x) inla.tmarginal(function(y) 100000*tasa.global*exp(y),x))
temporal.risk <- unlist(lapply(aux, function(x) inla.emarginal(function(y) y,x)))
q1 <- unlist(lapply(aux, function(x) inla.qmarginal(0.025,x)))
q2 <- unlist(lapply(aux, function(x) inla.qmarginal(0.5,x)))
q3 <- unlist(lapply(aux, function(x) inla.qmarginal(0.975,x)))

graphics.off()
jpeg(file = archivo)
x <- 1:T
plot(temporal.risk, type="l", ylim=range(temporal.risk)*c(1/1.1,1.1),
     xaxt="n", xlab="Año", ylab="Casos x 100,000 habitantes", main="Patrón temporal")
axis(1, at=seq(1,T,2), labels=seq(T.from,T.to,2))
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1, tail(q3, 1), rev(q3), q1[1])
polygon(X.Vec, Y.Vec, col="gray", border = NA)
lines(q2)
abline(h=100000*tasa.global, lty=2, col="blue")
dev.off()


# Mapas con las tasas ajustadas por año
TA <- matrix(model$summary.fitted.values$`0.5quant`*100000,S,T,byrow=F)

carto <- carto_gb
for(i in seq(T)){
  carto$var <- TA[,i]
  names(carto)[ncol(carto)] <- paste("Year",T.from+i-1,sep=".")
}

archivo <- "Figuras/Mapa_CancerPulmon_f_TA.jpeg"

n.color <- 6
paleta <- get_brewer_pal("-RdYlGn", n = n.color, contrast = 0.8, plot = F)
breaks <- 100000*seq(min(model$summary.fitted.values$`0.5quant`),max(model$summary.fitted.values$`0.5quant`),length.out=n.color+1)
breaks[1] <- floor(breaks[1])
breaks[n.color+1] <- floor(breaks[n.color+1])
breaks <- round(breaks)

graphics.off()
jpeg(file = archivo, width = 2500, height = 1500)
print(tm_shape(carto) + 
        tm_polygons(paste("Year",seq(T.from,T.to,2),sep="."), style = "fixed", breaks=breaks, 
                    lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")) + 
        tm_facets(nrow=2, ncol=5))
dev.off()


# Mapas con las probabilidades a posteriori
tasa.global <- exp(model$summary.fixed$`0.5quant`)
Prob <- unlist(lapply(model$marginals.fitted.values, function(x) 1-inla.pmarginal(tasa.global,x)))
Prob <- matrix(Prob,S,T,byrow=F)
  
carto <- carto_gb
for(i in seq(T)){
  carto$var <- Prob[,i]
  names(carto)[ncol(carto)] <- paste("Year",T.from+i-1,sep=".")
}

archivo <- "Figuras/Mapa_CancerPulmon_f_Prob.jpeg"

n.color <- 5
paleta <- get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)

graphics.off()
jpeg(file = archivo, width = 2500, height = 1500)
print(tm_shape(carto) + 
        tm_polygons(paste("Year",seq(T.from,T.to,2),sep="."), style = "fixed", interval.closure="left",
                    breaks=c(0,0.1,0.2,0.8,0.9,1), labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]"),
                    lwd = 2, border.col = "black", palette = paleta) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")) + 
        tm_facets(nrow=2, ncol=5))
dev.off()
