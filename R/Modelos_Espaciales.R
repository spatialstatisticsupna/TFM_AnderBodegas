rm(list = ls())

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

# Representamos el grafo de vecindad
par(mfrow=c(1,2), pty="s")

carto.nb1 <- poly2nb(carto_gb)
n.comp.nb(carto.nb1)$nc==1

plot(carto_gb$geometry, main="Grafo de vecindad original")
plot(carto.nb1, st_centroid(st_geometry(carto_gb), of_largest_polygon=TRUE),
     pch=19, cex=0.5, col="red", add=TRUE)

carto.nb2 <- mat2listw(adj_gb)$neighbours
n.comp.nb(carto.nb2)$nc==1 

plot(carto_gb$geometry, main="Grafo de vecindad modificado (conexo)")
plot(carto.nb2, st_centroid(st_geometry(carto_gb), of_largest_polygon=TRUE),
     pch=19, cex=0.5, col="red", add=TRUE)


#-------------------------------------------------------------------------------------#
# 2) Análisis espacial: Nuevos casos de cáncer de pulmón femenino/masculino (año 2019)
#-------------------------------------------------------------------------------------#

# Ordenamos por Región y Año
datos_gb_f <- datos_gb_f[order(datos_gb_f$Year,datos_gb_f$Code), ]
datos_gb_m <- datos_gb_m[order(datos_gb_m$Year,datos_gb_m$Code), ]
carto_gb <- carto_gb[order(carto_gb$Code), ]

# Seleccionamos las variables de interés
inc_lung_2019_f <- datos_gb_f[which(datos_gb_f$Year==2019), c("Code","Population","Count_Lung")]
inc_lung_2019_m <- datos_gb_m[which(datos_gb_m$Year==2019), c("Code","Population","Count_Lung")]

# Cambiamos el nombre de la variable para que cuadre con las fórmulas
names(inc_lung_2019_f)[3] <- "Obs"
names(inc_lung_2019_m)[3] <- "Obs"


#------------------------------------------------------------------------------#
# 3) TASA CRUDA: Casos/Población x 100.000
#------------------------------------------------------------------------------#
inc_lung_2019_f$TC <- inc_lung_2019_f$Obs/inc_lung_2019_f$Population*100000
carto <- merge(carto_gb,inc_lung_2019_f)

# Creamos los tres mapas y los guardamos en la carpeta de figuras
archivo <- "Figuras/Mapa_CancerPulmon_2019f.jpeg"

n.color <- 7
paleta <-  get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)
breaks.obs <- round(seq(min(carto$Obs),max(carto$Obs),length.out=n.color+1))
breaks.pop <- round(seq(min(carto$Population),max(carto$Population),length.out=n.color+1))
breaks.TC <- round(seq(min(carto$TC),max(carto$TC),length.out=n.color+1))
breaks.TC[1] <- floor(min(carto$TC))
breaks.TC[n.color+1] <- ceiling(max(carto$TC))

graphics.off()
jpeg(file = archivo, width = 2000, height = 1200)
print(tm_shape(carto) + 
        tm_polygons(c("Obs","Population","TC"), style = "fixed", n = n.color,
                    breaks=list(breaks.obs,breaks.pop,breaks.TC),
                    lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 2, legend.text.size = 1.5, legend.position = c("left", "top")))
dev.off()


#------------------------------------------------------------------------------#
# 4) MODELOS ESPACIALES 
#------------------------------------------------------------------------------#

## Distribuciones no-informativas para los hiperparámetros ##
sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"

lunif = "expression:
          a = 1;
          b = 1;
          beta = exp(theta)/(1+exp(theta));
          logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
          log_jacobian = log(beta*(1-beta));
          return(logdens+log_jacobian)"

# Ajustamos el modelo iCAR
codes <- unique(inc_lung_2019_f$Code)
form1 <- Obs ~ f(Code, model='besag', values=codes, graph=adj_gb, constr=TRUE, hyper=list(prec=list(prior=sdunif)))

model1 <- inla(formula=form1, family="poisson", data=inc_lung_2019_f, E=Population,
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
               control.inla=list(strategy="laplace"))
summary(model1)


# Representamos las distribuciones marginales a posteriori del valor base, del parámetro espacial correspondiente
# a la primera región (Barnsley) y del hiperparámetro. Lo guardamos en la carpeta de figuras
archivo <- "Figuras/Posteriores_SP.pdf"

graphics.off()
pdf(file = archivo, width = 10, height = 5)
par(mfrow = c(1, 3))
plot(inla.smarginal(model1$marginals.fixed$`(Intercept)`), type = "l", lwd = 1.5, main = "Valor base", xlab = expression(eta), ylab = expression(tilde(p)(paste(eta, "|", y))))
abline(v = model1$summary.fixed$`0.5quant`, col = "blue")
plot(inla.smarginal(model1$marginals.random$Code$index.1), type = "l", lwd = 1.5, main = "Parámetro espacial", xlab = expression(u[1]), ylab = expression(tilde(p)(paste(u[1], "|", y))))
abline(v = model1$summary.random$Code$`0.5quant`[1], col = "blue")
plot(inla.smarginal(model1$marginals.hyperpar$`Precision for Code`), type = "l", lwd = 1.5, main = "Hiperparámetro", xlab = expression(tau), ylab = expression(tilde(p)(paste(tau, "|", y))))
abline(v = model1$summary.hyperpar$`0.5quant`, col = "blue")
par(mfrow = c(1, 1))
dev.off()

# Valor base
model1$summary.fixed$`0.5quant`
inla.emarginal(function(x) 100000*exp(x), model1$marginals.fixed$`(Intercept)`)

# Parámetro espacial de Barnsley
model1$summary.random$Code$`0.5quant`[1]
inla.emarginal(function(x) exp(x), model1$marginals.random$Code[[1]])

# Barnsley (manualmente y con fitted values)
100000 * exp(model1$summary.fixed$`0.5quant`) * exp(model1$summary.random$Code$`0.5quant`[1])
100000 * model1$summary.fitted.values$`0.5quant`[1]

# Distribuciones a posteriori de los hiperparámetros
model1$summary.hyperpar


# Ajustamos todos los modelos espaciales

# i.i.d.
form0 <- Obs ~ f(Code, model="iid", values=codes, graph=adj_gb)
model0 <- inla(formula=form0, family="poisson", data=inc_lung_2019_f, E=Population,
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
               control.inla=list(strategy="laplace"))

# iCAR
form1 <- Obs ~ f(Code, model='besag', values=codes, graph=adj_gb, constr=TRUE, hyper=list(prec=list(prior=sdunif)))
model1 <- inla(formula=form1, family="poisson", data=inc_lung_2019_f, E=Population,
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
               control.inla=list(strategy="laplace"))

# BYM
form2 <- Obs ~ f(Code, model="bym", values=codes, graph=adj_gb, constr=TRUE, hyper=list(theta1=list(prior=sdunif), theta2=list(prior=sdunif)))
model2 <- inla(formula=form2, family="poisson", data=inc_lung_2019_f, E=Population,
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
               control.inla=list(strategy="laplace"))

# Leroux
Q <- diag(abs(colSums(adj_gb))) - adj_gb
C <- diag(rep(1, nrow(Q))) - Q

form3 <- Obs ~ f(Code, model="generic1", values=codes, Cmatrix=C, constr=TRUE, hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))
model3 <- inla(formula=form3, family="poisson", data=inc_lung_2019_f, E=Population,
               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
               control.inla=list(strategy="laplace"))


# # NOTA: Los modelos pueden ajustarse directamente utilizando la librería "bigDM"
# library(bigDM)
# model1 <- CAR_INLA(carto=carto, ID.area="Code", O="Obs", E="Population", W=adj_gb,
#                    prior="intrinsic", model="global", strategy="laplace")
# 
# model2 <- CAR_INLA(carto=carto, ID.area="Code", O="Obs", E="Population", W=adj_gb,
#                    prior="BYM", model="global", strategy="laplace")
# 
# model3 <- CAR_INLA(carto=carto, ID.area="Code", O="Obs", E="Population", W=adj_gb,
#                    prior="Leroux", model="global", strategy="laplace")


# Comparamos los modelos
MODELOS <- list(iid=model0,iCAR=model1,BYM=model2,Leroux=model3)

Tabla <- data.frame(mean.deviance=unlist(lapply(MODELOS, function(x) x$dic$mean.deviance)),
                    p.eff=unlist(lapply(MODELOS, function(x) x$dic$p.eff)),
                    DIC=unlist(lapply(MODELOS, function(x) x$dic$dic)),
                    WAIC=unlist(lapply(MODELOS, function(x) x$waic$waic)),
                    LS=unlist(lapply(MODELOS, function(x) -sum(log(x$cpo$cpo)))))
print(Tabla)


#------------------------------------------------------------------------------#
# 5) REPRESENTACIONES GRÁFICAS
#------------------------------------------------------------------------------#
model <- model3

# Representamos en un mapa las tasas ajustadas y la probabilidad de que supere a la tasa media de GB
tasa.global <- exp(model$summary.fixed$`0.5quant`)
  
carto$TA <- model$summary.fitted.values$`0.5quant`*100000
carto$Prob <- unlist(lapply(model$marginals.fitted.values, function(x) 1-inla.pmarginal(tasa.global,x)))


# Creamos los mapas y los guardamos en la carpeta figuras
archivo <- "Figuras/Mapa_CancerPulmon_2019f_TA.jpeg"

n.color <- 7
paleta <- get_brewer_pal("-RdYlGn", n = n.color, contrast = 0.8, plot = F)
values <- seq(min(model$summary.fitted.values$mean),max(model$summary.fitted.values$mean),length.out=n.color+1)*100000
values[1] <- floor(values[1])
values[n.color+1] <- ceiling(values[n.color+1])
values <- round(values)
  
graphics.off()
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("TA", style = "fixed", breaks=values, lwd = 2, border.col = "black", palette = paleta, legend.reverse = T) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()


archivo <- "Figuras/Mapa_CancerPulmon_2019f_Prob.jpeg"
n.color <- 5
paleta <- get_brewer_pal("Blues", n = n.color, contrast = 0.8, plot = F)
graphics.off()
jpeg(file = archivo, width = 1500, height = 2500)
print(tm_shape(carto) + 
        tm_polygons("Prob", style = "fixed", interval.closure="left",
                    breaks=c(0,0.1,0.2,0.8,0.9,1), labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]"),
                    lwd = 2, border.col = "black", palette = paleta) + 
        tm_layout(title = "", legend.title.size = 4, legend.text.size = 3.5, legend.position = c("left", "top")))
dev.off()
