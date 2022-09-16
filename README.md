# TFM_AnderBodegas
Este repositorio contiene el código R para reproducir y replicar el análisis de datos del Trabajo de Fin de Máster ["Modelos espaciales y espaciotemporales en disease mapping"](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/TFM_Ander_Bodegas.pdf) realizado por Ander Bodegas Díez (bajo la supervisión de Aritz Adin y Jaione Etxeberria) en el [Máster de Modelización e Investigación Matemática, Estadística y Computación](https://www.unavarra.es/sites/masteres/ciencias/modelizacion-invest-matematica/presentacion.html) de la Universidad Pública de Navarra.


## Índice

- [Datos](#Datos)
- [Código R](#Código-r)
- [Agradecimientos](#Agradecimientos)


# Datos

Esta carpeta contiene los ficheros con los que se ha realizado el tercer capítulo del trabajo, en el que se ilustra el funcionamiento de los modelos con los datos reales de cáncer.

- [**Datos_gb_f.csv**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/Datos/Datos_gb_f.csv)

  Esta base de datos contiene para cada región, año y tipo de cáncer el número de casos observados y el número de muertes en mujeres correspondiente.
  
- [**Datos_gb_m.csv**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/Datos/Datos_gb_m.csv)

  Esta base de datos contiene para cada región, año y tipo de cáncer el número de casos observados y el número de muertes en hombres correspondiente.
  
- [**Carto**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/Datos/Carto/)

  Esta carpeta contiene los cuatro archivos de la cartografía de Gran Bretaña, además de la matriz de adyacencia conexa.
  
- [**adj_bg.txt**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/Datos/Carto/adj_gb.txt)

  Este archivo es la matriz de adyacencia de las regiones bajo estudio, en la que se entiende que dos regiones son vecinas si comparten frontera. Además, se han añadido 9 conexiones extra para conectar diversas islas y hacer el grafo conexo. 



# Código R
El código de R correspondiente al tercer capítulo del trabajo, en el que se ilustra el funcionamiento de los modelos con los datos reales de cáncer.

- [**Modelos_Espaciales.R**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/R/Modelos_Espaciales.R)

  Este script realiza el análisis espacial (datos de incidencia por cáncer de pulmón en mujeres en el año 2019) ajustando los modelos en INLA, permitiendo reproducir todos los mapas y gráficos presentados en el trabajo. Para analizar otras causas disponibles en la base de datos original, es posible elegir el tipo de cáncer, el sexo, el tipo de dato (incidencia o mortalidad) y el umbral para las probabilidades de exceso, además de otras opciones al comienzo del script.

- [**Modelos_EspacioTemporales.R**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/R/Modelos_EspacioTemporales.R)

  Este script realiza el análisis espaciotemporal (datos de incidencia por cáncer de pulmón en mujeres durante el periodo 2002-2019) ajustando los modelos en INLA, permitiendo reproducir todos los mapas y gráficos presentados en el trabajo. Al igual que en el caso espacial, es posible seleccionar datos correspondientes a otras causas de la base de datos orginal para replicar el estudio realizado.
  
- [**Figuras**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/R/Figuras)

  En esta carpeta se guardan los diferentes mapas y gráficos que crean los dos scripts mencionados. En particular, en la carpeta [**Espacial**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/R/Figuras/Espacial) se guardarán los resultados de [**Modelos_Espaciales.R**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/R/Modelos_Espaciales.R) y en [**EspacioTemporal**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/R/Figuras/EspacioTemporal) se guardarán los de [**Modelos_EspacioTemporales.R**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/R/Modelos_EspacioTemporales.R).

  
# Agradecimientos
Este Trabajo Fin de Máster ha sido realizado bajo la financiación de las Ayudas de Iniciación a la Investigación de la Universidad Pública de Navarra en el ámbito de sus Institutos de Investigación durante el curso académico 2021/2022 ([resolución nº 602/2022](https://www2.unavarra.es/gesadj/centroJeronimoAyanz/JDA22Res.%20602%20TFM_INSTITUTOS.pdf)).
