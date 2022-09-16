# TFM_AnderBodegas
Este repositorio contiene el código R para reproducir y replicar el análisis de datos del Trabajo de Fin de Máster ["Modelos espaciales y espaciotemporales en disease mapping"](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/TFM_Ander_Bodegas.pdf) realizado por Ander Bodegas Díez (bajo la supervisión de Aritz Adin y Jaione Etxeberria) en el [Máster de Modelización e Investigación Matemática, Estadística y Computación](https://www.unavarra.es/sites/masteres/ciencias/modelizacion-invest-matematica/presentacion.html) de la Universidad Pública de Navarra.


## Índice

- [Datos](#Datos)
- [Código R](#Código-r)
- [Agradecimientos](#Agradecimientos)


# Datos

Esta carpeta contiene los ficheros con los que se ha realizado el tercer capítulo del trabajo, en el que se ilustra el funcionamiento de los modelos con los datos reales de cáncer.

- [**Datos_gb_f.csv**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/Datos/Datos_gb_f.csv) y [**Datos_gb_m.csv**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/Datos/Datos_gb_m.csv)

  Bases de datos de incidencia y mortalidad po cáncer en la isla de Gran Bretaña para 11 causas (leucemia, mama, cervix, melanoma, hígado, colorectal, pancreas, estómago, pulmón, vegija y esófago) desagregadas por área y año. Cada fichero contiene las siguientes variables:
  
    - **_Code_**: código identificador del área (S=142 regiones)
    - **_Year_**: identificador del año (periodo 2002-2019)
    - **_Population_**: población en riesgo
    - **_Count_All_**: número total de casos registrados
    - **_Deaths_All_**: número total de muertes registradas
    - **_Count_xxx_**: número total de casos registrados para cada una de las 11 causas 
    - **_Deaths_xxx_**: número total de muertes registradas para cada una de las 11 causas

  _*Fuente de datos_: Sistema nacional de salud para Inglaterra [(NHS England)](https://www.cancerdata.nhs.uk/incidence_and_mortality), Gales [(NHS Wales)](https://phw.nhs.wales/services-and-teams/welsh-cancer-intelligence-and-surveillance-unit-wcisu/) y Escocia [(NHS Scotland)](https://www.opendata.nhs.scot/dataset).
  
- [**Carto**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/Datos/Carto/)

  Esta carpeta contiene la cartografía (archivos _shapefile_) de las regiones de Gran Bretaña (106 _clinical commissioning group0s_ para Inglaterra, 22 _local authorities_ para Gales y 14 _health boards_ para Escocia), además de la matriz de adyacencia espacial.
  
- [**adj_bg.txt**](https://github.com/spatialstatisticsupna/TFM_AnderBodegas/blob/main/Datos/Carto/adj_gb.txt)

  Este archivo contiene la matriz de adyacencia (binaria) correspondiente al grafo de vecindad de las 142 regiones bajo estudio. El grafo original ha sido modificado (se han añadido 9 conexiones extra) para conectar las islas y así obtener un grafo conexo.



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
