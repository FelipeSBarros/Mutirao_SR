rm(list=ls())
# bibliotecas
library(raster)
library(rgdal)
library(RColorBrewer)
library(rgeos)
library(data.table)
library(ggplot2)
library(data.table)

source('./LULCC_fct.R')
args(ForIndex)
args(ForIndexZone)

# ForIndex ----
# Vários anos
satFolder <- list('./ModData/1985/', './ModData/1991/', './ModData/1995/', './ModData/2000/', './ModData/2005/', './ModData/2010/')

args(ForIndex)
teste <- lapply(satFolder, FUN=ForIndex, studyArea = readOGR(dsn = './shp/', layer='Mutirao_UTM_SIR'), idField = 'iis_id' , group='OBRA', idRef = NULL, stats='mean', rastLayer=c(8,9,10), respRatio = FALSE)
 
# teste <- rbindlist(teste)
teste2 <- na.omit(rbindlist(teste))
?grep

teste2$data2 <- gsub('./ModData/',"",teste2$data)
teste2$data2 <- gsub('/','',teste2$data2)

df <- teste2[which(teste2$OBRA==unique(teste2$OBRA)[6:8]),]
unique(df$OBRA)

library(ggplot2)
ggplot(df, aes(x=as.factor(data2), y=mean8, group=as.factor(iis_id), linetype=as.factor(iis_id))) + geom_line(size=1) + geom_point() + theme(text=element_text(size=16)) + facet_grid(. ~ OBRA)

#ggsave('./teste.png',dpi=300)

