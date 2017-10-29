rm(list=ls())
gc()

# TESTE Com projetos aleatorios com Indices de vegetacao
# bibliotecas
library(raster)
library(rgdal)
library(RColorBrewer)
library(rgeos)
library(data.table)
library(ggplot2)

## -----
source('./LULCC_fct.R')
args(ForIndex) #This function is based only on spectral distance or mean. It is not needed a R output from Segmentation_fct

args(ForIndexZone)  #This function is based on the unsupervised results. Must have a R output from Segmentation_fct
## ----

# ForIndex ----
# Vários anos
satFolder <- list('./ModData/1985/', './ModData/1991/', './ModData/1995/', './ModData/2000/', './ModData/2005/', './ModData/2010/')

# Loading project polygons
proj_pol <- readOGR(dsn = './shp/', layer='Mutirao_UTM_SIR')
proj_pol@data$AREA <- gArea(proj_pol, byid = TRUE)
head(proj_pol@data)

#Subsetting to a few restoration areas
# Morro do encontro 126, 127
# Rocinha 149 - 158
# Vale dos eucaliptos 243 - 248
#
#ids <- c(390:394, 126,127, 149:153, 156:158,243:245, 247,238:242,248) #read.csv('./ids.csv')
#class(ids)
#dim(ids)
Rocinha <- proj_pol[which(proj_pol@data$iis_id %in% c(149:153,156:158)), ]
unique(Rocinha@data$OBRA)
plot(Rocinha)
#text(Rocinha,(Rocinha@plotOrder))

# Running function
random_result <- lapply(satFolder, FUN=ForIndex, studyArea = Rocinha, idField = 'iis_id', group = "OBRA", idRef = NULL, stats='mean', rastLayer=c(8, 9, 10), respRatio = FALSE)


# removendo NA do resultado
random_result <- na.omit(rbindlist(random_result))
random_result$data <- gsub('./ModData/',"",random_result$data)
random_result$data <- gsub('/','',random_result$data)
head(random_result)

ggplot(random_result, aes(x=as.factor(data), y=mean8, group=as.factor(iis_id), linetype=as.factor(iis_id), colour=as.factor(iis_id))) + geom_line(size=1) + geom_point() + theme(text=element_text(size=16)) + facet_grid(. ~ OBRA) + facet_wrap(~ OBRA) 
#ggsave('./RocinhaMeanNDVI_Result.png')
# + guides(fill=guide_legend(title="Restoration area")) + scale_fill_hue(guide = guide_legend(title = NULL))

ggplot(random_result, aes(x=as.factor(data), y=mean9, group=as.factor(iis_id), linetype=as.factor(iis_id), colour=as.factor(iis_id))) + geom_line(size=1) + geom_point() + theme(text=element_text(size=16)) + facet_grid(. ~ OBRA) + facet_wrap(~ OBRA)
#ggsave('./RocinhaMeanEVI_Result.png')

ggplot(random_result, aes(x=as.factor(data), y=mean10, group=as.factor(iis_id), linetype=as.factor(iis_id), colour=as.factor(iis_id))) + geom_line(size=1) + geom_point() + theme(text=element_text(size=16)) + facet_grid(. ~ OBRA) + facet_wrap(~ OBRA)
#ggsave('./RocinhaMeanSAVI_Result.png')

# Tests with standart deviation ----
ggplot(random_result, aes(x=as.factor(data), y=mean10, group=as.factor(id_fct), linetype=as.factor(id_fct))) + geom_line(size=1) + geom_point() + theme(text=element_text(size=16)) + facet_grid(. ~ OBRA) + facet_wrap(~ OBRA) +geom_errorbar(aes(ymin=sd10-mean10, ymax=sd10+mean10), width=.2)

ggplot(random_result, aes(x=as.factor(data), y=mean8, group=as.factor(id_fct), linetype=as.factor(id_fct))) + geom_line(size=1) + geom_point() + theme(text=element_text(size=16)) + facet_grid(. ~ OBRA) + facet_wrap(~ OBRA) + geom_ribbon(aes(ymin=sd8-mean8, ymax=sd8+mean8), alpha=0.2)

# estimating temp variation ----
# Running function
temp_result <- lapply(satFolder, FUN=ForIndex, studyArea = Rocinha, idField = 'iis_id', group = "OBRA", idRef = NULL, stats='mean', rastLayer=c(6), respRatio = FALSE)

# removendo NA do resultado
temp_result <- na.omit(rbindlist(temp_result))
temp_result$data <- gsub('./ModData/',"",temp_result$data)
temp_result$data <- gsub('/','',temp_result$data)
head(temp_result)

ggplot(temp_result, aes(x=as.factor(data), y=mean6, group=as.factor(iis_id), linetype=as.factor(iis_id), colour=as.factor(iis_id))) + geom_line(size=1) + geom_point() + theme(text=element_text(size=16)) + facet_grid(. ~ OBRA) + facet_wrap(~ OBRA) 
#ggsave('./RocinhaMeanNDVI_Result.png')

# Using response ratio ----
RocinhaRR <- proj_pol[which(proj_pol@data$iis_id %in% c(149:153,156:158, 666:668)), ]
unique(RocinhaRR@data$OBRA)
plot(RocinhaRR)

RR_result <- lapply(satFolder, FUN=ForIndex, studyArea = RocinhaRR, idField = 'iis_id', group = "OBRA", idRef = 666, stats='mean', rastLayer=c(8), respRatio = T)

RR_result <- na.omit(rbindlist(RR_result))
RR_result$data <- gsub('./ModData/',"",RR_result$data)
RR_result$data <- gsub('/','',RR_result$data)
head(RR_result)

ggplot(RR_result, aes(x=as.factor(data), y=RR8, group=as.factor(iis_id), linetype=as.factor(iis_id), colour=as.factor(iis_id))) + geom_line(size=1) + geom_point() + theme(text=element_text(size=16)) + facet_grid(. ~ OBRA) + facet_wrap(~ OBRA) 
#ggsave('./RocinhaRespRatio_Result.png')

# Using distance ----
RocinhaRR <- proj_pol[which(proj_pol@data$iis_id %in% c(149:153,156:158, 666:668)), ]
unique(RocinhaRR@data$OBRA)
plot(RocinhaRR)

RR_result <- lapply(satFolder, FUN=ForIndex, studyArea = RocinhaRR, idField = 'iis_id', group = "OBRA", idRef = 666, stats='mean', rastLayer=c(8), respRatio = F)

RR_result <- na.omit(rbindlist(RR_result))
RR_result$data <- gsub('./ModData/',"",RR_result$data)
RR_result$data <- gsub('/','',RR_result$data)
head(RR_result)

ggplot(RR_result, aes(x=as.factor(data), y=Dist8, group=as.factor(iis_id), linetype=as.factor(iis_id), colour=as.factor(iis_id))) + geom_line(size=1) + geom_point() + theme(text=element_text(size=16)) + facet_grid(. ~ OBRA) + facet_wrap(~ OBRA) 
ggsave('./RocinhaDist_Result.png')

# ForIndexZone ----
satFolder <- list('./ModData/1985/', './ModData/1991/', './ModData/1995/', './ModData/2000/', './ModData/2005/', './ModData/2010/')

DoisIrmaos <- lapply(satFolder, FUN=ForIndexZone, studyArea = readOGR(dsn = './shp/', layer='DoisIrmaos'), idField = 'SETOR', idRef = NULL, stats='max')

#removendo NA do resultado
DoisIrmaos <- na.omit(rbindlist(DoisIrmaos))
DoisIrmaos <- DoisIrmaos[-which(DoisIrmaos$AREA<3600),] #4 pixels
DoisIrmaos$ano <- c(rep(1985, nrow(DoisIrmaos)/length(satFolder)), rep(1991, nrow(DoisIrmaos)/length(satFolder)), rep(1995, nrow(DoisIrmaos)/length(satFolder)), rep(2000, nrow(DoisIrmaos)/length(satFolder)), rep(2005, nrow(DoisIrmaos)/length(satFolder)), rep(2010, nrow(DoisIrmaos)/length(satFolder)))

ggplot(DoisIrmaos, aes(x=as.factor(ano), y=dist, group=as.factor(SETOR), colour=as.factor(SETOR), shape=as.factor(SETOR), linetype = as.factor(SETOR))) + geom_line(size=1) + geom_point() + scale_colour_brewer(palette="Set1") + theme(text=element_text(size=16))
ggsave('./DoisIrmaos.png',dpi=300)

#----
#library(scatterplot3d)
#s3d <- scatterplot3d(fuzzy$membership[,c(1:3)], color=fuzzy$cluster, type="h", angle=55, scale.y=0.7, #pch=16)
#plot(s3d)

#Teste CVA----
#library(RStoolbox)
#?rasterCVA
teste85 <- stack("./ModData/1985/stack_1985.tif")
teste2010 <- stack("./ModData/2010/stack_2010.tif")
cva85_10 <- rasterCVA(teste85[[3:4]], teste2010[[3:4]], tmf = 2)
writeRaster(cva85_10, './CVA.tif')
#par(mfrow=c(2,2))
#plotRGB(teste85, 3,2,1, stretch='lin')
#plotRGB(teste2010, 3,2,1, stretch='lin')
#plot(cva85_10[[1]])
#plot(cva85_10[[2]])
#dev.off()

#cva <- cva85_10[[1]]
#plot(cva==0)
#cva[cva==0] <- NA
#plot(cva)

#LandChange ----
#Biomass loss
#biomassloss <- cva<=90
#biomassloss[biomassloss==0]<-NA
#plot(biomassloss)

#plotRGB(teste2010, 4,2,1, stretch='lin')
#plot(biomassloss, add=TRUE, legend=FALSE)

# Forest cleaning
#forestclearing <- ((cva<=180)-biomassloss)
#forestclearing[forestclearing==0] <- NA
#plot(forestclearing)

#plotRGB(teste2010, 4,2,1, stretch='lin')
#plot(forestclearing, add=TRUE, legend=FALSE)

# biomass gain
#biomassgain <- ((cva<=270)-(cva<=180))
#biomassgain[biomassgain==0]<-NA
#plot(biomassgain)

#plotRGB(teste2010, 4,2,1, stretch='lin')
#plot(biomassgain, add=TRUE, legend=FALSE)

# Regrowth
#regrowth <-  ((cva<=360)-(cva<=270))
#regrowth[regrowth==0]<-NA
#plot(regrowth)

#plotRGB(teste2010, 4,2,1, stretch='lin')
#plot(regrowth, add=TRUE, legend=FALSE)



library(RColorBrewer)
pallete <- brewer.pal(8, 'Accent')

par(mfrow=c(2,2))
plotRGB(teste85, 4, 3, 2, stretch='hist')
plot(fuzzy85, col=pallete, axes=FALSE, box=FALSE)
plot(fuzzy91, col=pallete, axes=FALSE, box=FALSE, legend=FALSE)
plot(fuzzy95, col=pallete, axes=FALSE, box=FALSE, legend=FALSE)
dev.off()

par(mfrow=c(2,2))
plotRGB(teste2010, 4, 3, 2, stretch='hist')
plot(fuzzy10, col=pallete, axes=FALSE, box=FALSE, legend=FALSE)
plot(fuzzy10, col=pallete, axes=FALSE, box=FALSE)
dev.off()
#----
# Testes 
library(scatterplot3d)
s3d <- scatterplot3d(fuzzy$membership[,c(1:3)], color=fuzzy$cluster, type="h", angle=55, scale.y=0.7, pch=16)
plot(s3d)

#Teste CVA----
library(RStoolbox)
?rasterCVA
cva85_10 <- rasterCVA(teste85[[3:4]], teste2010[[3:4]], tmf = 2)
par(mfrow=c(2,2))
plotRGB(teste85, 3,2,1, stretch='lin')
plotRGB(teste2010, 3,2,1, stretch='lin')
plot(cva85_10[[1]])
plot(cva85_10[[2]])
dev.off()

cva <- cva85_10[[1]]
plot(cva==0)
cva[cva==0] <- NA
plot(cva)

#LandChange ----
#Biomass loss
biomassloss <- cva<=90
biomassloss[biomassloss==0]<-NA
plot(biomassloss)

plotRGB(teste2010, 4,2,1, stretch='lin')
plot(biomassloss, add=TRUE, legend=FALSE)

# Forest cleaning
forestclearing <- ((cva<=180)-biomassloss)
forestclearing[forestclearing==0] <- NA
plot(forestclearing)

plotRGB(teste2010, 4,2,1, stretch='lin')
plot(forestclearing, add=TRUE, legend=FALSE)

# biomass gain
biomassgain <- ((cva<=270)-(cva<=180))
biomassgain[biomassgain==0]<-NA
plot(biomassgain)

plotRGB(teste2010, 4,2,1, stretch='lin')
plot(biomassgain, add=TRUE, legend=FALSE)

# Regrowth
regrowth <-  ((cva<=360)-(cva<=270))
regrowth[regrowth==0]<-NA
plot(regrowth)

plotRGB(teste2010, 4,2,1, stretch='lin')
plot(regrowth, add=TRUE, legend=FALSE)