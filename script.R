#rm(list=ls())
#gc()

library(raster)
library(rgdal)

# Script used to:
# 1) Prepare the Stelite Images and Indexes (NDVI)
# 2) Run unsupervised land cover classification
# For Restoration mutirão project

# Infos about LAndSat 5 -----
# source (http://www.dgi.inpe.br/siteDgi/ATUS_LandSat.php)
# Source (http://landsat.usgs.gov/documents/si_product_guide.pdf)
# Banda faixa espectral
#1 0,45 a 0,52 µm - azul
#2 0,52 a 0,60 µm - verde
#3 0,63 a 0,69 µm - vermelho
#4 0,76 a 0,90 µm - infravermelho próximo
#5 1,55 a 1,75 µm - infravermelho médio
#6 10,4 a 12,5 µm - infravermelho termal
#7 2,08 a 2,35 µm - infravermelho distante

# Loading and defining study area ----
rj <- readOGR(dsn='/home/felipe/Projetos/R_RS/shape/', layer='33MUE250GC_SIR')
#rj <- readOGR(dsn='./', layer='Pol_teste')
rj <- rj[rj@data$NM_MUNICIP=='RIO DE JANEIRO',]

# 1) Preparing landsat imagery ----
source('../SegmentationFCT/LULCC_fct.R')
args(satPrep)

# 1985 
#Image85 <- satPrep(rasterFolder = "./ModData/1985/", crop.Ext = rj, ndvi = TRUE, evi = TRUE, savi = TRUE, nameOutput = "stack_1985")
Image85 <- stack('./ModData/1985/stack_1985.tif')
# plotRGB(Image85, r=8,g=2,b=1, stretch='lin')

# 1991
#Image91 <- satPrep(rasterFolder = "./ModData/1991/", crop.Ext = rj, ndvi = TRUE, evi = TRUE, savi = TRUE, nameOutput = "stack_1991")
Image91 <- stack('./ModData/1991/stack_1991.tif')

# 1995
#Image95 <- satPrep(rasterFolder = "./ModData/1995/", crop.Ext = rj, ndvi = TRUE, evi = TRUE, savi = TRUE, nameOutput = "stack_1995")
Image95 <- stack('./ModData/1995/stack_1995.tif')

# 2000
#Image00 <- satPrep(rasterFolder = "./ModData/2000/", crop.Ext = rj, ndvi = TRUE, evi = TRUE, savi = TRUE, nameOutput = "stack_2000")
Image00 <- stack('./ModData/2000/stack_2000.tif')

# 2005
#Image05 <- satPrep(rasterFolder = "./ModData/2005/", crop.Ext = rj, ndvi = TRUE, evi = TRUE, savi = TRUE, nameOutput = "stack_2005")
Image05 <- stack('./ModData/2005/stack_2005.tif')

# 2010
#Image2010 <- satPrep(rasterFolder = "./ModData/2010/", crop.Ext = rj, ndvi = TRUE, evi = TRUE, savi = TRUE, nameOutput = "stack_2010")
Image2010 <- stack('./ModData/2010/stack_2010.tif')

#par(mfrow=c(2,1))
#plotRGB(teste85, r=8,g=2,b=1, stretch='lin')
#plotRGB(Image2010, r=8,g=2,b=1, stretch='lin')

# 1.2 Running PCA ----
source('~/Projetos/ARF_spatial_planning_old/ENM/fct/eigenvariables.fct.R')
args(eigenvariables.fct)
#View(eigenvariables.fct)

#pca85 <- eigenvariables.fct(Image85, 'pca85', 95)
pca85 <- stack('./env/pca85.eigenvariables.grd')
#pca91 <- eigenvariables.fct(Image91, 'pca91', 95)
pca91 <- stack('./env/pca91.eigenvariables.grd')
#pca95 <- eigenvariables.fct(Image95, 'pca95', 95)
pca95 <- stack('./env/pca95.eigenvariables.grd')
#pca00 <- eigenvariables.fct(Image00, 'pca00', 95)
pca00 <- stack('./env/pca00.eigenvariables.grd')
#pca05 <- eigenvariables.fct(Image05, 'pca05', 95)
pca05 <- stack('./env/pca05.eigenvariables.grd')
#pca10 <- eigenvariables.fct(Image2010, 'pca10', 95)
pca10 <- stack('./env/pca10.eigenvariables.grd')

# 2)Unsupervised classification ----
source('~/Projetos/SegmentationFCT/segmentation.R')
#View(segmentation)
args(segmentation)

# Exploratory Analysis of whithin sums of squares -----
#segmentation(envLayer = pca85, studyArea = rj, projName = "RJ_1985_PCA", folder = './ModData/1985/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = FALSE, random.pt = NULL, ngroup = 10, save.shp = FALSE, save.raster = FALSE, explore = TRUE, h.life=TRUE, save.fit = TRUE, seed = 123)

#segmentation(envLayer = pca91, studyArea = rj, projName = "RJ_1991_PCA", folder = './ModData/1991/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = FALSE, random.pt = NULL, ngroup = 10, save.shp = FALSE, save.raster = FALSE, explore = TRUE, h.life=TRUE, save.fit = TRUE, seed = 123)

#segmentation(envLayer = pca95, studyArea = rj, projName = "RJ_1995_PCA", folder = './ModData/1995/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = FALSE, random.pt = NULL, ngroup = 10, save.shp = FALSE, save.raster = FALSE, explore = TRUE, h.life=TRUE, save.fit = TRUE, seed = 123)

#segmentation(envLayer = pca00, studyArea = rj, projName = "RJ_2000_PCA", folder = './ModData/2000/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = FALSE, random.pt = NULL, ngroup = 10, save.shp = FALSE, save.raster = FALSE, explore = TRUE, h.life=TRUE, save.fit = TRUE, seed = 123)

#segmentation(envLayer = pca05, studyArea = rj, projName = "RJ_2005_PCA", folder = './ModData/2005/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = FALSE, random.pt = NULL, ngroup = 10, save.shp = FALSE, save.raster = FALSE, explore = TRUE, h.life=TRUE, save.fit = TRUE, seed = 123)

#segmentation(envLayer = pca10, studyArea = rj, projName = "RJ_2010_PCA", folder = './ModData/2010/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = FALSE, random.pt = NULL, ngroup = 10, save.shp = FALSE, save.raster = FALSE, explore = TRUE, h.life=TRUE, save.fit = TRUE, seed = 123)

# Segmentation considering exploratory anaçysis previusly run ----

#fuzzy85 <- segmentation(envLayer = pca85, studyArea = rj, projName = "RJ_1985_PCA", folder = './ModData/1985/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = TRUE, random.pt = NULL, ngroup = 18, save.shp = FALSE, save.raster = TRUE, explore = FALSE, h.life = FALSE, save.fit = TRUE, seed = 123)

fuzzy85 <- stack('./ModData/1985/fuzzy_segmentation_RJ_1985.tif') 
load('./ModData/1985/RJ_1985_Fuzzy_.RData')
fit85 <- fuzzy
str(fit85)
fit85$center
# writting distance between cluster centers
write.csv(as.matrix(dist(fit85$centers, method = "euclidean")), './ModData/1985/1985Dists.csv')

#fuzzy91 <- segmentation(envLayer = pca91, studyArea = rj, projName = "RJ_1991_PCA", folder = './ModData/1991/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = TRUE, random.pt = NULL, ngroup = 16, save.shp = FALSE, save.raster = TRUE, explore = FALSE, h.life = FALSE, save.fit = TRUE, seed = 123) 

fuzzy91 <- stack('./ModData/1991/fuzzy_segmentation_RJ_1991.tif') 
load('./ModData/1991/RJ_1991_Fuzzy_.RData')
fit91 <- fuzzy
# writting distance between cluster centers
write.csv(as.matrix(dist(fit91$centers, method = "euclidean")), './ModData/1991/1991Dists.csv')

#fuzzy95 <- segmentation(envLayer = pca95, studyArea = rj, projName = "RJ_1995_PCA", folder = './ModData/1995/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = TRUE, random.pt = NULL, ngroup = 18, save.shp = FALSE, save.raster = TRUE, explore = FALSE, h.life = FALSE, save.fit = TRUE, seed = 123) 
fuzzy95 <- stack('./ModData/1995/fuzzy_segmentation_RJ_1995.tif') 
load('./ModData/1995/RJ_1995_Fuzzy_.RData')
fit95 <- fuzzy
# writting distance between cluster centers
write.csv(as.matrix(dist(fit95$centers, method = "euclidean")), './ModData/1995/1995Dists.csv')

#fuzzy00 <- segmentation(envLayer = pca00, studyArea = rj, projName = "RJ_2000_PCA", folder = './ModData/2000/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = TRUE, random.pt = NULL, ngroup = 19, save.shp = FALSE, save.raster = TRUE, explore = FALSE, h.life = FALSE, save.fit = TRUE, seed = 123) 
fuzzy00 <- stack('./ModData/2000/fuzzy_segmentation_RJ_2000.tif') 
load('./ModData/2000/RJ_2000_Fuzzy_.RData')
fit00 <- fuzzy
# writting distance between cluster centers
write.csv(as.matrix(dist(fit00$centers, method = "euclidean")), './ModData/2000/2000Dists.csv')

#fuzzy05 <- segmentation(envLayer = pca05, studyArea = rj, projName = "RJ_2005_PCA", folder = './ModData/2005/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = TRUE, random.pt = NULL, ngroup = 15, save.shp = FALSE, save.raster = TRUE, explore = FALSE, h.life = FALSE, save.fit = TRUE, seed = 123) 
fuzzy05 <- stack('./ModData/2005/fuzzy_segmentation_RJ_2005.tif') 
load('./ModData/2005/RJ_2005_Fuzzy_.RData')
fit05 <- fuzzy
# writting distance between cluster centers
write.csv(as.matrix(dist(fit05$centers, method = "euclidean")), './ModData/2005/2005Dists.csv')

#fuzzy10 <- segmentation(envLayer = pca10, studyArea = rj, projName = "RJ_2010_PCA", folder = './ModData/2010/', randomforest = FALSE, Kmeans = FALSE, fuzzy.cluster = TRUE, random.pt = NULL, ngroup = 16, save.shp = FALSE, save.raster = TRUE, explore = FALSE, h.life = FALSE, save.fit = TRUE, seed = 123) 
fuzzy10 <- stack('./ModData/2010/fuzzy_segmentation_RJ_2010.tif') 
load('./ModData/2010/RJ_2010_Fuzzy_.RData')
fit10 <- fuzzy

# writting distance between cluster centers
write.csv(as.matrix(dist(fit10$centers, method = "euclidean")), './ModData/2010/2010Dists.csv')

# ----