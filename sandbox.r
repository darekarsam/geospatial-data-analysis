# WGS84 proj4string
wgs84 = "+init=epsg:4326"

# load the point data as a dataframe object
bosbiz.df = read.csv(file = "bosbiz.csv", as.is = T)
#2939 missing values in SIC codes
sum(is.na(bosbiz.df$sic))
#ggplot removes rows with missing values hence replace missing values with 0
bosbiz.df[is.na(bosbiz.df)] <- 0

head(bosbiz.df)

bosbiz.df[duplicated(bosbiz.df),]
bosbiz.df <- bosbiz.df[!duplicated(bosbiz.df),]

counts <- table(bosbiz.df$sic)
counts <- as.data.frame(counts)
filtered.counts <- counts[counts$Freq > 100,]
filtered.counts

# convert to a Spatial Points dataframe
library(sp)
library(rgdal)
bosbiz.sp = SpatialPointsDataFrame(coords = bosbiz.df[,c("lon", "lat")], data = bosbiz.df[,1:2], proj4string = CRS(wgs84))
head(bosbiz.sp)
summary(bosbiz.sp)

cov.sp(coords = bosbiz.sp, sp.type = "exponential", sp.par = c(2, 1),
       error.var = 1)
#distance matrix distGeo function in https://cran.r-project.org/web/packages/geosphere/geosphere.pdf
library(geosphere)
library(dplyr)

mdist <- distm(bosbiz.sp)
typeof(mdist)
#converting to Km
mdist <- mdist/1000
mdist[1,]
#mdist <- as.dist(mdist)

save(mdist, file="mdist.rda")
rm(mdist)
rm(bosbiz.sp)
************************************
load("mdist.rda")
library(dbscan)
dbscan::kNNdistplot(mdist, k = 6)
abline(h = 9, lty = 2)
6 - 10
8 -10
10 -11  
clusters <- dbscan(mdist, eps=9, minPts = 6)


save(clusters1, file="dbscanClusters.rda")
rm(clusters1)

load("dbscan.rda")

clusters1
mdist[1]

mat <- select(bosbiz.df, lon, lat)
hullplot(mat, clusters1$cluster)

wmdist[2]

mat <- select(bosbiz.df, lon, lat)
head(mat)
distGeo(mat, mat, a=6378137, f=1/298.257223563)

dists <- spDists(bosbiz.sp,bosbiz.sp,longlat = TRUE)
dists.mat <- spDists(bosbiz.sp, longlat = TRUE)
dim(dists.mat) <- c(length(bosbiz.sp), length(bosbiz.sp))



dist.matrix <- matrix(dists.mat[1,], length(bosbiz.sp), byrow=TRUE)
dists.mat[3]

bosbiz.sp[1,]


library(dbscan)
library(dplyr)
dbscan::kNNdistplot(select(bosbiz.df, lat, lon), k =  6)

clusters1 <- dbscan(dists.mat, eps=2, minPts = 10 )
clusters1

mat <- select(bosbiz.df, lon, lat)
hullplot(dists.mat, clusters1$cluster)

maxLat <- max(bosbiz.df$lat)
maxLon <- max(bosbiz.df$lon)
minLat <- min(bosbiz.df$lat)
minLon <- min(bosbiz.df$lon)

#install.packages('ggmap')
#install.packages('rgdal')

#proj4string(bosbiz.sp) <- CRS(wgs84)                   # set original projection
#bosbiz.sp <- spTransform(bosbiz.sp, CRS("+proj=longlat +datum=WGS84"))

# other helpful libraries
#libs = c("rgdal", "rgeos", "maptools", "spatstat", "dbscan")
#lapply(libs, library, character.only = T)

head(bosbiz.sp)
plot(bosbiz.sp)

head(bosbiz.df)

library(ggmap)
theme_update(plot.title = element_text(hjust = 0.5)) #for bringing tile to center

boston.bounds <- c(left = -71.19, bottom = 42.22, right = -70.93, top = 42.45)
map <- get_map(boston.bounds, zoom = 12,source="google",maptype="roadmap")
ggmap(map) + geom_point(data=as.data.frame(bosbiz.df), aes(x=lon, y=lat),alpha=0.2, color="red") +
  ggtitle("Boston City Map")

#Dynamic bounds with varying color depth as per density
boston.bounds <- c(left = minLon, bottom = minLat, right = maxLon, top = maxLat)
map <- get_map(boston.bounds, zoom = 13, scale= 2,source="google",maptype="roadmap")
ggmap(map) + geom_point(data=as.data.frame(bosbiz.df), aes(x=lon, y=lat), alpha=0.2, color='red') + 
  ggtitle("Boston Map boundary adjusted")

library(dbscan)
library(dplyr)
#getting optimal eps
dbscan::kNNdistplot(select(bosbiz.df, lat, lon), k =  6)+
abline(h = 0.0015, lty = 2) + title("Optimal value of eps for DBSCAN using KNN")

mat <- select(bosbiz.df, lon, lat)
clusters <- dbscan(mat, eps=0.0015, minPts = 6 )
clusters
ggmap(map)+ geom_point(aes(x=lon, y=lat), data=mat, col = clusters$cluster + 1, alpha=0.5, pch = 20)

plot(lat ~ lon, data = mat, col = clusters$cluster + 1, pch = 20)

hullplot(mat, clusters$cluster)

#Calculate LOF (local outlier factor) and visualize (larger bubbles in the visualization have a larger LOF)
lof <- lof(mat, k = 6)
lof
plot(mat, pch = ".", main = "LOF (k=6)")
points(mat, cex = (lof-1), pch = 1, col="red")
text(mat[lof>2,], labels = round(lof, 1)[lof>2], pos = 3)

kmeans.result <- kmeans(mat, 6)
kmeans
plot(mat, col=kmeans.result)

#Run OPTICS
opt <- optics(mat, eps=0.0015, minPts = 6)
#Extract DBSCAN-like clustering from OPTICS and create a reachability plot (extracted DBSCAN clusters at eps_cl=.4 are colored)
opt <- extractDBSCAN(opt, eps_cl = 0.0015)
opt
plot(opt)

library(dbscan)
#Extract a hierarchical clustering using the Xi method (captures clusters of varying density)
opt <- extractXi(opt, xi = 0.0015)
opt
plot(opt)

#Run HDBSCAN (captures stable clusters)
hdb <- hdbscan(mat, minPts = 15)

plot(hdb, show_flat = T)

hdb$membership_prob

#See how well each point corresponds to the clusters found by the model used

colors <- mapply(function(col, i) adjustcolor(col, alpha.f = hdb$membership_prob[i]), 
                 palette()[hdb$cluster+1], seq_along(hdb$cluster))
plot(mat, col=colors, pch=20)

#library(ggplot2)
#library(magrittr)
plot(select(bosbiz.df, lat, lon), clusters$cluster)

dbscan::kNNdistplot(select(bosbiz.df, lat, lon), k =  6)
abline(h = 0.0015, lty = 2)

groups <- clusters %>% filter(cluster != 0) #clusters[clusters$cluster!=0]
clusters
n <- c(1,2,3)
n

library(dplyr)
bosbiz.filt <- bosbiz.df[bosbiz.df$sic==80119904,]
bosbiz.filt
mat.filt <-  select(bosbiz.filt, lon, lat)
mat.filt

library(ggmap)
boston.bounds <- c(left = minLon, bottom = minLat, right = maxLon, top = maxLat)
map <- get_map(boston.bounds, zoom = 13, scale= 2,source="google",maptype="roadmap")
ggmap(map) + geom_point(data=as.data.frame(bosbiz.filt), aes(x=lon, y=lat), alpha=0.2, color='red') + 
  ggtitle("Boston Map boundary adjusted")
