# this code was created in R version 3.6.0 
## Cecilia O'Leary
## 16 July 2019

## this code takes a shapefile and converts that into a grid to be used in VAST to create your own extrapolation grid
## input: minimum is shapefile that covers region of observations with latitude and longitudes of strata/regions
## output: a single data frame with columns (Lat, Lon, area of each grid cell (row)) 

## generate raster from shapefile
library(sp) #version sp_1.3-1    
library(maps) #version maps_3.3.0
library(mapdata) #version mapdata_2.3.0
library(raster) #version  raster_2.8-19
library(sf) #version sf_0.7-4
library(rgdal) #version rgdal_1.4-3
library(rgeos)

getwd()
working_directory <- "D:/Projects/VAST/V3.0/GeoDBs"

##load shapefile
fgdb <- "C:/Users/jon.richar/Work/GitRepos/VASTModels/GeoDB/VAST_final_geo.gdb"
subset(ogrDrivers(), grepl("GDB", name))
fc_list <- ogrListLayers(fgdb)
print(fc_list)



strata <- readOGR(dsn=fgdb,layer="Pribilof_District_BKC_wEast9_No_Land")
##transform shapefile into decimal degrees
strata <- spTransform(strata, CRS("+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))

plot(strata,axes=T)
strata@data

lon <- sum(bbox(strata)[1,])/2 #retrieves spatial bounding box from spatial data [,1] is longitude
utmzone <- floor((lon + 180) / 6) + 1 #convert decimal degrees to utm zone for average longitude, use for new CRS

crs1 <- strata@proj4string
crs1 <- CRS('+proj=longlat +ellps=WGS84 +no_defs') 
strata@proj4string <- crs1

crs2 <- CRS(paste0("+proj=utm +zone=",utmzone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
strata1 <- spTransform(strata, crs2) #convert strata to utm

plot(strata1,axes=T)


CELL_SIZE = 5000 #size of grid in meters (if working in UTM)
#####sf method###########
#########################
st_bbox(strata1) #get bounding box info
grid <- st_make_grid(strata1,cellsize = CELL_SIZE, what = "centers") #cellsize in meters, creates grid using sf() package = create a square or hexagonal grid over the bounding box of an sf or sfc object
plot(strata1,axes=T)
plot(grid, axes = TRUE, add = TRUE, col = "red")
grid_spatial <- as(grid, "Spatial") #convert grid to Spatial Points
grid_spatial_2 <- as(grid_spatial, "SpatialPointsDataFrame") #convert Spatial Points to SpatialPointsDataFrame
grid_spatial_2@data <- over(grid_spatial,strata1) #combine shapefile data (strata1) with Spatial Points (grid_spatial) & place in SpatialPointsDataFrame data 
names(grid_spatial_2@data)

# (this provides you with your strata identifier (here called Id) in your data frame) 

#put grid in correct/original CRS
grid1 <- spTransform(grid_spatial_2, crs1) 
grid1 <- as.data.frame(grid1) #confirm it's still a data frame
names(grid1)
grid2 <- with(grid1,data.frame(Lon=-(180-coords.x1),Lat=coords.x2,Id=Dissolve_field, Area_km2=(CELL_SIZE/1000*CELL_SIZE/1000),1:nrow(grid1)))
grid2 <- subset(grid2,!is.na(Id)) #this is to filter out the grid that does not overlap with your region. REGION can be

#replaced by whatever your shapefile@data identifier is for your areas/stratas/regions that contain 
#observations

names(grid2) <- c("Lon","Lat","id","Area_km2","row_number")
names(grid2)
saveRDS(grid2, file = "C:/Users/jon.richar/Work/GitRepos/VASTModels/ExtrapolationAreas/Prib_BKC_extrapolation_areas/Prib_BKC_5km_grid.rds")
