rm(list=ls())#clear environment

#Load libraries
library(raster)
library(sf)
library(sp)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(lwgeom)
library(reshape)
library(stars)
library(broom)
library(viridis)
library(ggpubr)
library(ggnewscale)
library(Metrics)
library(terra)

#Parameters
Imagery_filepath <- "C:/Users/OCRONING/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/TickTickBloom/Inputs/CyAN_S3_imagery/Aug_2019/"
AOI_filepath <- "C:/Users/OCRONING/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/TickTickBloom/Inputs/Shapefiles/"
Figure_filepath <- "C:/Users/OCRONING/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/TickTickBloom/Outputs/Figures/"
Table_filepath <- "C:/Users/OCRONING/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/TickTickBloom/Outputs/Tables/"
S3_row <- 5
S3_col <- 7
n_samples <- 200

#Set up
S3_Cyan_files <- list.files(path = Imagery_filepath, full.names = T)
S3_Cyan_filenames <- list.files(path = Imagery_filepath)
Year <- substring(S3_Cyan_filenames[1],2,5)
JDates <- as.numeric(substring(S3_Cyan_filenames,6,8)) #first two weeks of august
Dates <- as.Date(JDates, origin=as.Date(paste0(Year,"-01-01")))

#Select Target Sentinel Tile
S3_tiles <- read_sf((paste0(AOI_filepath, "CONUS_Tiles.shp")), layer='CONUS_Tiles', quiet = T)
FL_S3_tile <-  S3_tiles %>% dplyr::filter(row_1 == S3_row & col_1 == S3_col) # Filter for target S3 tile
FL_S3_tile <- st_transform(FL_S3_tile,  4326) # Reproject tile to WGS84

#Read in Resolvable MERIS Lakes Shapefile
MERIS_lakes <- read_sf(paste0(AOI_filepath, "MERIS_OLCI_Lakes.shp"), layer = "MERIS_OLCI_Lakes", quiet = T)
MERIS_lakes <- st_transform(MERIS_lakes,  4326) # Reproject tile to WGS84

#Read in Florida Shapefile
FL_counties <- read_sf(paste0(AOI_filepath, "Florida_Counties.shp"), layer = "Florida_Counties", quiet = T)
FL_counties <- st_transform(FL_counties,  4326) # Reproject tile to WGS84
FL <- st_union(FL_counties) # Dissolve counties to make the shapefile be the boundaries of just Florida state

#Filter and merge MERIS lakes that fall within Florida
sf::sf_use_s2(FALSE)
FL_lakes <- st_intersection(st_sf(MERIS_lakes), st_sf(FL)) #Filter and merge MERIS lakes that fall inside Florida
FL_lakes <- st_union(FL_lakes) #Dissolve lakes to make one polygon

#Filter MERIS lakes that fall inside S3 Tile 7-5, extract the centroid of each lake, and merge lakes
S3_lakes <- st_intersection(st_sf(MERIS_lakes), st_sf(FL_S3_tile))
S3_lake_centroids <- st_centroid(S3_lakes)
S3_lakes <- st_union(S3_lakes)

#Create centroid dataframe
S3_lake_centroids_df <- as.data.frame(coordinates(as_Spatial(S3_lake_centroids)))
names(S3_lake_centroids_df) <- c("longitude", "latitude")

#Extract stratified random point sample to supplement the centroids so number of samples matches target
Extra_pts <- as.data.frame(st_coordinates(st_transform(st_sample(S3_lakes, (n_samples - nrow(S3_lake_centroids_df)), "stratified"), 4326)))
names(Extra_pts) <- c("longitude", "latitude")

#Merge centroid and extra points to make final sampling dataframe
S3_lake_centroids_df <- as.data.frame(rbind(S3_lake_centroids_df, Extra_pts))

############################Mapping##############################

#Read in One CyAN raster
CyAN_1 <- raster(S3_Cyan_files[2], band=1)

#Reclassify raster to reflect flags descibed on CyAN site
m <- c(-Inf, 0, 1,  0, 253, 2,  253, 254, 3,  254, 255, 4)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(CyAN_1, rclmat)
pol <- sf::as_Spatial(sf::st_as_sf(stars::st_as_stars(rc), as_points = FALSE, merge = TRUE))
colnames(pol@data)[1] <- 'raster'

TooLow <- st_as_sf(pol[pol$raster == 1,])
Land <- st_as_sf(pol[pol$raster == 3,])
NoData <- st_as_sf(pol[pol$raster == 4,])

#Convert raster to df for mapping
CyAN_2 <- as.data.frame(CyAN_1, xy = TRUE) 
names(CyAN_2) <- c("long", "lat", "CyAN")
CyAN_2$Cell_Values <- seq(1,nrow(CyAN_2),1)

#Extract water values
Water <- na.omit(CyAN_2[CyAN_2$CyAN < 254 & CyAN_2$CyAN > 0,])

#Project AOI files to CyAN projection   
FL_S3_tile <- st_transform(FL_S3_tile, crs(CyAN_1))
FL <- st_transform(FL, crs(CyAN_1))

######################Extract CyAN Values############################

#Create sf object and df with same projection as CyAn raster out of sampling schema df
coordinates(S3_lake_centroids_df) <- ~ longitude + latitude
S3_lake_centroids_df <- st_as_sf(S3_lake_centroids_df)
S3_lake_centroids_df <- st_set_crs(S3_lake_centroids_df, 4326 )

S3_lake_centroids_pts <- st_transform(st_as_sf(S3_lake_centroids_df), crs(CyAN_1))
S3_lake_centroids_pts <- as.data.frame(coordinates(as_Spatial(S3_lake_centroids_pts)))
names(S3_lake_centroids_pts) <- c("longitude", "latitude")

coords <- as.data.frame(cbind(S3_lake_centroids_pts$longitude, S3_lake_centroids_pts$latitude))

#Extract CyAN values from all rasters
CyAN_extract <- list()
for(i in 1:length(S3_Cyan_files)){
  CyAN_raster <- raster(S3_Cyan_files[i], band = 1)
  CyAN_extract[[i]] <- raster::extract(CyAN_raster, coords, fun = mean)
}

#Turn onto melted df
CyAN <- as.data.frame(do.call(cbind, CyAN_extract))
CyAN <- cbind(as.data.frame(coords), CyAN)
names(CyAN) <- c("long", "lat", JDates)

CyAN_melt <- na.omit(melt(CyAN, id=c("long","lat")))
names(CyAN_melt) <- c("long", "lat", "JDate", "CyAN_DN")
CyAN_melt$JDate <- as.numeric(as.character(CyAN_melt$JDate))
CyAN_melt$Date <- as.Date(CyAN_melt$JDate, origin=as.Date(paste0(Year,"-01-01")))
CyAN_melt <- CyAN_melt[,c("long", "lat", "Date", "CyAN_DN")]

#Classify each point with its flag
CyAN_melt$Category <-
  ifelse(CyAN_melt$CyAN_DN == 0, "Below threshold of CI detection limits",
         ifelse(CyAN_melt$CyAN_DN == 254, "Land",
                ifelse(CyAN_melt$CyAN_DN == 255, "No Data",
                       ifelse(CyAN_melt$CyAN_DN > 0 & CyAN_melt$CyAN_DN < 254, "CyAN Data", NA)
                )
         )
  )

#Convert DN of water values to CI
CyAN_DN_Coversion <- function(DN){ (10^((DN*0.011714)-4.1870866)) * 1.0E+8}

CyAN_melt$CyAN_CI_cells_ml <-  ifelse(CyAN_melt$CyAN_DN == 0, NA,
                                      ifelse(CyAN_melt$CyAN_DN == 254, NA,
                                             ifelse(CyAN_melt$CyAN_DN == 255, NA,
                                                    ifelse(CyAN_melt$CyAN_DN > 0 & CyAN_melt$CyAN_DN < 254, CyAN_DN_Coversion(CyAN_melt$CyAN_DN) , NA)
                                             )
                                      )
)

Water$CyAN_CI_cells_ml <- CyAN_DN_Coversion(Water$CyAN)

CyAN_melt <- CyAN_melt[,c("long", "lat", "Date", "Category", "CyAN_DN", "CyAN_CI_cells_ml")]
names(CyAN_melt) <- c("longitude", "latitude", "date", "Category", "CyAN_DN", "CyAN_CI_cells_ml")

#Create summary table of number of each cell type 
CyAN_Summary <- as.data.frame(table(CyAN_melt$Category))
names(CyAN_Summary) <- c("Flag", "CyAN_Count")
write.csv(CyAN_Summary, paste0(Table_filepath, "CyAN_Flag_Summary.csv"), row.names = F)

#Export and map points for CyFi processing
CyFi_pts <- na.omit(CyAN_melt)
CyFi_pts <- CyFi_pts[c("latitude", "longitude", "date")]
CyFi_data <- CyFi_pts[c(2,1,3)]

my.sf.point <- st_as_sf(x = CyFi_data, 
                        coords = c("longitude", "latitude"),
                        crs = proj4string(CyAN_raster))

sf_trans <- st_transform(my.sf.point, 4326)
sf_trans <- as.data.frame(st_coordinates(sf_trans))
sf_trans <- as.data.frame(cbind(sf_trans, CyFi_data$date))
names(sf_trans) <- c("longitude","latitude", "date")
sf_trans <- sf_trans[,c(2,1,3)]

write.csv(sf_trans, paste0(Table_filepath, "CyAN_to_CyFi_Points.csv"), row.names = F)

Flag_df <- (rbind(TooLow, Land, NoData))
names(Flag_df) <- c("Flag", "geometry")

#Plot Study Area with points
site_map <- ggplot() +
  geom_sf(data = Flag_df, aes(fill = factor(Flag)), color = NA) + 
  scale_fill_manual(values = c("1" = "grey30", "3" = "#9b8d73", "4" = "grey80"), 
                    name= "Flag", labels = c("Too Low", "Land", "No data")) +
  # scale_color_manual(values = c("1" = "green", "3" = NA, "4" = NA)) +
  geom_sf(data = FL_S3_tile, color = "black", 
          fill = NA, lwd = 0.5, alpha = 0.5) +
  geom_sf(data = FL, color = "black", fill = NA, lwd = 0.4)+
  ggnewscale::new_scale_fill() +
  geom_tile(data = Water,
            aes(x = long, y = lat, fill = (CyAN))) +
  scale_fill_viridis(name="Cyano DN",  limits=c(1, 253),
                     values=c(0,1)) +
  ggnewscale::new_scale_fill()+
  geom_point(data = CyFi_pts, 
             aes(x = longitude, y = latitude), shape = 21,
             color = "darkred", fill = NA, stroke = 0.25,
             size = 0.75, alpha = 0.7) +
  coord_sf(datum = sf::st_crs(4326)) +
  ggtitle('Sample Site Locations') + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 1,
                                   hjust=1),
        axis.title = element_blank(),
        legend.direction = "vertical", 
        legend.box = "vertical")
site_map
ggsave(filename="Site_Sampling_Map.pdf", plot = site_map, device = "pdf", path = Figure_filepath, width = 6, height = 5, dpi = 300, units = "in")

#Number of points per date
CyAN_Date_Counts <- as.data.frame(table(CyFi_pts$date))
CyAN_Dates_All <- as.data.frame(table(CyAN_melt$date))
CyAN_Dates_All$Freq <- NA
CyAN_Date_Counts <- merge(CyAN_Dates_All, CyAN_Date_Counts, by="Var1", all.x = T)
CyAN_Date_Counts <- CyAN_Date_Counts[,c(1,3)]
CyAN_Date_Counts[is.na(CyAN_Date_Counts)] <- 0
names(CyAN_Date_Counts) <- c("Date", "Frequency")

write.csv(CyAN_Date_Counts, paste0(Table_filepath, "CyAN_Date_Counts.csv"), row.names = F)

rm(CyAN_Dates_All)

#Create CyAn df
CyAN_Data <- na.omit(CyAN_melt)
CyAN_Data <- CyAN_Data[c("latitude", "longitude", "date", "CyAN_CI_cells_ml")]

my.sf.point <- st_as_sf(x = CyAN_Data, 
                        coords = c("longitude", "latitude"),
                        crs = proj4string(CyAN_raster))

sf_trans <- st_transform(my.sf.point, 4623)
sf_trans <- as.data.frame(st_coordinates(sf_trans))
sf_trans <- as.data.frame(cbind(sf_trans, CyAN_Data$date, CyAN_Data$CyAN_CI_cells_ml))
names(sf_trans) <- c("longitude","latitude", "date", "CyFi_cells_ml")
CyAN_Data <- sf_trans[,c(2,1,3,4)]
write.csv(CyAN_Data, paste0(Table_filepath, "CyAN_Cyano_Estimates.csv"), row.names = F)
