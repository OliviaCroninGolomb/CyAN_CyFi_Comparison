rm(list=ls())#clear environment
#args <- commandArgs(trailingOnly = TRUE)

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
#library(reticulate)

#Parameters
Imagery_filepath <- "C:/Users/OCRONING/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/TickTickBloom/Inputs/CyAN_S3_imagery/Aug_2019/"
CyFi_filepath <- "C:/Users/OCRONING/OneDrive - Environmental Protection Agency (EPA)/Profile/Documents/TickTickBloom/Inputs/CyFi_Data/"
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

############################################################################

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

Flag_df <- (rbind(TooLow, Land, NoData))
names(Flag_df) <- c("Flag", "geometry")

#Convert raster to df for mapping
CyAN_2 <- as.data.frame(CyAN_1, xy = TRUE) 
names(CyAN_2) <- c("long", "lat", "CyAN")
CyAN_2$Cell_Values <- seq(1,nrow(CyAN_2),1)

#Extract water values
Water <- na.omit(CyAN_2[CyAN_2$CyAN < 254 & CyAN_2$CyAN > 0,])

#Project AOI files to CyAN projection   
FL_S3_tile <- st_transform(FL_S3_tile, crs(CyAN_1))
FL <- st_transform(FL, crs(CyAN_1))

#################################CyFi###########################################

#Create CyFi Dataset
CyFi_Data <- read.csv(paste0(CyFi_filepath, "preds.csv") )
CyFi_meta <- read.csv(paste0(CyFi_filepath, "sentinel_metadata.csv") )
CyFi_meta <- CyFi_meta[CyFi_meta$days_before_sample <= 7,] #select only samples within a week on a sentinel 2 image
CyFi_meta <- as.data.frame(CyFi_meta$sample_id)
names(CyFi_meta) <- "sample_id"
CyFi_Data <- merge(CyFi_meta, CyFi_Data, all = T, by=c("sample_id"))

CyFi_Data <- as.data.frame(na.omit(CyFi_Data[,c(3,4,2,5)]))
names(CyFi_Data) <- c("latitude", "longitude", "date", "CyFi_cells_ml")

#################################CyAN###########################################

#Read in CyAN Dataset
CyAN_Data <- read.csv(paste0(Table_filepath, "CyAN_Cyano_Estimates.csv"))

#################################CyFi_CyAN###########################################

#Merge CyAN and CyFi dfs
Cyano_Data <- na.omit(merge(CyAN_Data, CyFi_Data, all = T, by=c("latitude", "longitude", "date")))
names(Cyano_Data) <- c("latitude", "longitude", "date", "CyAN_CI_cells_ml", "CyFi_cells_ml")
write.csv(Cyano_Data, paste0(Table_filepath, "Cyano_Data.csv"), row.names = F)

#Number of points per date
CyAN_Date_Counts <- as.data.frame(table(Cyano_Data$date))
CyAN_Dates_All <- as.data.frame(table(CyAN_Data$date))
CyAN_Dates_All$Freq <- NA
CyAN_Date_Counts <- merge(CyAN_Dates_All, CyAN_Date_Counts, by="Var1", all.x = T)
CyAN_Date_Counts <- CyAN_Date_Counts[,c(1,3)]
CyAN_Date_Counts[is.na(CyAN_Date_Counts)] <- 0
names(CyAN_Date_Counts) <- c("Date", "Frequency")
write.csv(CyAN_Date_Counts, paste0(Table_filepath, "CyAN_Date_Counts.csv"), row.names = F)

#Linear Regression Model
model <- lm(Cyano_Data$CyFi_cells_ml~Cyano_Data$CyAN_CI_cells_ml)
summary <- summary(model)
ar2 <- round(summary$adj.r.squared,2)
b <- round(summary$coefficients[1],2)
m <- round(summary$coefficients[2],2)
MAE <- mae(Cyano_Data$CyAN_CI_cells_ml,Cyano_Data$CyFi_cells_ml)

#Plot Regression
cc_reg <- ggplot(Cyano_Data, aes(x=CyAN_CI_cells_ml, y=CyFi_cells_ml)) +
  geom_point() +
  stat_smooth(method = lm, fullrange = T) +
  geom_abline() +
  annotate(x = 1000000, y = 1700000, label = bquote("y" == .(m) ~ "x + " ~ .(b)), geom = "text", size = 4) +
  annotate(x = 1000000, y = 1600000, label = bquote(R^2 == .(ar2)), geom = "text", size = 4) +
  annotate(x = 1000000, y = 1450000, label = bquote("MAE" == .(MAE)), geom = "text", size = 4) +
  xlim(0,1750000) +
  ylim(0,1750000) +
  ggtitle("Cyanobacteria Abundance in S3 7-5") +
  xlab("CyAN (cells/ml)") +
  ylab("CyFi (cells/ml)") +
  theme(plot.title = element_text(hjust = 0.5, size=15))
cc_reg
ggsave(filename="CyAN_CyFi_Regression.pdf", plot = cc_reg, device = "pdf", path = Figure_filepath, width = 6, height = 5, dpi = 300, units = "in")

#Residuals
resid_df <- broom::augment(model)
cc_resid <- ggplot(resid_df, aes(x = .fitted, y = .resid)) + 
  geom_point() +
  geom_hline(yintercept=0) +
  xlab("Predicted CyFi Values") +
  ylab("Residuals") +
  ggtitle("Cyanobacteria Abundance LM Residuals in S3 7-5") +
  theme(plot.title = element_text(hjust = 0.5, size=15))
  
cc_resid
ggsave(filename="CyAN_CyFi_Residual.pdf", plot = cc_resid, device = "pdf", path = Figure_filepath, width = 6, height = 5, dpi = 300, units = "in")

#Difference Histogram
Cyano_Data$Difference <- (Cyano_Data$CyAN_CI_cells_ml - Cyano_Data$CyFi_cells_ml)

cc_diff <- ggplot(Cyano_Data, aes(x=Difference)) + 
  geom_histogram(aes(y =..density..), binwidth=50000, fill="grey20", color="#e9ecef", alpha=0.9) +
  geom_vline(aes(xintercept = mean(Difference)),col='red',lwd=0.75)+
  geom_vline(aes(xintercept = median(Difference)),col='blue',lwd=0.75)+
  stat_function(fun = dnorm, args = list(mean = mean(Cyano_Data$Difference), sd = sd(Cyano_Data$Difference))) +
  annotate(x = 0.75E+6, y = 19E-7, label = paste0("Mean =", round(mean(Cyano_Data$Difference)),0), geom = "text", color = "red", size = 4) +
  annotate(x = 0.75E+6, y = 17E-7, label = paste0("Mean =", round(median(Cyano_Data$Difference)),0), geom = "text", color = "blue", size = 4) +
  ggtitle("Difference between CyAN and CyFi") +
  theme(plot.title = element_text(hjust = 0.5, size = 15))
cc_diff
ggsave(filename="CyAN_CyFi_Diff_Hist.pdf", plot = cc_diff, device = "pdf", path = Figure_filepath, width = 6, height = 5, dpi = 300, units = "in") 

CyAN_Point_Counts_coords <- as.data.frame(cbind(Cyano_Data$latitude, Cyano_Data$longitude))
names(CyAN_Point_Counts_coords) <- c("lat", "long")

CyAN_Point_Counts <- as.data.frame(table(CyAN_Point_Counts_coords$lat, CyAN_Point_Counts_coords$long))
CyAN_Point_Counts$Var1 <- as.numeric(as.character(CyAN_Point_Counts$Var1))
CyAN_Point_Counts$Var2 <- as.numeric(as.character(CyAN_Point_Counts$Var2))
names(CyAN_Point_Counts) <- c("latitude", "longitude", "Freq")
CyAN_Point_Counts$Rank <- rank(-CyAN_Point_Counts$Freq)

Target_TS <- tail(CyAN_Point_Counts[order(CyAN_Point_Counts$Freq),],3)
Target_TS$Group <- seq(1, nrow(Target_TS), 1)

CyAN_merge <- na.omit(merge(Cyano_Data, Target_TS, all = T, by=c("longitude", "latitude")))
CyAN_merge$date <- as.Date(CyAN_merge$date)

cc_ts <- ggplot(data=CyAN_merge) + 
  geom_line(aes(x=date, y=CyAN_CI_cells_ml, group = factor(Group), 
                colour = factor(Group), linetype = "solid")) + 
  geom_point(aes(x=date, y=CyAN_CI_cells_ml, group = factor(Group), 
                 colour = factor(Group)), size = 1) +
  geom_line(aes(x=date, y=CyFi_cells_ml, group = factor(Group), 
                colour = factor(Group), linetype = "dashed")) + 
  geom_point(aes(x=date, y=CyFi_cells_ml, group = factor(Group), 
                 colour = factor(Group)), size = 1) +
  facet_wrap(~Group, ncol=1, strip.position = "right") +
  labs(y= "Cyanobacteria Abundance [cells/ml]", x = NULL) +
  scale_x_date(date_breaks = "day", date_labels = "%m-%d-%Y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        legend.position = "right", 
        strip.text = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  guides(colour = F) +
  scale_linetype_discrete(name = "", labels = c("CyAN", "CyFi")) +
  ggtitle("CyAN vs CyFi Cyanobacteria Estimates Time-Series")
cc_ts

ggsave(filename="CyAN_CyFi_TimeSeries.pdf", plot = cc_ts, device = "pdf", path = Figure_filepath, width = 6, height = 5, dpi = 300, units = "in")

TS_pts <- unique(as.data.frame(cbind(CyAN_merge$longitude, CyAN_merge$latitude)))
names(TS_pts) <- c("longitude", "latitude")

coordinates(TS_pts) <- ~ longitude + latitude
TS_df <- st_as_sf(TS_pts)
TS_df <- st_set_crs(TS_df, crs(CyAN_1))
TS_df <- as.data.frame(coordinates(as_Spatial(TS_df)))
names(TS_df) <- c("longitude", "latitude")
TS_df <- unique(merge(TS_df, CyAN_merge[,c(1,2,9)], by = c("latitude", "longitude"), all.x = T))

my.sf.point <- st_as_sf(x = TS_df, 
                        coords = c("longitude", "latitude"),
                        crs = 4326 )

sf_trans <- st_transform(my.sf.point, proj4string(CyAN_1))
sf_trans <- as.data.frame(st_coordinates(sf_trans))
sf_trans <- as.data.frame(cbind(sf_trans, TS_df$Group))
names(sf_trans) <- c("longitude","latitude", "Group")
TS_df <- sf_trans

cc_tsm <- ggplot() +
  geom_sf(data = Flag_df, aes(fill = factor(Flag)), color = NA, show.legend = FALSE) + 
  scale_fill_manual(values = c("1" = "grey30", "3" = "#9b8d73", "4" = "grey80"), name= "Flag", labels = c("Too Low", "Land", "No data")) +
  geom_sf(data = FL_S3_tile, color = "black", 
          fill = NA, lwd = 0.5, alpha = 0.5) +
  geom_sf(data = FL, color = "black", fill = NA, lwd = 0.4)+
  ggnewscale::new_scale_fill() +
  geom_tile(data = Water,
            aes(x = long, y = lat, fill = (CyAN)), show.legend = FALSE) +
  scale_fill_viridis(name="Cyano DN",  limits=c(1, 253), values=c(0,1)) +
  ggnewscale::new_scale_fill() +
  geom_point(data = TS_df, 
             aes(x = longitude, y = latitude, fill = factor(Group), 
                 group = factor(Group)), 
             colour = "black",
             shape = 21,
             size = 2) +
  labs(fill = "Time Series\nLocations") +
  coord_sf(datum = sf::st_crs(4326)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, vjust = 1,
                                   hjust=1),
        axis.title = element_blank(),
        legend.direction = "vertical", 
        legend.box = "vertical")
cc_tsm
ggsave(filename="CyAN_CyFi_TimeSeries_Map.pdf", plot = cc_tsm, device = "pdf", path = Figure_filepath, width = 6, height = 5, dpi = 300, units = "in")
