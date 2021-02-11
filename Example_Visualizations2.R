library(tidyverse)
library(lubridate)
library(sp)
library(sf)
library(ggmap)
library(stringr)
library(scales)
library(ggspatial)
library(survival)
library(ggfortify)
library(PerformanceAnalytics)

######################################
####### STIP Acoustic Telemetry ######
######################################

receivers <- read_csv("C:/Users/tspaulding/Downloads/stip_20200630_acoustictelemetry_deployment.csv")
detections <- read_csv("C:/Users/tspaulding/Downloads/stip_20200630_acoustictelemetry_detections.csv")
fish <- read_csv("C:/Users/tspaulding/Downloads/stip_20200630_organismsampling_fish.csv")


reciever_list <- receivers %>%
  arrange(RKM) %>%
  distinct(STATION_ID)

receivers <- receivers %>%
  mutate(STATION_ID = factor(STATION_ID, levels = reciever_list$STATION_ID))

# Receiver Map
# set the limits of the map (based on guess and check with final map)
xmin <- -122.05
xmax <- -121.25
ymin <- 37.6
ymax <- 38.25

# Read in SHP file of streams
waterbodies <- st_read("H:/GROUPS/Fisheries/GIS/GIS_BASELAYERS/DMVS_Locations_Poly.shp")

waterbodies_proj <- st_transform(waterbodies, crs = 4326)

#Check the Map
ggplot()+
  geom_sf(data = waterbodies_proj, fill = "skyblue", color = "skyblue4")+
  geom_point(data = receivers, aes(x = DEPLOY_LONG, y = DEPLOY_LAT), shape = 21, fill = "red", color = "black", size = 2)+
  coord_sf(crs = st_crs(4326))+
  lims(x = c(xmin,xmax), y = c(ymin, ymax))+
  labs(x = "Longitude", y = "Latitude")

#Receiver Detection Data GANTT
ggplot(detections)+
  geom_point(aes(x = DETECT_DATETIME, y = STATION_ID))


#Fish Length-Weight
ggplot(fish)+
  geom_point(aes(x = LENGTH, y = WEIGHT))+
  geom_smooth(aes(x = LENGTH, y = WEIGHT))+
  scale_x_log10()+
  scale_y_log10()+
  labs(x = "Length (mm)", y = "Weight (g)")

#Survival Waterfall Plot
det_surv <- detections %>%
  left_join(receivers, by = "STATION_ID") %>%
  distinct(TRANS_ID, RKM) %>%
  arrange(TRANS_ID, desc(RKM))


km <- with(det_surv, Surv(RKM))
head(km,1000)

km_fit <- survfit(Surv(RKM)~1, data = det_surv)
summary(km_fit, times = c(seq(1, 1550, 100)))

autoplot(km_fit) 


###########################################
####### ADCP DATA #########################
###########################################

ew <- read_csv("C:/Users/tspaulding/Downloads/ADCP/ADCP/Turner_Cut/TurnerCut_2014/TurnerCut/2014-SL/1.7_SJTC_Database_SLADCP_QAQC/E-WVelocities.csv", 
               col_names = FALSE,
               na = c("","NaN"))
ns <- read_csv("C:/Users/tspaulding/Downloads/ADCP/ADCP/Turner_Cut/TurnerCut_2014/TurnerCut/2014-SL/1.7_SJTC_Database_SLADCP_QAQC/N-SVelocities.csv", 
               col_names = FALSE,
               na = c("","NaN"))
pos <- read_csv("C:/Users/tspaulding/Downloads/ADCP/ADCP/Turner_Cut/TurnerCut_2014/TurnerCut/2014-SL/1.7_SJTC_Database_SLADCP_QAQC/Position.csv", 
                col_names = FALSE,
                na = c("","NaN"))
time <- read_csv("C:/Users/tspaulding/Downloads/ADCP/ADCP/Turner_Cut/TurnerCut_2014/TurnerCut/2014-SL/1.7_SJTC_Database_SLADCP_QAQC/Timestamp.csv", 
                 col_names = FALSE,
                 na = c("","NaN"))

comb = NULL
for(m in 1:ncol(ns)){
  for(n in 1:nrow(ns)){
    
    DateTime <- time[n,1,drop = TRUE]
    Easting <- pos[1,m,drop = TRUE]
    Northing <- pos[2,m,drop = TRUE]
    ns_veloc <- ns[n,m,drop = TRUE]
    ew_veloc <- ew[n,m,drop = TRUE]
    
    comb = rbind(comb, data.frame(DateTime,Easting,Northing,ns_veloc,ew_veloc))
    
  }
}

comb <- comb %>%
  mutate(vec_mag = sqrt((ns_veloc^2)+(ew_veloc^2)),
         vec_dir = atan(ns_veloc/ew_veloc)*180/pi,
         ns_vec_end = Northing + (ns_veloc*200),
         ew_vec_end = Easting + (ew_veloc*200))


ggplot(data = comb %>% filter(DateTime == comb$DateTime[1]))+
  geom_point(aes(x = Easting, y = Northing, color = vec_mag))+
  geom_segment(aes(x = Easting, xend = ew_vec_end, y = Northing, yend = ns_vec_end))

ew <- read_csv("C:/Users/tspaulding/Downloads/ADCP/ADCP/Cornell_Nov2017_Iselton_Sutter/November 2017 with Cornell Work (GSB-03)/November 2017 with Cornell Work (GSB-03)/Task1.1_Isleton/SLADCPs/QAQC_Data/SLADCP_towardBank_E-W_Velocities.csv", 
               col_names = FALSE,
               na = c("","NaN"))
ns <- read_csv("C:/Users/tspaulding/Downloads/ADCP/ADCP/Cornell_Nov2017_Iselton_Sutter/November 2017 with Cornell Work (GSB-03)/November 2017 with Cornell Work (GSB-03)/Task1.1_Isleton/SLADCPs/QAQC_Data/SLADCP_towardBank_N-S_Velocities.csv", 
               col_names = FALSE,
               na = c("","NaN"))
pos <- read_csv("C:/Users/tspaulding/Downloads/ADCP/ADCP/Cornell_Nov2017_Iselton_Sutter/November 2017 with Cornell Work (GSB-03)/November 2017 with Cornell Work (GSB-03)/Task1.1_Isleton/SLADCPs/QAQC_Data/SLADCP_towardBank_Horizontal_Position.csv", 
                col_names = FALSE,
                na = c("","NaN"))
time <- read_csv("C:/Users/tspaulding/Downloads/ADCP/ADCP/Cornell_Nov2017_Iselton_Sutter/November 2017 with Cornell Work (GSB-03)/November 2017 with Cornell Work (GSB-03)/Task1.1_Isleton/SLADCPs/QAQC_Data/SLADCP_towardBank_Timestamp.csv", 
                 col_names = FALSE,
                 na = c("","NaN")) %>%
  mutate(DateTime = parse_datetime(X1, format = "%d-%b-%Y %H:%M:%S"))

comb = NULL
for(m in 1:ncol(ns)){
  for(n in 1:nrow(ns)){
    
    DateTime <- time[m,2,drop = TRUE]
    Easting <- pos[n,1,drop = TRUE]
    Northing <- pos[n,2,drop = TRUE]
    ns_veloc <- ns[n,m,drop = TRUE]
    ew_veloc <- ew[n,m,drop = TRUE]
    
    comb = rbind(comb, data.frame(DateTime,Easting,Northing,ns_veloc,ew_veloc))
    
  }
}

comb <- comb %>%
  mutate(vec_mag = sqrt((ns_veloc^2)+(ew_veloc^2)),
         vec_dir = atan(ns_veloc/ew_veloc)*180/pi,
         ns_vec_end = Northing + (ns_veloc),
         ew_vec_end = Easting + (ew_veloc))


ggplot(data = comb)+
  geom_point(aes(x = Easting, y = Northing, color = vec_mag))+
  geom_segment(aes(x = Easting, xend = ew_vec_end, y = Northing, yend = ns_vec_end, color = vec_mag))


#####################
#### TIME SERIES ####
#####################

env <- read_csv(file = "CCF_Environmental.csv")

ggplot(env)+
  geom_histogram(aes(x = `EL COND`))+
  labs(y = "Count", x = "Electrical Conductivity")

ggplot(env)+
  geom_point(aes(x = DATETIME_UTCminus8, y = `EL COND`))+
  labs(y = "Electrical Conductivity", x = "Date")

##########################
#### Observation Date ####
##########################

avian <- read_csv(file = "C:/Users/tspaulding/Downloads/avian_speciesactivity.csv") %>%
  pivot_longer(cols = SpeciesCount:FlyOver, values_to = "Count", names_to = "Behavior") %>%
  mutate(AvianSurveyDate = mdy(`AvianSurveyDate`))

ggplot(avian %>% filter(AvianLocation %in% c("Radial Gates NE", "Radial Gates SW", "Trash Racks")))+
  geom_bar(aes(x = AvianSurveyDate, y = Count, fill = SpeciesName), stat = "identity")+
  theme(
    axis.text = element_text(angle = 45, hjust = 1, vjust =-0.5)
  )

data = avian2 %>% select(WindSpeed_mph, TempF, Humidity_Perc, CloudCover_Perc, Precip, Feeding) %>%
  mutate(`WindSpeed_mph` = as.numeric(`WindSpeed_mph`),
         `TempF` = as.numeric(`TempF`),
         `Humidity_Perc` = as.numeric(`Humidity_Perc`),
         `CloudCover_Perc` = as.numeric(`CloudCover_Perc`),
         `Precip` = as.numeric(`Precip`))


chart.Correlation(data, histogram=TRUE, pch=19)
