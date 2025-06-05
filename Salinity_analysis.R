####Salinity analysis####
new_lib <- "~/new_R_library"
.libPaths(c(new_lib, .libPaths()))

#analysis for salinity data
#code for sample map

library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(celestial)
library(ggspatial)
library(ggrepel)
library(scales)
library(raster)
library(terra)
library(dplyr)
library(maps)
library(stars)
library(paletteer)
library(reshape2)

world <- ne_countries(scale = "medium", returnclass = "sf")

setwd("~/Documents/Experiments/Salinity")

#download data from E.U. Copernicus Marine Service Information (CMEMS)
#Global Ocean Physics Reanalysis
#https://doi.org/10.48670/moi-00021
raster_path <- "cmems_mod_glo_phy_my_0.083deg_P1M-m_1740971824938.nc"
raster_path_SW08 <- "SW08_cmems_mod_glo_phy_my_0.083deg_P1D-m_1742300511722.nc" #for daily means in SW08 area
raster_path_WH <- "WH_cmems_mod_glo_phy_my_0.083deg_P1D-m_1742300158566.nc" #for daily means in WH area

r_brick <- brick(raster_path, varname="so")
r_brick_SW08 <- brick(raster_path_SW08, varname="so")
r_brick_WH <- brick(raster_path_WH, varname="so")

#### salinity map ####
# Convert the entire raster brick to a dataframe
r_df <- raster::as.data.frame(r_brick, xy = TRUE)

# Reshape the data to a long format
r_df_long <- r_df %>%
  pivot_longer(cols = -c(x, y), names_to = "time", values_to = "salinity")

r_df_long$time <- as.Date(gsub("X", "", r_df_long$time), format="%Y.%m.%d")

r_df_long <- r_df_long %>%
  dplyr::group_by(x, y) %>%
  mutate(mean_salinity = mean(salinity))

#create data frame with sampling locations
coords<-data.frame(site = c("SW08", "WH"),
                    lon=c(11.617000,8.149000),
                    lat=c(54.413000,53.513000))

#plot sampling coordinates on mean salinity map

transparent_palette <- alpha(paletteer::paletteer_c("grDevices::Blue-Red", n = 256), 0.55)

#plot map
map <- ggplot() +
  geom_raster(data = r_df_long, aes(x = x, y = y, fill = mean_salinity)) +
  scale_fill_gradientn(
    colors = transparent_palette,
    name = "Salinity",
    limits = range(r_df_long$mean_salinity, na.rm = TRUE)
  ) +
  geom_sf(data = world, fill = "white", color = "black") +
  theme_light(base_size = 17) +
  annotation_scale(location = "bl") +
  coord_sf(xlim = c(5,19), ylim = c(52,59), expand = FALSE)+
  geom_point(data = coords, aes(x = lon, y = lat, color = site) , size = 3, 
             shape = 16)+
  ylab("Latitude")+ xlab("Logitude")


#### stats sampling locations ####
#only look at sampling locations
#r_df_WH <- raster::as.data.frame(r_brick_WH, xy = TRUE)#to use all data from WH area
#r_df_SW08 <- raster::as.data.frame(r_brick_SW08, xy = TRUE)#to use all data from SW08 area

coord<-data.frame(lon=c(11.617000,8.149000),
                   lat=c(54.413000,53.513000))#specific locations

SW08_sal <- terra::extract(x = r_brick_SW08, y = coord[1,])
WH_sal <- terra::extract(x = r_brick_WH, y = coord[2,])

#subset the extracted data to get only summer months
layer_dates_WH <- as.Date(gsub("X", "", names(r_brick_WH)), format="%Y.%m.%d")
layer_dates_SW08 <- as.Date(gsub("X", "", names(r_brick_SW08)), format="%Y.%m.%d")
summer_indices_WH <- which(months(layer_dates_WH) %in% month.name[6:10])
summer_indices_SW08 <- which(months(layer_dates_SW08) %in% month.name[6:10])

WH_summer <- WH_sal[, summer_indices_WH]
SW08_summer <- SW08_sal[, summer_indices_SW08]

#transform for plotting
SW08_df <- data.frame(salinity = as.numeric(t(SW08_summer)))
SW08_df$site <- "SW08"
SW08_df$year_month <- gsub("X","", names(SW08_summer))

WH_df <- data.frame(salinity = as.numeric(t(WH_summer)))
WH_df$site <- "WH"
WH_df$year_month <- gsub("X","", names(WH_summer))

#merge dfs
combined_df <- rbind(WH_df, SW08_df)
combined_df$site <- factor(combined_df$site, levels = c("WH", "SW08"))
combined_long <- melt(combined_df, id.vars = c("site", "year_month"), 
                      value.name = "salinity")

#violin plot of daily SST per sampling locations (1993 - 2021)
ggplot(combined_long, aes(x = site, y = salinity, fill = site, col = site)) +
  geom_violin(alpha = 0.2) +
  scale_fill_manual(values = c("#DE3C22", "#83B3C2"))+
  scale_color_manual(values = c("#DE3C22", "#83B3C2"))+
  theme_light(base_size = 14)+
  xlab("Sampling location")+
  ylab("Salinity [PSU]")+
  stat_summary(fun = "mean",
               geom = "point",
               color = c("#DE3C22", "#83B3C2"))

#mean salinity per sampling location
mean_sal <- combined_df %>%
  dplyr::group_by(site) %>%
  summarise(Mean = mean(salinity),
            SD = sd(salinity),
            Max = max(salinity),
            Min = min(salinity))

hist(combined_df$salinity, breaks = 50)
wilcox.test(salinity ~ site, data = combined_df)


