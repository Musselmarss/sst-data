
setwd("C:\\Users\\Martin Cheng\\Desktop\\Martin\\sst data")


###########################################################
#################### Download SST data ####################
###########################################################

# generate url of the temporal temperature data
# choose period
dt = seq(as.Date('2022-01-01'), as.Date('2022-01-16'), by="days")

# convert date to new format
dt = format(dt, '%Y%m%d')


for (i in dt) {
 dtt = strtrim(i, 6) 
 # assemble url to query NOAA database
url_base = paste0("https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/", dtt,"/")
data_file = paste0("oisst-avhrr-v02r01.",i,".nc")
# define data url
data_url = paste0(url_base, data_file)
#download netcdf
suppressPackageStartupMessages(library(ocedata))
if(!file.exists(data_file)){
  download.file(url = data_url, destfile = data_file, mode = "wb")
} else {
  message('SST data already downloaded! Located at:\n', data_file)
} 
}


###########################################################
##################### Extract SST data ####################
###########################################################

#specific coordinates to create polygons
pa<-rbind(c(-8.33056241071365, 137.53527841519727),
          c(-6.674997121768877, 141.24865732144727),
          c(-5.932449052213845, 144.80822763394727),
          c(-8.55877378817969, 147.68664560269727),
          c(-10.865159283071419, 150.87268075894727),
          c(-12.280452468217703, 149.25219736050977),
          c(-11.802319380707218, 145.69812021207227),
          c(-22.58816876983026, 154.33886728238477),
          c(-34.71409231104185, 159.81555185269727),
          c(-38.10821380556104, 149.88391122769727),
          c(-11.002884109462219, 142.16307351001464),
          c(-10.253692762713492, 139.69536736885595))
pa <- pa[,2:1]

NZ <- rbind(c(-33.92276529900936, 172.38989453967466),
            c(-35.89564193945088, 173.97192578967466),
            c(-37.58586079433415, 175.40014844592466),
            c(-40.75342115716705, 175.55395703967466),
            c(-41.36642567586272, 174.35644727404966),
            c(-41.84288994530583, 173.15893750842466),
            c(-45.28061213067501, 168.41284375842466),
            c(-46.03313111702218, 166.69897657092466),
            c(-46.10935022479369, 164.72143750842466))
NZ <- NZ[,2:1]

#create and rasterize polygons
library(raster);library(ncdf4)
pols <- spPolygons(pa,NZ, crs = ("+proj=longlat +datum=WGS84"))
ra<-raster("avhrr-only-v2.20200101.nc", varname="sst")
points <- rasterize(pols, ra)

#Provide longitude & latitude of specific area
co <- data.frame(rasterToPoints(points))
colnames(co)[1:2]<- c("Lon", "Lat")
co <- co[,-3]

#Make points as coordinates
spts <- SpatialPoints(co, proj4string=CRS("+proj=longlat +datum=WGS84"))

#Rasterize the map
library(stringr)

#Set function for coastal pixel extraction
coast<-function(x){
    if(any(is.na(x))) return(x[25]) else return(NA)
}

#Extract and output selected sst data
for (i in dir(pattern = "oisst")) {
  r <- raster(i, varname="sst")
  e<-extent(110,180,-49,-1)
  r_cropped<-crop(r,e)
  
#Select three pixels near the coast
r_surround <- focal(r_cropped,
      matrix(1, ncol=7, nrow=7), fun=coast) 

#Extract raster values
  ex <- extract(r_surround, spts, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
  temp<-cbind(co,ex)
#temp<-temp[,c(3,1,2,4)]
  temp<-na.omit(temp[,c(3,1,2,4)])
  colnames(temp)<- c("ID","Lon", "Lat","SST")
#Write the extracted values in .csv file for further processing
  write.csv(temp,paste0("C:\\Users\\Martin Cheng\\Desktop\\Martin\\sst data\\Extracted\\",str_remove(i,".nc"),".csv"), row.names = FALSE) 
}

#Plot those points to make sure the locations are selected correctly
#plot(spts, add = TRUE)


###########################################################
########### Plot SST data check if the grids ##############
################ extracted are correct ####################
###########################################################

#Finding centroids of each country
library(sf);library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")
sf::sf_use_s2(FALSE)
#class(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))

#Correction of centroid of Fiji
Fi<-subset(world_points, world_points$name=="Fiji")
Fi$X <-178.318
world_points[71,] <-Fi

#Mapping temperature on the map with colour gradient
library(ggplot2)
library(viridis)
library(viridisLite)
sst_plot<-ggplot(data = temp) +
  geom_tile(aes(x=Lon, y=Lat, fill=SST))+ 
  scale_fill_viridis_c(option = "magma") + borders(colour = "black", fill = "grey", alpha= 1) + 
  geom_text(data=world_points,aes(x=X, y=Y, label=name,family="Comic Sans MS"), size = 5,
            color = "navy", fontface = "bold", check_overlap = FALSE) +
  coord_sf(ylim = c(-48,-1), xlim = c(110, 180), expand = FALSE)+
  labs(y="Latitude", x = "Longitude")+
  theme_classic() + theme(text = element_text(family="Comic Sans MS", size=20),
                          axis.text.x = element_text(size = 20, margin = margin(t=15), colour = "black"),
                          axis.text.y = element_text(size = 20, margin = margin(r=15), colour = "black"),
                          axis.title.x = element_text(size = 20 , margin = margin(t = 15)),
                          axis.title.y = element_text(size = 20 , margin = margin(r = 15)),
                          plot.title = element_text(size = 20))
#c(-48.914227,-1.01865), xlim = c(120.769626, 180.143358)

#Plot the selected map with chosen pixels
library(ggspatial)
sst_plot<-ggplot(data = temp) +
  geom_tile(aes(x=Lon, y=Lat, fill=SST), fill="lightgreen")+ 
  borders(colour = "black", fill = "grey95", alpha= 1) + 
  geom_text(data=world_points,aes(x=X, y=Y, label=name), size = 15,
            color = "black", fontface = "bold", check_overlap = FALSE) +
  annotation_scale(location = "bl", text_cex = 2.5, height = unit(0.5,"cm"), pad_x = unit(1.5,"cm"), pad_y = unit(2,"cm")) +
  annotation_north_arrow(location = "bl", 
                         which_north = "true",
                         height = unit(5.25, "cm"),
                         width = unit(5.25, "cm"),
                         pad_x = unit(2.3, "cm"), pad_y = unit(3.5,"cm"),
                         style = north_arrow_fancy_orienteering(text_size = 35)) +
  coord_sf(ylim = c(-48,-1), xlim = c(110, 181), expand = FALSE)+
  labs(y="Latitude", x = "Longitude")+
  theme_bw(base_size = 15)+ theme(axis.title.x = element_text(size = 45, margin = margin(t = 25)),
                                  axis.text.x = element_text(size = 42, margin = margin(t=15), colour = "black"),
                                  axis.title.y = element_text(size = 45, margin = margin(r = 25)),
                                  axis.text.y = element_text(size = 42,margin = margin(r=15), colour = "black"),
                                  plot.title = element_text(size = 45))

#Output the mapped plot
library(ragg)
agg_png(filename = "map2.png", width = 20, height = 20, units = "cm", res = 500)
plot(sst_plot)
invisible(dev.off())

###########################################################
############### Tidy up sst data by grid ##################
###########################################################

setwd("C:/Users/Martin Cheng/Desktop/Martin/sst data/Extracted")

library(plyr)
id<-read.csv("Extracted/coordinate_ID.csv")

#Take very long time
#for (y in id$ID) {
 # rr<-data.frame()
  #for (g in dir(pattern = "avhrr-only-v2.")) {
   # b<-read.csv(g)
    #a<-data.frame(subset(b, ID==y))
    #a[,5]<-str_remove(str_remove(g, "avhrr-only-v2."), ".csv")
    #colnames(a)[5]<-"Date"
    #a<-a[,c(5,1,2,3,4)]
    #rr<-rbind(rr,a)
    #write.csv(rr, paste("D:/Present SST/extracted SST/Tidied/", y,"present.csv"), row.names = FALSE)
  #}
#}

list_file<-list.files(pattern = "oisst")
all_data<-ldply(list_file, read.csv)

dt<- seq(as.Date("2020-01-01"), as.Date("2022-01-01"), "days")
dtt<- format(dt, "%Y%m%d")

date<-data.frame(dtt)
colnames(date)<- "Date"

for (g in id$ID) {
  a<- subset(all_data, ID==g)
  a<- cbind(a, date)
  a<- a[,c(5,1,2,3,4)]
  write.csv(a, paste("C:\\Users\\Martin Cheng\\Desktop\\Martin\\sst data\\Extracted\\Tidied\\", g,"present.csv"), row.names = FALSE)
}
