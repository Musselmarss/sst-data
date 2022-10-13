
setwd("D:/Environmental data for mapping/Predicted SST/RCP2.6")

nc<-"tos_Omon_MPI-ESM-MR_rcp26_r1i1p1_200601-210012.nc"
month<-seq(as.Date('2030-01-01'), as.Date('2033-12-31'), by="month") # generate sequence of month
mon<-strtrim(month, 7) #select first 7 characters


####################################################################
#### extract coordinates and sst in each month into a dataframe ####
####################################################################

library(ncdf4);library(ncdf4.helpers)


t<-nc_open(nc) # open nc file
lat<-ncvar_get(t,"lat") # extract latitude
lon<-ncvar_get(t,"lon") # extract longitude
time1<-as.character(nc.get.time.series(t, v="tos",time.dim.name="time")) # extract all time
time2<-which(strtrim(time1,7) %in% mon) # select desired time

coord<-as.data.frame(cbind(as.numeric(lon), as.numeric(lat))) # set up dataframe for coordinates
colnames(coord)[1:2]<- c("x", "y") # rename columns
 
list.sst<-list() # set up an empty list
count1<-1

# create dataframe for each month sst data
for (i in time2) {
  tos<-as.numeric(nc.get.var.subset.by.axes(t, "tos", axis.indices = list(T=i)))
  new<-cbind(coord, tos,i)
  colnames(new)[3:4]<-c("sst","time")
 list.sst[[count1]]<-new
 count1<-count1+1
}
nc_close(t)  #need to close nc file otherwise it keeps consuming RAM and slows down PC
names(list.sst)<-mon # chnage name of items in the list


####################################################
#### rasterize and crop dataframe of each month ####
####################################################

library(raster)
r<-raster("D:/Environmental data for mapping/Present SST/NOAA/v2.1/oisst-avhrr-v02r01.20160101.nc", varname="sst") # input noaa data to get the raster template resolution 0.25 x 0.25?

ex<-extent(110,180,-50,0) # set up extent

list.cropped<-list() # set up an empty list for storing cropped raster below
count2<-1

# turn the dataframes with sst data into rasters which is then cropped
for (g in 1:length(list.sst)) {
  tt<-rasterize(coord[,1:2], r, list.sst[[g]][["sst"]], fun=mean)
  cropped_tt<-crop(tt,ex)
  projection(cropped_tt)<- "+proj=longlat +datum=WGS84"
  list.cropped[[count2]]<-cropped_tt
  count2<-count2+1
}
names(list.cropped)<-paste("cropped.", mon) # rename items in this list

##########################################
#### Interpolate the regridded raster ####
##########################################

# function for interpolating value of empty grids
interp<-function(x){
  return(mean(x, na.rm=TRUE))
}

list.inter<-list() # set up list for storing interpolated rasters
count3<-1

# interpolate empty rasters and add projection to rasters
for (u in 1:length(list.cropped)) {
  r.inter<- focal(list.cropped[[u]], matrix(1, ncol=5, nrow=5), fun=interp, pad=T, NAonly=T) 
  projection(r.inter)<- "+proj=longlat +datum=WGS84"
  list.inter[[count3]]<-r.inter
  count3<-count3+1
}
names(list.inter)<-paste("inter.", mon) # rename the list items

land <- r <- crop(r, ex)
values(land)[which(!is.na(r[]))] <- 1
values(land)[which(is.na(r[]))] <- 0

for (u in 1:length(list.inter)) {
  values(list.inter[[u]])[which(values(land) == 0)] <- NA
}

plot(list.inter[[1]])
################################
#### Select coastal 3 grids ####
################################

# function for extract coastal 3 grids
coast<-function(y){
  if(any(is.na(y))) return(y[25]) else return(NA)
}

list.surround<-list()
count4<-1

# extract coastal 3 grids
for (v in 1:length(list.inter)) {
  r.surround <- focal(list.inter[[v]], matrix(1, ncol=7, nrow=7), fun=coast)
  list.surround[[count4]]<-r.surround
  count4<-count4+1
}
names(list.surround)<-paste("surround.", mon) # rename list items

a<-data.frame(rasterToPoints(list.surround[[1]]))
##############################################
#### extract grids from highlighted areas ####
##############################################
 
# specific coordinates to create polygons 
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

pols <- spPolygons(pa,NZ, crs = ("+proj=longlat +datum=WGS84")) # create polygons
points <- rasterize(pols, list.surround[[1]]) # rasterize polygons

co <- data.frame(rasterToPoints(points)) # create dataframe with extracted coordinates of areas masked by polygons
colnames(co)[1:2]<- c("Lon", "Lat")
co <- co[,-3]

spts <- SpatialPoints(co, proj4string=CRS("+proj=longlat +datum=WGS84")) # make points as coordinates

count5<-1

# extract sst from areas masked by polygons
for (s in 1:length(list.surround)) {
  ex <- extract(list.surround[[s]], spts, fun='mean', na.rm=TRUE, df=TRUE, weights = TRUE)
temp<-cbind(co,ex)
temp<-na.omit(temp[,c(3,1,2,4)])
colnames(temp)<- c("ID","Lon", "Lat","SST")
# write the extracted vaues in .csv file for further processing
write.csv(temp,paste0("RCP26_",mon[count5],".csv"), row.names = FALSE) 
count5<-count5+1
}


##################################
#### tidy up sst data by grid ####
##################################
setwd("D:/")
library(plyr)
list_file<-list.files(pattern = "RCP85_") # list all files with specified words
all_data<-ldply(list_file, read.csv) # stack data from all files with specified words

date<-data.frame(mon)
colnames(date)<- "Date"

id<-read.csv(choose.files()) # choose file containing id for the coordinates

for (g in id$ID) {
  a<- subset(all_data, ID==g)
  a<- cbind(a, date)
  a<- a[,c(5,1,2,3,4)]
  write.csv(a, paste("D:/", g,"present.csv"), row.names = FALSE)
}

##########################################
#### To copy monthly files into daily ####
##########################################
setwd("D:/Environmental data for mapping/Predicted SST/RCP2.6/RCP_26")
file<-list.files(pattern= "RCP")

library(stringr)
for (o in file) {
  file.rename(o, paste0(str_remove(o, ".csv"), "-01", ".csv"))
}

file<-list.files(pattern= "RCP")

library(sjmisc);library(lubridate)
for (y in file) {
  a<-read.csv(y)
  if(str_contains(y, pattern= c("-01-","-03-","-05-","-07-","-08-","-10-","-12-"), logic = "or")) {
    for (t in sprintf("%02d", 2:31)) {
      write.csv(a, paste0(str_replace(y,"01.csv", ""),t,".csv"), row.names = F) }
  }
  else if(str_contains(y, pattern=c("-04-","-06-","-09-","-11-"), logic = "or")) {
    for (h in sprintf("%02d", 2:30)) {
      write.csv(a, paste0(str_replace(y,"01.csv", ""),h,".csv"), row.names = F)}
  }
  else {
    if(leap_year(as.numeric(str_remove(strtrim(y,10),"RCP26_")))){
      for (m in sprintf("%02d", 2:29)) {
        write.csv(a, paste0(str_replace(y,"01.csv", ""),m,".csv"), row.names = F)}
    }
    else {
      for (s in sprintf("%02d", 2:28)) {
        write.csv(a, paste0(str_replace(y,"01.csv", ""),s,".csv"), row.names = F)}
    }
  }
}


#Finding centroids of each country
library(sf);library(rnaturalearth)

world <- ne_countries(scale = "medium", returnclass = "sf")
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
ggplot(data = a) +
  geom_tile(aes(x=x, y=y, fill=layer))+ 
  scale_fill_viridis_c() + borders(colour = "black", fill = "grey") + 
  coord_sf(ylim = c(-49,0), xlim = c(110, 180), expand = FALSE)+
  labs(y="Latitude", x = "Longitude")+
  theme_classic() + theme(text = element_text(size=20),
                          axis.text.x = element_text(size = 20, margin = margin(t=15), colour = "black"),
                          axis.text.y = element_text(size = 20, margin = margin(r=15), colour = "black"),
                          axis.title.x = element_text(size = 20 , margin = margin(t = 15)),
                          axis.title.y = element_text(size = 20 , margin = margin(r = 15)),
                          plot.title = element_text(size = 20))



