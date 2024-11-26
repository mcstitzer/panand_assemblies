library(data.table)
library(ggmap)

#library(devtools)
#remotes::install_github("8Ginette8/gbif.range")

library(gbif.range)
library(terra)
library(rnaturalearth)
library(tidyterra)
library(ggspatial)
library(ggnewscale)


ploidycolors=c( '#FFC857', '#A997DF', '#E5323B', '#2E4052', '#97cddf')
names(ploidycolors)=c('Diploid', 'Tetraploid', 'Hexaploid', 'Octaploid', 'Paleotetraploid')


asize=read.table('../general_summaries/panand_assembly_sizes.txt', header=T, sep='\t')

ll=fread('../general_summaries/panand_latlong.txt', fill=T)
asize$latloose=as.numeric(ll$latloose[match(asize$V2, ll$genome)])
asize$longloose=as.numeric(ll$longloose[match(asize$V2, ll$genome)])

## full andropogoneae distribution, from Taylor's paper
andro=read.csv('AuBuchonElder2022_Tropicos_Accepted_LatLon.csv')


bbox <- c(left = -170, bottom = -60, right = 170, top = 80)
world=ggmap(get_stadiamap(bbox, zoom =5, maptype="stamen_toner_lines"), extent = "device")

world + geom_point(data=andro, aes(x=Longitude, y=Latitude), color='gray', size=1, alpha=0.2) + geom_point(data=asize, aes(x=longloose, y=latloose, color=ploidy), size=2) + scale_color_manual(values=ploidycolors) + theme(legend.position='NULL')

pdf('~/transfer/panand_map_ploidy.pdf', 8,5)
world + geom_text(data=asize, aes(x=longloose, y=latloose, color=ploidy, label=V2)) + scale_color_manual(values=ploidycolors)
world + geom_point(data=asize, aes(x=longloose, y=latloose, color=ploidy, label=V2)) + scale_color_manual(values=ploidycolors)
dev.off()




# Plot species records
countries = vect(ne_countries(type = "countries",returnclass = "sf"))
plot(countries,col = "#bcbddc")
points(asize[,c("longloose","latloose")],pch=20,col="#99340470",cex=1.5)

# Download ecoregion and read
eco_terra = read_bioreg(bioreg_name = "eco_terra", save_dir = NULL)

andro_decimal=andro[,c('Longitude', 'Latitude')]
colnames(andro_decimal)=paste0('decimal', colnames(andro_decimal))
andro_decimal=rbind(andro_decimal, data.frame('decimalLongitude'=asize$longloose, 'decimalLatitude'=asize$latloose))
# Range
range_andro = get_range(occ_coord = andro_decimal[!is.na(andro_decimal$decimalLongitude),],
                        bioreg = eco_terra,
                        bioreg_name = "ECO_NAME")

range_andro_do20 = get_range(occ_coord = andro_decimal[!is.na(andro_decimal$decimalLongitude),],
                        bioreg = eco_terra,
                        bioreg_name = "ECO_NAME", degrees_outlier = 20)


range_andro_do100 = get_range(occ_coord = andro_decimal[!is.na(andro_decimal$decimalLongitude),],
                             bioreg = eco_terra,
                             bioreg_name = "ECO_NAME", degrees_outlier = 100)


plot(countries,col = "gray90", border='gray70', lwd=0.5)
plot(range_andro,col = "#238b45",add = TRUE,axes = FALSE,legend = FALSE)
points(asize[,c("longloose","latloose")],pch=21,bg=ploidycolors[asize$ploidy],col='black',cex=1, alpha=0.8)

ggplot() +geom_sf(data = ne_countries(type = "countries",returnclass = "sf"), fill='gray90', color='gray70', lwd=0.5) + 
  geom_spatraster(data = range_andro, alpha = 0.7)
  

# Crop the raster to the extent of the countries
range_andro_cropped <- crop(range_andro, countries)

# Mask the raster using the country boundaries
range_andro_masked <- mask(range_andro_cropped, countries)
names(range_andro_masked) <- "value"

# Plot using ggplot and tidyterra
ggplot() +
  geom_sf(data = ne_countries(type = "countries",returnclass = "sf"), fill = "gray90", color = "gray80", lwd = 0.5) +
  geom_spatraster(data = range_andro_masked, aes(fill = value), alpha = 0.9) +
  scale_fill_gradient(low = "white", high = "#238b45", na.value = NA) +
  new_scale_fill() +
#  geom_point(data=asize, aes(x=longloose, y=latloose), color='black', size=2.2) +
  geom_point(data=asize, aes(x=longloose, y=latloose, fill=ploidy), pch=21, color='black',size=2) + scale_fill_manual(values=ploidycolors) + theme(legend.position='NULL')+
  xlab('')+ylab('')


## fewer outliers
# Plot using ggplot and tidyterra
names(range_andro_do20)
ggplot() +
  geom_sf(data = ne_countries(type = "countries",returnclass = "sf"), fill = "gray90", color = "gray80", lwd = 0.5) +
  geom_spatraster(data = range_andro_do20, aes(fill = layer), alpha = 0.9) +
  scale_fill_gradient(low = "white", high = "#238b45", na.value = NA) +
  new_scale_fill() +
  #  geom_point(data=asize, aes(x=longloose, y=latloose), color='black', size=2.2) +
  geom_point(data=asize, aes(x=longloose, y=latloose, fill=ploidy), pch=21, color='black',size=2) + scale_fill_manual(values=ploidycolors) + theme(legend.position='NULL')+
  xlab('')+ylab('')

## see the text
names(range_andro_do20)
ggplot() +
  geom_sf(data = ne_countries(type = "countries",returnclass = "sf"), fill = "gray90", color = "gray80", lwd = 0.5) +
  geom_spatraster(data = range_andro_do20, aes(fill = layer), alpha = 0.9) +
  scale_fill_gradient(low = "white", high = "#238b45", na.value = NA) +
  new_scale_fill() +
  #  geom_point(data=asize, aes(x=longloose, y=latloose), color='black', size=2.2) +
  geom_point(data=asize, aes(x=longloose, y=latloose, fill=ploidy), pch=21, color='black',size=2) + scale_fill_manual(values=ploidycolors) + theme(legend.position='NULL')+
  xlab('')+ylab('')+ geom_text(data=asize, aes(x=longloose, y=latloose, label=V2, color=ploidy)) + scale_color_manual(values=ploidycolors)







## no outliers
rangemap=ggplot() +
  geom_sf(data = ne_countries(type = "countries",returnclass = "sf"), fill = "gray90", color = "gray80", lwd = 0.5) +
  geom_spatraster(data = range_andro_do100, aes(fill = layer), alpha = 0.9) +
  scale_fill_gradient(low = "white", high = "#238b45", na.value = NA) +
  new_scale_fill() +
  #  geom_point(data=asize, aes(x=longloose, y=latloose), color='black', size=2.2) +
  geom_point(data=asize, aes(x=longloose, y=latloose, fill=ploidy), pch=21, color='black',size=2) + scale_fill_manual(values=ploidycolors) + theme(legend.position='NULL')+
  xlab('')+ylab('')


ggplot() +
  geom_sf(data = ne_countries(type = "countries",returnclass = "sf"), fill = "gray90", color = "gray80", lwd = 0.5) +
  geom_spatraster(data = range_andro_do20, aes(fill = layer), alpha = 0.9) +
  scale_fill_gradient(low = "white", high = "#238b45", na.value = NA) +
  new_scale_fill() +
  #  geom_point(data=asize, aes(x=longloose, y=latloose), color='black', size=2.2) +
  geom_point(data=andro_decimal, aes(x=decimalLongitude, y=decimalLatitude), color='gold', size=0.5, alpha=0.9)+ theme(legend.position='NULL')+
  xlab('')+ylab('')

