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




## which occurrences go to each species>??

andro$genome=NA
andro$genome[andro$AcceptedSpeciesCorrected=='Andropogon chinensis']='achine'
andro$genome[andro$AcceptedSpeciesCorrected=='Andropogon gerardi']='agerar'
andro$genome[andro$AcceptedSpeciesCorrected=='Andropogon burmanicus']='atenui'
andro$genome[andro$AcceptedSpeciesCorrected=='Andropogon virginicus']='avirgi' ## had not been switched for taylor's paper
andro$genome[andro$AcceptedSpeciesCorrected=='Bothriochloa laguroides']='blagur'
andro$genome[andro$AcceptedSpeciesCorrected=='Chrysopogon serrulatus']='cserru'
andro$genome[andro$AcceptedSpeciesCorrected=='Cymbopogon citratus']='ccitra'
andro$genome[andro$AcceptedSpeciesCorrected=='Cymbopogon refractus']='crefra'
andro$genome[andro$AcceptedSpeciesCorrected=='Elionurus tripsacoides']='etrips'
andro$genome[andro$AcceptedSpeciesCorrected=='Hemarthria compressa']='hcompr'
andro$genome[andro$AcceptedSpeciesCorrected=='Heteropogon contortus']='hconto'
andro$genome[andro$AcceptedSpeciesCorrected=='Ischaemum rugosum']='irugos'
andro$genome[andro$AcceptedSpeciesCorrected=='Pogonatherum paniceum']='ppanic'
andro$genome[andro$AcceptedSpeciesCorrected=='Rhytachne rottboellioides']='rrottb'
andro$genome[andro$AcceptedSpeciesCorrected=='Rottboellia tuberculosa']='rtuber'
andro$genome[andro$AcceptedSpeciesCorrected=='Schizachyrium microstachyum']='smicro'
andro$genome[andro$AcceptedSpeciesCorrected=='Schizachyrium scoparium']='sscopa'
andro$genome[andro$AcceptedSpeciesCorrected=='Sorghastrum nutans']='snutan'
andro$genome[andro$AcceptedSpeciesCorrected=='Thelepogon elegans']='telega'
andro$genome[andro$AcceptedSpeciesCorrected=='Themeda triandra']='ttrian'
andro$genome[andro$AcceptedSpeciesCorrected=='Tripsacum dactyloides']='tdac'
andro$genome[andro$AcceptedSpeciesCorrected=='Urelytrum digitatum']='udigit'
andro$genome[andro$AcceptedSpeciesCorrected=='Vossia cuspidata']='vcuspi'
andro$genome[andro$AcceptedSpeciesCorrected=='Zea diploperennis']='zdip'
andro$genome[andro$AcceptedSpeciesCorrected=='Zea luxurians']='zluxur'
andro$genome[andro$AcceptedSpeciesCorrected=='Zea mays subsp. huehuetenangensis']='zmhuet'
andro$genome[andro$AcceptedSpeciesCorrected=='Zea mexicana']='zTILm'
andro$genome[andro$AcceptedSpeciesCorrected=='Zea mays subsp. parviglumis']='zTILp'
andro$genome[andro$AcceptedSpeciesCorrected=='Zea nicaraguensis']='znicar'
andro$genome[andro$AcceptedSpeciesCorrected=='Sorghum bicolor']='sbicol'
andro$genome[andro$AcceptedSpeciesCorrected=='Zea mays subsp. mays']='zmB735'



andro$genusMap=NA
andro$genusMap[andro$genus=='Andropogon']='andro'
andro$genusMap[andro$genus=='Andropogon']='andro'
andro$genusMap[andro$genus=='Andropogon']='andro'
andro$genusMap[andro$genus=='Andropogon']='andro' ## had not been switched for taylor's paper
andro$genusMap[andro$genus=='Bothriochloa']='blagur'
andro$genusMap[andro$genus=='Chrysopogon']='cserru'
andro$genusMap[andro$genus=='Cymbopogon']='cymbo'
andro$genusMap[andro$genus=='Cymbopogon']='cymbo'
andro$genusMap[andro$genus=='Elionurus']='etrips'
andro$genusMap[andro$genus=='Hemarthria']='hcompr'
andro$genusMap[andro$genus=='Heteropogon']='hconto'
andro$genusMap[andro$genus=='Ischaemum']='irugos'
andro$genusMap[andro$genus=='Pogonatherum']='ppanic'
andro$genusMap[andro$genus=='Rhytachne']='rrottb'
andro$genusMap[andro$genus=='Rottboellia']='rtuber'
andro$genusMap[andro$genus=='Schizachyrium']='schiza'
andro$genusMap[andro$genus=='Schizachyrium']='schiza'
andro$genusMap[andro$genus=='Sorghastrum']='snutan'
andro$genusMap[andro$genus=='Thelepogon']='telega'
andro$genusMap[andro$genus=='Themeda']='ttrian'
andro$genusMap[andro$genus=='Tripsacum']='tdac'
andro$genusMap[andro$genus=='Urelytrum']='udigit'
andro$genusMap[andro$genus=='Vossia']='vcuspi'
andro$genusMap[andro$genus=='Zea']='zea'
andro$genusMap[andro$genus=='Zea']='zea'
andro$genusMap[andro$genus=='Zea']='zea'
andro$genusMap[andro$genus=='Zea']='zea'
andro$genusMap[andro$genus=='Zea']='zea'
andro$genusMap[andro$genus=='Zea ']='zea'
andro$genusMap[andro$genus=='Sorghum']='sbicol'
andro$genusMap[andro$genus=='Zea']='zea'



## now, plot these occurrences for each species!!

for(i in asize$V2){
  pdf(paste0('maps_by_sp/', i, '_range.pdf'),8,4)
  species=i
  if(i%in%c('tdacn1', 'tdacs1')){
    species='tdac'
  }
  if(i%in%c('zTIL25', 'zTIL18')){
    species='zTILm'
  }
  if(i%in%c('zdmomo', 'zdgigi')){
    species='zdip'
  }
  if(i%in%c('zTIL01', 'zTIL11')){
    species='zTILp'
  }
  genus=i
  if(substr(i,1,1)=='z'){
    genus='zea'
  }
  if(i%in%c('tdacn1', 'tdacs1')){
    genus='tdac'
  }
  if(i%in%c('smicro', 'sscopa')){
    genus='schiza'
  }
  if(i%in%c('achine', 'agerar', 'atenui', 'avirgi')){
    genus='andro'
  }
  if(i%in%c('ccitra', 'crefra')){
    genus='cymbo'
  }
  genusdim=nrow(andro[andro$genusMap==genus,])
 # print(genusdim)
  print(ggplot() +
          geom_sf(data = ne_countries(type = "countries",returnclass = "sf"), fill = "gray90", color = "gray80", lwd = 0.5) +
          geom_spatraster(data = range_andro_do20, aes(fill = layer), alpha = 0.9) +
          scale_fill_gradient(low = "white", high = "#238b45", na.value = NA) +
          new_scale_fill() +
          #  geom_point(data=asize, aes(x=longloose, y=latloose), color='black', size=2.2) +
  #        geom_point(data=andro[andro$genusMap==genus,], aes(x=Longitude, y=Latitude), color='gray50',size=0.4, alpha=ifelse(genusdim>2000,0.1,0.3))+
          geom_point(data=andro[andro$genome==species,], aes(x=Longitude, y=Latitude), color='black',size=1)+
          geom_point(data=asize[asize$V2==i,], aes(x=longloose, y=latloose, fill=ploidy), pch=21, color='black',size=2) + scale_fill_manual(values=ploidycolors) + theme(legend.position='NULL')+
          xlab('')+ylab(''))
  dev.off()
}



