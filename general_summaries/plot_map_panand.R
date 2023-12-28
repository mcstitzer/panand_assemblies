library(data.table)
library(ggmap)


asize=read.table('../general_summaries/panand_assembly_sizes.txt', header=F)

ll=fread('../panand_latlong.txt', fill=T)
colnames(ll)=c('genome', 'latitude', 'longitude', 'latloose', 'longloose', 'nope', 'notreal')
asize$latloose=ll$latloose[match(asize$V2, ll$genome)]
asize$longloose=ll$longloose[match(asize$V2, ll$genome)]


bbox <- c(left = -170, bottom = -60, right = 170, top = 80)
world=ggmap(get_stadiamap(bbox, zoom =3, maptype="stamen_toner_lines"), extent = "device")

pdf('~/transfer/panand_map_ploidy.pdf', 8,5)
world + geom_text(data=asize, aes(x=longloose, y=latloose, color=ploidy, label=V2)) + scale_color_manual(values=ploidycolors)
world + geom_point(data=asize, aes(x=longloose, y=latloose, color=ploidy, label=V2)) + scale_color_manual(values=ploidycolors)
dev.off()
