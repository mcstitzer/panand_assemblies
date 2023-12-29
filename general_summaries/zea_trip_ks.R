




ksd$pgene=substr(ksd$V4,8,21)
ksp$pgene=substr(ksp$V4,8,21)

ksp$zea=substr(ksp$genome,1,1)=='z'
ksp$trip=substr(ksp$genome,1,4)=='tdac'

ksp$tripdup=ksp$pgene %in% ksd$pgene[ksd$trip & ksd$npgene==1]
ksp$zeadup=ksp$pgene %in% ksd$pgene[ksd$zea & ksd$npgene==1]


ksp %>% group_by(genome) %>% filter(zea|trip) %>%summarize(ks=median(V17, na.rm=T), dnds=median(V19,na.rm=T), n=n()) %>% data.frame
summary(ksp$V17)
ksp %>% group_by(genome) %>% filter(zea|trip) %>%summarize(ks=median(V17, na.rm=T), dnds=median(V19,na.rm=T), n=n()) %>% arrange(-ks)%>% data.frame
ksp %>% group_by(genome) %>% filter((zea|trip)|(zeadup|tripdup)) %>%summarize(ks=median(V17, na.rm=T), dnds=median(V19,na.rm=T), n=n()) %>% arrange(-ks)%>% data.frame
ksp %>% group_by(genome) %>% filter((zea|trip)&(zeadup|tripdup)) %>%summarize(ks=median(V17, na.rm=T), dnds=median(V19,na.rm=T), n=n()) %>% arrange(-ks)%>% data.frame
head(ksp)
ksp %>% group_by(genome,intraspecificdup) %>% filter((zea|trip)) %>%summarize(ks=median(V17, na.rm=T), dnds=median(V19,na.rm=T), n=n()) %>% arrange(-ks)%>% data.frame
ksp %>% group_by(genome,intraspecificdup) %>% filter((zea|trip)&!intraspecificdup) %>%summarize(ks=median(V17, na.rm=T), dnds=median(V19,na.rm=T), n=n()) %>% arrange(-ks)%>% data.frame


ksd$tripdup=ksd$pgene %in% ksd$pgene[ksd$trip & ksd$npgene==1]
ksd$zeadup=ksd$pgene %in% ksd$pgene[ksd$zea & ksd$npgene==1]
head(ksd) %>% data.frame
ksd %>% group_by(genome) %>% filter(trip|zea) %>% group_by(genome, tripdup, zeadup) %>% summarize(ks=mean(V17,na.rm=T), dnds=mean(V19,na.rm=T))
ksd %>% group_by(genome) %>% filter(trip|zea) %>% group_by(genome, tripdup, zeadup) %>% summarize(ks=mean(V17,na.rm=T), dnds=mean(V19,na.rm=T)) %>% arrange(-ks)
ksd %>% group_by(genome) %>% filter(trip|zea) %>% group_by(genome, tripdup, zeadup) %>% summarize(ks=mean(V17,na.rm=T), dnds=mean(V19,na.rm=T)) %>% filter(!genome %in% c('tdactm', 'tzopot')) %>% arrange(-ks)
ksd %>% group_by(genome) %>% filter(trip|zea) %>% group_by(genome, tripdup, zeadup) %>% summarize(ks=mean(V17,na.rm=T), dnds=mean(V19,na.rm=T)) %>% filter(!genome %in% c('tdactm', 'tzopot')) %>% arrange(-ks) %>% data.frame
ksd %>% group_by(genome) %>% filter(trip|zea) %>% group_by(genome, tripdup, zeadup) %>% summarize(ks=mean(V17,na.rm=T), dnds=mean(V19,na.rm=T), n=n()) %>% filter(!genome %in% c('tdactm', 'tzopot')) %>% arrange(-ks) %>% data.frame
ksd %>% group_by(genome) %>% filter(trip|zea) %>% group_by(genome, tripdup, zeadup) %>% summarize(ks=mean(V17,na.rm=T), dnds=mean(V19,na.rm=T), n=n()) %>% filter(!genome %in% c('tdactm', 'tzopot') & (tripdup|zeadup)) %>% arrange(-ks) %>% data.frame
