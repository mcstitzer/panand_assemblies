
library(dplyr)

a=read.csv('AllAssembliesTelomeres_formatted_10TR_withLengths.csv')
a$LeftTelomere=as.logical(a$LeftTelomere)
a$RightTelomere=as.logical(a$RightTelomere)

a %>% group_by(Code,Contig) %>% dplyr::summarize(left=sum(LeftTelomere), right=sum(RightTelomere), both=sum(LeftTelomere&RightTelomere)) %>% group_by(Code) %>% summarize(l=sum(left), r=sum(right), b=sum(both), one=(sum(left)+sum(right))-(sum(both)*2)) %>% data.frame() %>% summarize(range(b, na.rm=T))



sort(table(a$Code[a$LeftTelomere & a$RightTelomere]))
