library(rtracklayer)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(dplyr)
library(stringr)
library(tidyverse)
library(RColorBrewer)
library(viridis)

genomebase='smicroHap1_aggressivecorrection'
genomesize=1043 ## size in mb flwo, cutoff for tr will be 10% of this?? arbitrary i know

genomebase='irugosBothHaps_aggressivecorrection'
genomesize=633*2 ## multiply for both haplotypes???
genomebase='rtuberHap1_aggressivecorrection'
genomesize=1751
genomebase='snutan_aggressivecorrection'
genomesize=2541*2
genomebase='achineBothHaps_aggressivecorrection'
genomesize=1449*2
genomebase='etripsHap1_aggressivecorrection'
genomesize=4600


## chrlengths
a=read.table(paste0(genomebase, '.FINAL.fa.fai'), header=F)
seqs=a$V1[a$V2>1e6]

## gc content
gc=read.table(paste0(genomebase, '.FINAL.1mbnuccontent.bed'), header=F)
gc=gc[gc$V1%in%seqs,]
gc$gc=gc$V5
gc$start=gc$V2
gc$n=gc$V10

## telomeres
tel=read.table(paste0(genomebase, '.tidk_telomeric_repeat_windows.tsv'), header=T)
tel=tel[tel$id%in%seqs,]
tel_long <- tel %>%
  pivot_longer(
    cols = c(forward_repeat_number, reverse_repeat_number), # Columns to make longer
    names_to = "direction",                                 # New column for the original column names
    values_to = "repeat_count"                              # New column for the values (the counts)
  ) %>%
  # Clean up the 'direction' names for better plotting (optional, but recommended)
  mutate(
    direction = case_when(
      direction == "forward_repeat_number" ~ "Forward",
      direction == "reverse_repeat_number" ~ "Reverse",
      TRUE ~ direction # Keep others as is
    )
  )


## ribosomes
rib=import.gff3(paste0(genomebase, '.FINAL.rrna.gff3'))
rib=rib[seqnames(rib)%in%seqs,]

## helixer genes
genes=import.gff3(paste0(genomebase, '.FINAL.helixer.gff3'))
genes=genes[seqnames(genes)%in%seqs,]
genes=genes[genes$type=='gene',]

## tandem repeats
#tr=read.table(paste0(genomebase, '_TRASH2/',genomebase, '.FINAL.fa_repeats_with_seq.csv'), header=T, sep=',')
tr=read.table(paste0(genomebase, '.FINAL.fa_repeats_with_seq.csv'), header=T, sep=',')
tr=tr[tr$seqID%in%seqs,]
#trash=import.gff3(paste0(genomebase, '_TRASH2/', genomebase, '.FINAL.fa_repeats.gff'))
#trash=trash[seqnames(trash)%in%seqs,]
## filter out short sequences and short arrays
tr$conswidth=as.numeric(str_split_fixed(tr$class, '_', 2)[,1])
tr%>%filter(conswidth>100)%>%group_by(arrayID)%>% mutate(n=n())%>%filter(n>20)
tr%>%group_by(arrayID)%>% mutate(n=n())%>%filter(n>20)%>%head%>%data.frame

## scale n by genome size???
trnumb=genomesize*0.1

trf=tr%>%group_by(arrayID, conswidth, seqID)%>% mutate(n=n())%>%filter(n>trnumb)

### ready to plot!!!

ribcols=viridis(4)
names(ribcols)=unique(rib$Name)

trfcols=plasma(length(unique(trf$class)))
names(trfcols)=unique(trf$class)



pdf(paste0('~/transfer/', genomebase, '_assembly.pdf'),20,length(seqs))
ggplot(gc, aes(x=start, y=gc))+geom_point()+facet_wrap(~factor(V1, levels=seqs), ncol=1)
ggplot(gc, aes(x=start, y=n))+geom_point()+facet_wrap(~factor(V1, levels=seqs), ncol=1)

ggplot(tel_long, aes(x=window, y=repeat_count, color=direction))+geom_point()+geom_line()+facet_wrap(~factor(id, levels=seqs), ncol=1)

 ggplot(trf, aes(x=start, y=class, color=class))+geom_point()+facet_wrap(~factor(seqID, levels=seqs), ncol=1)+scale_fill_manual(values=trfcols)
 ggplot(trf, aes(x=start, fill=class))+geom_histogram(binwidth=1e5, position='dodge')+facet_wrap(~factor(seqID, levels=seqs), ncol=1)+scale_fill_manual(values=trfcols)
 ggplot(trf, aes(x=start, fill=class))+geom_histogram(binwidth=1e6, position='dodge')+facet_wrap(~factor(seqID, levels=seqs), ncol=1)+scale_fill_manual(values=trfcols)

ggplot(data.frame(rib), aes(x=start, y=Name, color=Name))+geom_point()+facet_wrap(~factor(seqnames, levels=seqs), ncol=1)+scale_fill_manual(values=ribcols)
ggplot(data.frame(rib), aes(x=start, fill=Name))+geom_histogram(binwidth=1e5, position='dodge')+facet_wrap(~factor(seqnames, levels=seqs), ncol=1)+scale_fill_manual(values=ribcols)
ggplot(data.frame(rib), aes(x=start, fill=Name))+geom_histogram(binwidth=1e6, position='dodge')+facet_wrap(~factor(seqnames, levels=seqs), ncol=1)+scale_fill_manual(values=ribcols)

ggplot(data.frame(genes), aes(x=start))+geom_histogram(binwidth=1e6)+facet_wrap(~factor(seqnames, levels=seqs), ncol=1)


maxygenes=as.numeric(data.frame(genes)%>%group_by(round(start,-6), seqnames)%>%summarize(ngene=n())%>%ungroup()%>%summarize(max(ngene)))
maxyrib=as.numeric(data.frame(rib)%>%group_by(round(start,-6), seqnames)%>%summarize(ngene=n())%>%ungroup()%>%summarize(max(ngene)))
maxytrf=as.numeric(data.frame(trf)%>%group_by(round(start,-6), seqID)%>%summarize(ngene=n())%>%ungroup()%>%summarize(max(ngene)))
maxytel=max(tel_long$repeat_count)
maxns=max(gc$n)
mingc=min(gc$gc)
maxgc=max(gc$gc)

for(i in seqs){
print(
plot_grid(
ggplot(gc[gc$V1==i,], aes(x=start, y=gc))+geom_point()+ggtitle(i)+xlim(0,a$V2[a$V1==i])+ylim(mingc,maxgc),
ggplot(gc[gc$V1==i,], aes(x=start, y=n))+geom_point()+facet_wrap(~factor(V1, levels=seqs), ncol=1)+ylim(0,maxns),

ggplot(data.frame(genes[seqnames(genes)==i,]), aes(x=start))+geom_histogram(binwidth=1e6)+xlim(0,a$V2[a$V1==i])+ylim(0,maxygenes),
ggplot(data.frame(rib[seqnames(rib)==i,]), aes(x=start, fill=Name))+geom_histogram(binwidth=1e6, position='dodge')+xlim(0,a$V2[a$V1==i])+ylim(0,maxyrib)+scale_fill_manual(values=ribcols),
 ggplot(trf[trf$seqID==i,], aes(x=start, fill=class))+geom_histogram(binwidth=1e6, position='dodge')+xlim(0,a$V2[a$V1==i])+ylim(0,maxytrf)+scale_fill_manual(values=trfcols),
ggplot(tel_long[tel_long$id==i,], aes(x=window, y=repeat_count, color=direction))+geom_point()+geom_line()+xlim(0,a$V2[a$V1==i])+ylim(0,maxytel),
align='hv',
axis='b',
ncol=1,
rel_heights=c(0.5,0.5,1,0.5,0.5,0.5)

)
)
}

dev.off()