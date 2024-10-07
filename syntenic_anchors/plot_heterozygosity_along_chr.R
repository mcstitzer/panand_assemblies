

tdhet=read.table('~/Downloads/tdacs1_snp_summary_100kb.bed')
ggplot(tdhet[substr(tdhet$V1,1,1)=='c',], aes(x=V2, y=V4/100000)) + geom_point() + facet_wrap(~V1, ncol = 1)
ggplot(tdhet[substr(tdhet$V1,1,1)=='c',], aes(x=V2, y=V4/100000)) + geom_point() + facet_wrap(~V1, ncol = 2)


cents=fread(cmd='grep -E "156bp_repeat|155bp_repeat" ~/Downloads/repeatmask_tandems/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a_EDTAandTandemRepeat.gff3')
knobs=fread(cmd='grep -E "179bp_repeat|180bp_repeat" ~/Downloads/repeatmask_tandems/Td-FL_9056069_6-REFERENCE-PanAnd-2.0a_EDTAandTandemRepeat.gff3')


cents <- cents %>%
  mutate(V4 = as.integer(V4), V5 = as.integer(V5)) %>%  # Ensure start and end are integers
  arrange(V1, V4) %>%  # Sort by chromosome and start position
  group_by(V1) %>%  # Group by chromosome
  mutate(
    distance_to_prev = V4 - lag(V5, default = first(V4)),  # Calculate the distance to the previous entry
    group_id = cumsum(distance_to_prev > 100000)  # Create a group ID, start a new group if the distance > 100,000
  ) %>%
  group_by(V1, group_id) %>%
  summarize(
    V4 = min(V4),  # Take the smallest start in the group
    V5 = max(V5),  # Take the largest end in the group
    V6 = mean(V6),  # Optionally average similarity (or adjust as needed)
    V9 = paste(unique(V9), collapse = "; "),  # Combine the motif information
    .groups = 'drop'
  ) %>% filter(substr(V1,1,1)=='c')

cents_extended <- cents %>%
  mutate(
    V4 = pmax(0, V4 - 100000),  # Extend the start position by 100,000 bp, ensuring it doesn't go below 0
    V5 = V5 + 100000            # Extend the end position by 100,000 bp
  )


knobs <- knobs %>%
  mutate(V4 = as.integer(V4), V5 = as.integer(V5)) %>%  # Ensure start and end are integers
  arrange(V1, V4) %>%  # Sort by chromosome and start position
  group_by(V1) %>%  # Group by chromosome
  mutate(
    distance_to_prev = V4 - lag(V5, default = first(V4)),  # Calculate the distance to the previous entry
    group_id = cumsum(distance_to_prev > 100000)  # Create a group ID, start a new group if the distance > 100,000
  ) %>%
  group_by(V1, group_id) %>%
  summarize(
    V4 = min(V4),  # Take the smallest start in the group
    V5 = max(V5),  # Take the largest end in the group
    V6 = mean(V6),  # Optionally average similarity (or adjust as needed)
    V9 = paste(unique(V9), collapse = "; "),  # Combine the motif information
    .groups = 'drop'
  ) %>% filter(substr(V1,1,1)=='c', V5-V4 > 100000) 

knobs_extended <- knobs %>%
  mutate(
    V4 = pmax(0, V4 - 100000),  # Extend the start position by 100,000 bp, ensuring it doesn't go below 0
    V5 = V5 + 100000            # Extend the end position by 100,000 bp
  )

ggplot(tdhet[substr(tdhet$V1,1,1)=='c',], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~factor(V1, levels=paste0('chr', c(1:18))), ncol = 2) + geom_rect(data=cents_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkblue', alpha=0.9)+ geom_rect(data=knobs_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkgreen', alpha=0.3)  + ylab('Heterozygosity') + xlab('T. dactyloides FL Hap 1') + ggtitle('T. dactyloides Southern')
## plot 0s for jeff
ggplot(tdhet[substr(tdhet$V1,1,1)=='c',], aes(x=V2, y=V4/100000, color=V4==0)) + geom_point(alpha=0.2) + scale_color_manual(values=c('black', 'red'))+ facet_wrap(~factor(V1, levels=paste0('chr', c(1:18))), ncol = 2) + geom_rect(data=cents_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), color=NA, fill='darkblue', alpha=0.9)+ geom_rect(data=knobs_extended, color=NA, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkgreen', alpha=0.3)  + ylab('Heterozygosity') + xlab('T. dactyloides FL Hap 1') + ggtitle('T. dactyloides Southern' )+ theme(legend.position='NULL')


tdhet=read.table('~/Downloads/tdacn1_snp_summary_100kb.bed')
ggplot(tdhet[substr(tdhet$V1,1,1)=='c',], aes(x=V2, y=V4/100000)) + geom_point() + facet_wrap(~V1, ncol = 1)
ggplot(tdhet[substr(tdhet$V1,1,1)=='c',], aes(x=V2, y=V4/100000)) + geom_point() + facet_wrap(~V1, ncol = 2)


cents=fread(cmd='grep -E "156bp_repeat|155bp_repeat" ~/Downloads/repeatmask_tandems/Td-KS_B6_1-REFERENCE-PanAnd-2.0a_EDTAandTandemRepeat.gff3')
knobs=fread(cmd='grep -E "179bp_repeat|180bp_repeat" ~/Downloads/repeatmask_tandems/Td-KS_B6_1-REFERENCE-PanAnd-2.0a_EDTAandTandemRepeat.gff3')


cents <- cents %>%
  mutate(V4 = as.integer(V4), V5 = as.integer(V5)) %>%  # Ensure start and end are integers
  arrange(V1, V4) %>%  # Sort by chromosome and start position
  group_by(V1) %>%  # Group by chromosome
  mutate(
    distance_to_prev = V4 - lag(V5, default = first(V4)),  # Calculate the distance to the previous entry
    group_id = cumsum(distance_to_prev > 100000)  # Create a group ID, start a new group if the distance > 100,000
  ) %>%
  group_by(V1, group_id) %>%
  summarize(
    V4 = min(V4),  # Take the smallest start in the group
    V5 = max(V5),  # Take the largest end in the group
    V6 = mean(V6),  # Optionally average similarity (or adjust as needed)
    V9 = paste(unique(V9), collapse = "; "),  # Combine the motif information
    .groups = 'drop'
  ) %>% filter(substr(V1,1,1)=='c')

cents_extended <- cents %>%
  mutate(
    V4 = pmax(0, V4 - 100000),  # Extend the start position by 100,000 bp, ensuring it doesn't go below 0
    V5 = V5 + 100000            # Extend the end position by 100,000 bp
  )

knobs <- knobs %>%
  mutate(V4 = as.integer(V4), V5 = as.integer(V5)) %>%  # Ensure start and end are integers
  arrange(V1, V4) %>%  # Sort by chromosome and start position
  group_by(V1) %>%  # Group by chromosome
  mutate(
    distance_to_prev = V4 - lag(V5, default = first(V4)),  # Calculate the distance to the previous entry
    group_id = cumsum(distance_to_prev > 100000)  # Create a group ID, start a new group if the distance > 100,000
  ) %>%
  group_by(V1, group_id) %>%
  summarize(
    V4 = min(V4),  # Take the smallest start in the group
    V5 = max(V5),  # Take the largest end in the group
    V6 = mean(V6),  # Optionally average similarity (or adjust as needed)
    V9 = paste(unique(V9), collapse = "; "),  # Combine the motif information
    .groups = 'drop'
  ) %>% filter(substr(V1,1,1)=='c', V5-V4 > 100000)

knobs_extended <- knobs %>%
  mutate(
    V4 = pmax(0, V4 - 100000),  # Extend the start position by 100,000 bp, ensuring it doesn't go below 0
    V5 = V5 + 100000            # Extend the end position by 100,000 bp
  )

ggplot(tdhet[substr(tdhet$V1,1,1)=='c',], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~factor(V1, levels=paste0('chr', c(1:18))), ncol = 2) + geom_rect(data=cents_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkblue', alpha=0.9) + geom_rect(data=knobs_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkgreen', alpha=0.3) + ylab('Heterozygosity') + xlab('T. dactyloides KS Hap 1') + ggtitle('T. dactyloides Northern')

## plot 0s for jeff
ggplot(tdhet[substr(tdhet$V1,1,1)=='c',], aes(x=V2, y=V4/100000, color=V4==0)) + geom_point(alpha=0.2) + scale_color_manual(values=c('black', 'red'))+ facet_wrap(~factor(V1, levels=paste0('chr', c(1:18))), ncol = 2) + geom_rect(data=cents_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), color=NA, fill='darkblue', alpha=0.9)+ geom_rect(data=knobs_extended, color=NA, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkgreen', alpha=0.3)  + xlab('T. dactyloides KS Hap 1') + ggtitle('T. dactyloides Northern') + ggtitle('T. dactyloides Southern' )+ theme(legend.position='NULL')



#####
tdhet=read.table('~/Downloads/avirgi_snp_summary_100kb.bed')
ggplot(tdhet[substr(tdhet$V1,1,1)=='c',], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~V1, ncol = 2)

tdhet=read.table('~/Downloads/ppanic_snp_summary_100kb.bed') %>% group_by(V1)%>% mutate(max=max(V3))
ggplot(tdhet[tdhet$max>10e6,], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~V1, ncol = 2)

tdhet=read.table('~/Downloads/rrottb_snp_summary_100kb.bed') %>% group_by(V1)%>% mutate(max=max(V3))
ggplot(tdhet[tdhet$max>10e6,], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~V1, ncol = 2)

tdhet=read.table('~/Downloads/achine_snp_summary_100kb.bed') %>% group_by(V1)%>% mutate(max=max(V3))
ggplot(tdhet[tdhet$max>10e6,], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~V1, ncol = 2)

tdhet=read.table('~/Downloads/smicro_snp_summary_100kb.bed') %>% group_by(V1)%>% mutate(max=max(V3))
ggplot(tdhet[tdhet$max>10e6,], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~V1, ncol = 2)

tdhet=read.table('~/Downloads/udigit_snp_summary_100kb.bed') %>% group_by(V1)%>% mutate(max=max(V3))
ggplot(tdhet[tdhet$max>10e6,], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~V1, ncol = 2)


tdhet=read.table('~/Downloads/zmhuet_snp_summary_100kb.bed') %>% group_by(V1)%>% mutate(max=max(V3))
ggplot(tdhet[tdhet$max>10e6,], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~V1, ncol = 2)

cents=fread(cmd='grep -E "156bp_repeat|155bp_repeat" ~/Downloads/repeatmask_tandems/Zh-RIMHU001-REFERENCE-PanAnd-1.0_EDTAandTandemRepeat.gff3')
knobs=fread(cmd='grep -E "179bp_repeat|180bp_repeat" ~/Downloads/repeatmask_tandems/Zh-RIMHU001-REFERENCE-PanAnd-1.0_EDTAandTandemRepeat.gff3')


cents <- cents %>%
  mutate(V4 = as.integer(V4), V5 = as.integer(V5)) %>%  # Ensure start and end are integers
  arrange(V1, V4) %>%  # Sort by chromosome and start position
  group_by(V1) %>%  # Group by chromosome
  mutate(
    distance_to_prev = V4 - lag(V5, default = first(V4)),  # Calculate the distance to the previous entry
    group_id = cumsum(distance_to_prev > 100000)  # Create a group ID, start a new group if the distance > 100,000
  ) %>%
  group_by(V1, group_id) %>%
  summarize(
    V4 = min(V4),  # Take the smallest start in the group
    V5 = max(V5),  # Take the largest end in the group
    V6 = mean(V6),  # Optionally average similarity (or adjust as needed)
    V9 = paste(unique(V9), collapse = "; "),  # Combine the motif information
    .groups = 'drop'
  ) %>% filter(substr(V1,1,3)=='chr')

cents_extended <- cents %>%
  mutate(
    V4 = pmax(0, V4 - 100000),  # Extend the start position by 100,000 bp, ensuring it doesn't go below 0
    V5 = V5 + 100000            # Extend the end position by 100,000 bp
  )


knobs <- knobs %>%
  mutate(V4 = as.integer(V4), V5 = as.integer(V5)) %>%  # Ensure start and end are integers
  arrange(V1, V4) %>%  # Sort by chromosome and start position
  group_by(V1) %>%  # Group by chromosome
  mutate(
    distance_to_prev = V4 - lag(V5, default = first(V4)),  # Calculate the distance to the previous entry
    group_id = cumsum(distance_to_prev > 100000)  # Create a group ID, start a new group if the distance > 100,000
  ) %>%
  group_by(V1, group_id) %>%
  summarize(
    V4 = min(V4),  # Take the smallest start in the group
    V5 = max(V5),  # Take the largest end in the group
    V6 = mean(V6),  # Optionally average similarity (or adjust as needed)
    V9 = paste(unique(V9), collapse = "; "),  # Combine the motif information
    .groups = 'drop'
  ) %>% filter(substr(V1,1,3)=='chr', V5-V4 > 100000) 

knobs_extended <- knobs %>%
  mutate(
    V4 = pmax(0, V4 - 100000),  # Extend the start position by 100,000 bp, ensuring it doesn't go below 0
    V5 = V5 + 100000            # Extend the end position by 100,000 bp
  )

ggplot(tdhet[substr(tdhet$V1,1,3)=='chr',], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~factor(V1, levels=paste0('chr', c(1:10))), ncol = 2) + geom_rect(data=cents_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkblue', alpha=0.9)+ geom_rect(data=knobs_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkgreen', alpha=0.3)  + ylab('Heterozygosity') + xlab('Z. m. subsp. huehue.') 
## plot 0s for jeff
ggplot(tdhet[substr(tdhet$V1,1,3)=='chr',], aes(x=V2, y=V4/100000, color=V4==0)) + geom_point(alpha=0.2) + scale_color_manual(values=c('black', 'red'))+ facet_wrap(~factor(V1, levels=paste0('chr', c(1:18))), ncol = 2) + geom_rect(data=cents_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), color=NA, fill='darkblue', alpha=0.9)+ geom_rect(data=knobs_extended, color=NA, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkgreen', alpha=0.3)  + ylab('Heterozygosity') + xlab('Z. m. subsp. huehue.') + theme(legend.position='NULL')





###

tdhet=read.table('~/Downloads/zTIL25_snp_summary_100kb.bed') %>% group_by(V1)%>% mutate(max=max(V3))
ggplot(tdhet[tdhet$max>10e6,], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~V1, ncol = 2)

cents=fread(cmd='grep -E "156bp_repeat|155bp_repeat" ~/Downloads/repeatmask_tandems/Zx-TIL25-REFERENCE-PanAnd-1.0_EDTA.gff3')
knobs=fread(cmd='grep -E "179bp_repeat|180bp_repeat" ~/Downloads/repeatmask_tandems/Zx-TIL25-REFERENCE-PanAnd-1.0_EDTA.gff3)


cents <- cents %>%
  mutate(V4 = as.integer(V4), V5 = as.integer(V5)) %>%  # Ensure start and end are integers
  arrange(V1, V4) %>%  # Sort by chromosome and start position
  group_by(V1) %>%  # Group by chromosome
  mutate(
    distance_to_prev = V4 - lag(V5, default = first(V4)),  # Calculate the distance to the previous entry
    group_id = cumsum(distance_to_prev > 100000)  # Create a group ID, start a new group if the distance > 100,000
  ) %>%
  group_by(V1, group_id) %>%
  summarize(
    V4 = min(V4),  # Take the smallest start in the group
    V5 = max(V5),  # Take the largest end in the group
    V6 = mean(V6),  # Optionally average similarity (or adjust as needed)
    V9 = paste(unique(V9), collapse = "; "),  # Combine the motif information
    .groups = 'drop'
  ) %>% filter(substr(V1,1,3)=='chr')

cents_extended <- cents %>%
  mutate(
    V4 = pmax(0, V4 - 100000),  # Extend the start position by 100,000 bp, ensuring it doesn't go below 0
    V5 = V5 + 100000            # Extend the end position by 100,000 bp
  )


knobs <- knobs %>%
  mutate(V4 = as.integer(V4), V5 = as.integer(V5)) %>%  # Ensure start and end are integers
  arrange(V1, V4) %>%  # Sort by chromosome and start position
  group_by(V1) %>%  # Group by chromosome
  mutate(
    distance_to_prev = V4 - lag(V5, default = first(V4)),  # Calculate the distance to the previous entry
    group_id = cumsum(distance_to_prev > 100000)  # Create a group ID, start a new group if the distance > 100,000
  ) %>%
  group_by(V1, group_id) %>%
  summarize(
    V4 = min(V4),  # Take the smallest start in the group
    V5 = max(V5),  # Take the largest end in the group
    V6 = mean(V6),  # Optionally average similarity (or adjust as needed)
    V9 = paste(unique(V9), collapse = "; "),  # Combine the motif information
    .groups = 'drop'
  ) %>% filter(substr(V1,1,3)=='chr', V5-V4 > 100000) 

knobs_extended <- knobs %>%
  mutate(
    V4 = pmax(0, V4 - 100000),  # Extend the start position by 100,000 bp, ensuring it doesn't go below 0
    V5 = V5 + 100000            # Extend the end position by 100,000 bp
  )

ggplot(tdhet[substr(tdhet$V1,1,3)=='chr',], aes(x=V2, y=V4/100000)) + geom_point(alpha=0.2) + facet_wrap(~factor(V1, levels=paste0('chr', c(1:10))), ncol = 2) + geom_rect(data=cents_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkblue', alpha=0.9)+ geom_rect(data=knobs_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkgreen', alpha=0.3)  + ylab('Heterozygosity') + xlab('Z. m. subsp. mexicana TIL25') 
## plot 0s for jeff
ggplot(tdhet[substr(tdhet$V1,1,3)=='chr',], aes(x=V2, y=V4/100000, color=V4==0)) + geom_point(alpha=0.2) + scale_color_manual(values=c('black', 'red'))+ facet_wrap(~factor(V1, levels=paste0('chr', c(1:18))), ncol = 2) + geom_rect(data=cents_extended, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), color=NA, fill='darkblue', alpha=0.9)+ geom_rect(data=knobs_extended, color=NA, aes(x=V4, y=0,xmin=V4, xmax=V5, ymin=-Inf, ymax=Inf), fill='darkgreen', alpha=0.3)  + ylab('Heterozygosity') + xlab('Z. m. subsp. mexicana TIL25') + theme(legend.position='NULL')
