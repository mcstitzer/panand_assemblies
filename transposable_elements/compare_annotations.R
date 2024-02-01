
all=read.table('../panand_sp_ploidy.txt', header=F)
fil=Sys.glob('earlGreyOutputs/*/*/*.high*.txt')
all$eg=NA
all$bp=NA
all$ltr=NA
for(i in 1:length(fil)){
 b=read.table(fil[i], sep='\t', header=T)
all$eg[all$V2==substr(fil[i],17,22)]=sum(b$proportion)
all$bp[all$V2==substr(fil[i],17,22)]=sum(b$cov)
all$ltr[all$V2==substr(fil[i],17,22)]=b$cov[b$tclassif=='LTR']
}


dput(all)
structure(list(V1 = c("Ac-Pasquet1232-DRAFT-PanAnd-1.0", "Ag-CAM1351-DRAFT-PanAnd-1.0", 
"Ab-Traiperm_572-DRAFT-PanAnd-1.0", "Av-Kellogg1287_8-REFERENCE-PanAnd-1.0", 
"Bl-K1279B-DRAFT-PanAnd-1.0", "Cs-KelloggPI219580-DRAFT-PanAnd-1.0", 
"Cc-PI314907-DRAFT-PanAnd-1.0", "Cr-AUB069-DRAFT-PanAnd-1.0", 
"Et-Layton_Zhong168-DRAFT-PanAnd-1.0", "Hp-KelloggPI404118-DRAFT-PanAnd-1.0", 
"Hc-AUB53_1-DRAFT-PanAnd-1.0", "Ir-Pasquet1136-DRAFT-PanAnd-1.0", 
"Pp-Kellogg1297-DRAFT-PanAnd-1.0", "Pi-Clark-DRAFT-PanAnd-1.0", 
"Rr-Malcomber3106-DRAFT-PanAnd-1.0", "Rt-Layton_Zhong169-DRAFT-PanAnd-1.0", 
"Sm-PI203595-DRAFT-PanAnd-1.0", "Ss-CAM1384-DRAFT-PanAnd-1.0", 
"Sn-CAM1369-DRAFT-PanAnd-1.0", "Te-Pasquet1246-DRAFT-PanAnd-1.0", 
"Tt-AUB21_1-DRAFT-PanAnd-1.0", "Td-McKain334_5-DRAFT-PanAnd-1.0", 
"Td-KS_B6_1-REFERENCE-PanAnd-2.0a", "Td-KS_B6_1-REFERENCE-PanAnd-2.0b", 
"Td-FL_9056069_6-REFERENCE-PanAnd-2.0a", "Td-FL_9056069_6-REFERENCE-PanAnd-2.0b", 
"Tz-DC_05_58_3A-DRAFT-PanAnd-1.0", "Ud-Pasquet1171-DRAFT-PanAnd-1.0", 
"Vc-Pasquet1098-DRAFT-PanAnd-1.0", "Zd-Gigi-REFERENCE-PanAnd-1.0", 
"Zd-Momo-REFERENCE-PanAnd-1.0", "Zl-RIL003-REFERENCE-PanAnd-1.0", 
"Zh-RIMHU001-REFERENCE-PanAnd-1.0", "Zx-TIL18-REFERENCE-PanAnd-1.0", 
"Zx-TIL25-REFERENCE-PanAnd-1.0", "Zv-TIL01-REFERENCE-PanAnd-1.0", 
"Zv-TIL11-REFERENCE-PanAnd-1.0", "Zn-PI615697-REFERENCE-PanAnd-1.0", 
"Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel", "Oryza_sativa.IRGSP-1.0.dna.toplevel", 
"Brachypodium_distachyon.Brachypodium_distachyon_v3.0.dna.toplevel", 
"Zm-B73-REFERENCE-NAM-5.0", "Andropogon_gerardii_var_Kellogg_1272.mainGenome.bothhaps", 
"Sviridis_726_v4.0", "Eremochloa_ophiuroides_genome_rc_oneline_rmNN"
), V2 = c("achine", "agerar", "atenui", "avirgi", "blagur", "cserru", 
"ccitra", "crefra", "etrips", "hcompr", "hconto", "irugos", "pprate", 
"ppanic", "rrottb", "rtuber", "smicro", "sscopa", "snutan", "telega", 
"ttrian", "tdactm", "tdacn1", "tdacn2", "tdacs1", "tdacs2", "tzopol", 
"udigit", "vcuspi", "zdgigi", "zdmomo", "zluxur", "zmhuet", "zTIL18", 
"zTIL25", "zTIL01", "zTIL11", "znicar", "sbicol", "osativ", "bdista", 
"zmB735", "agerjg", "svirid", "eophiu"), V3 = c(2L, 3L, 2L, 1L, 
3L, 1L, 3L, 1L, 2L, 3L, 2L, 1L, 4L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 
1L, 4L, 1L, 1L, 1L, 1L, 2L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 1L, 1L, 1L, 2L, 3L, 1L, 1L), eg = c(NA, NA, 0.602890037642115, 
0.702702665489812, NA, NA, NA, 0.69355983773986, NA, NA, 0.506618373396119, 
0.600423435073552, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
NA, NA, NA, NA, NA, NA, 0.827628578988442, 0.834867564866408, 
0.863062980685043, 0.841361942795097, NA, 0.850131057654226, 
NA, 0.855058207943679, 0.867982757381689, 0.669235302592893, 
0.500297199606712, 0.383365965746287, 0.845802670977004, NA, 
NA, NA), bp = c(NA, NA, 1173470913, 636418685, NA, NA, NA, 563815189, 
NA, NA, 1060842769, 423206173, NA, NA, NA, NA, NA, NA, NA, NA, 
NA, NA, NA, NA, NA, NA, NA, NA, NA, 2723660942, 2092415432, 2719377560, 
2021297883, NA, 1795758863, NA, 2044764362, 3173071490, 474310695, 
187636107, 103954826, 1845605704, NA, NA, NA)), row.names = c(NA, 
-45L), class = "data.frame")

## over on cbsu
te=read.table('../transposable_elements/total_repeat_bp.txt', header=T, sep='\t')
eg$edta=te$repeatbp[match(eg$V2, te$genome)]
cor(eg$bp, eg$edta, use='complete')




#"/local/workdir/mcs368/panand_htt/transposable_elements"
#write.table(all, 'all_trying_labmeeting.txt', quote=F, sep='\t', row.names=F, col.names=T)
