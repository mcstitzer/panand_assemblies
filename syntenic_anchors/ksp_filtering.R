# scp cbsublfs1:/data1/users/shengkaihsu/p_panAndOGASR/output/anchorAln_omega_v2/* .
### putting on ceres!!! /project/buckler_lab_panand/michelle.stitzer/anchorAln_omega_v2

library(data.table)

ksp=fread('cat ../anchorAln_omega_v2/Pavag*.fasta')
ksp=ksp[substr(ksp$V3,1,6)=='pvagin' | substr(ksp$V4,1,6)=='pvagin',] ## to paspalum
ksp$genome=substr(ksp$V3,1,6)                       


ksp=fread('cat Pavag*.fasta')
ksp=ksp[substr(ksp$V3,1,6)=='pvagin' | substr(ksp$V4,1,6)=='pvagin',] ## to paspalum
ksp$genome=substr(ksp$V3,1,6)                       

taxonnames=c("Z. mays ssp. parviglumis TIL11", "Z. mays ssp. mays B73v5", "Z. mays ssp. parviglumis TIL01", "Z. mays ssp. mexicana TIL25", "Z. mays ssp. mexicana TIL18", "Z. mays ssp. huehuetengensis", 
"Z. luxurians", "Z. nicaraguensis", "Z. diploperennis Momo", "Z. diploperennis Gigi", "T. zoloptense", "T. dactyloides FL", "T. dactyloides Southern Hap2", 
"T. dactyloides Northern Hap2", "T. dactyloides KS", "T. dactyloides tetraploid", "U. digitatum", "V. cuspidata", "R. rottboellioides", "R. tuberculosa", 
"H. compressa", "E. tripsacoides", "S. scoparium", "S. microstachyum", "A. virginicum", "A. chinensis", "A. gerardi", 
"C. refractus", "C. citratus", "H. contortus", "T. triandra", "B. laguroides", "P. paniceum", "S. bicolor", 
"I. rugosum", "S. nutans", '"A." burmanicus', "T. elegans", "C. serrulatus", "P. vaginatum")
names(taxonnames)=c("zTIL11", "zmB735", "zTIL01", "zTIL25", "zTIL18", "zmhuet", 
"zluxur", "znicar", "zdmomo", "zdgigi", "tzopol", "tdacs1", "tdacs2", 
"tdacn2", "tdacn1", "tdactm", "udigit", "vcuspi", "rrottb", "rtuber", 
"hcompr", "etrips", "sscopa", "smicro", "avirgi", "achine", "agerar", 
"crefra", "ccitra", "hconto", "ttrian", "blagur", "ppanic", "sbicol", 
"irugos", "snutan", "atenui", "telega", "cserru", "pvagin")
head(ksp)
summary(ksp$V19)
ksp$gene=str_split_fixed(ksp$V4, '_', 3)[,2]
library(stringr)
ksp$gene=str_split_fixed(ksp$V4, '_', 3)[,2]


write.table(ksp, 'ksp_to_pasp.txt', quote=F, sep='\t', row.names=F, col.names=T)


staydup=c("Pavag01G006300", "Pavag01G025500", "Pavag01G030900", "Pavag01G036800", 
"Pavag01G041500", "Pavag01G041700", "Pavag01G055800", "Pavag01G062000", 
"Pavag01G072100", "Pavag01G085300", "Pavag01G085400", "Pavag01G088500", 
"Pavag01G090500", "Pavag01G107500", "Pavag01G107600", "Pavag01G109500", 
"Pavag01G118700", "Pavag01G119300", "Pavag01G127800", "Pavag01G134900", 
"Pavag01G137000", "Pavag01G138500", "Pavag01G145000", "Pavag01G145100", 
"Pavag01G149200", "Pavag01G155300", "Pavag01G155500", "Pavag01G168900", 
"Pavag01G178700", "Pavag01G180600", "Pavag01G181400", "Pavag01G189700", 
"Pavag01G198200", "Pavag01G198300", "Pavag01G199600", "Pavag01G226000", 
"Pavag01G272700", "Pavag01G273000", "Pavag01G278900", "Pavag01G279000", 
"Pavag01G281400", "Pavag01G300000", "Pavag01G300500", "Pavag01G317100", 
"Pavag01G319000", "Pavag01G321400", "Pavag01G350200", "Pavag01G350400", 
"Pavag01G350600", "Pavag01G392300", "Pavag01G396400", "Pavag01G397300", 
"Pavag01G405400", "Pavag01G405500", "Pavag01G424600", "Pavag01G429400", 
"Pavag01G429600", "Pavag02G067200", "Pavag02G067400", "Pavag02G088000", 
"Pavag02G143200", "Pavag02G195800", "Pavag02G202100", "Pavag02G202600", 
"Pavag02G208300", "Pavag02G220100", "Pavag02G250900", "Pavag02G251300", 
"Pavag02G291300", "Pavag02G306400", "Pavag02G322700", "Pavag02G324100", 
"Pavag02G325300", "Pavag02G326700", "Pavag02G327600", "Pavag02G331500", 
"Pavag02G331600", "Pavag02G331700", "Pavag02G331800", "Pavag02G335600", 
"Pavag02G339700", "Pavag02G339800", "Pavag02G358500", "Pavag02G361300", 
"Pavag02G362600", "Pavag03G002000", "Pavag03G005200", "Pavag03G006800", 
"Pavag03G006900", "Pavag03G008800", "Pavag03G008900", "Pavag03G011000", 
"Pavag03G026500", "Pavag03G048800", "Pavag03G050400", "Pavag03G050500", 
"Pavag03G055100", "Pavag03G076900", "Pavag03G079400", "Pavag03G079800", 
"Pavag03G081200", "Pavag03G086200", "Pavag03G125300", "Pavag04G030000", 
"Pavag04G034500", "Pavag04G061000", "Pavag04G066300", "Pavag04G073900", 
"Pavag04G075700", "Pavag04G110900", "Pavag04G238900", "Pavag04G239800", 
"Pavag04G244000", "Pavag04G291400", "Pavag04G294000", "Pavag04G295200", 
"Pavag04G300700", "Pavag04G302100", "Pavag04G322400", "Pavag05G142100", 
"Pavag06G095300", "Pavag06G104800", "Pavag06G109600", "Pavag06G111100", 
"Pavag06G111500", "Pavag06G112700", "Pavag06G130100", "Pavag06G131100", 
"Pavag06G140600", "Pavag06G147500", "Pavag06G208300", "Pavag06G225000", 
"Pavag06G225300", "Pavag06G236000", "Pavag06G236100", "Pavag06G273900", 
"Pavag06G276400", "Pavag06G276600", "Pavag06G276900", "Pavag07G062200", 
"Pavag07G074300", "Pavag08G054700", "Pavag08G055200", "Pavag08G104900", 
"Pavag08G146200", "Pavag09G229200", "Pavag09G229300", "Pavag10G169300", 
"Pavag10G169400", "Pavag10G197400", "Pavag10G197900", "Pavag10G201100", 
"Pavag10G201200", "Pavag10G230400", "Pavag10G251700", "Pavag10G279100"
)