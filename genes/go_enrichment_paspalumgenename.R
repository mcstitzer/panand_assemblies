

library(topGO)


below=c('snutan', 'vcuspi', 'hcompr', 'agerar', 'sscopa', 'ccitra', 'atenui')
above=c('achine', 'blagur', 'hconto', 'etrips', 'udigit', 'rrottb')

lowerbelow=gsub('.v3.1', '', b$gene[rowSums(b[,below])<rowSums(b[,above]) ])
b$gene[rowSums(b[,below])<rowSums(b[,above]) ]


## from shengkai



GODB <- readMappings('~/Downloads/Pv_GOTable.txt',IDsep = ",")
background = names(GODB)

goi = lowerbelow
# goi = strsplit2(goi,"...v3")[,1]

tmp=factor(as.integer(background%in%goi))
names(tmp)=background
tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.gene2GO, gene2GO = GODB)
resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
GO_res_table=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)


allGO=genesInTerm(tgd1)
lapply(allGO[GO_res_table$GO.ID[1:3]] , function(x) x[x%in%goi])

## get paths to get the genes!
paste0('scp mcs368@cbsublfs1.tc.cornell.edu:/data1/users/mcs368/panand_gene_trees/aln_with_paspalumCDS_NOEMPTY/',
        substr(unlist(lapply(allGO[GO_res_table$GO.ID[1:3]] , function(x) x[x%in%goi])), 1,14), '.withPvCDS.aln.fa .')






