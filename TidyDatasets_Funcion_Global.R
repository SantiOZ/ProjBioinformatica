
#PREPARACION PREVIA DE DATASETS 

#Opening datasets
Shchetynsky <- read.table("C:/Users/tt_kb/Documents/Bioinformatica/Art_Shitinksy/Art.txt", header = TRUE)
orderShchet <- arrange(Shchetynsky, desc(EA1))
uniqueShchet <- orderShchet[!duplicated(orderShchet[,1]),]
rownames(uniqueShchet) <- uniqueShchet[,1]
uniqueShchet <- uniqueShchet[,-1]

Cheng <- read.table("C:/Users/tt_kb/Documents/Bioinformatica/Lupus_Cheng/Lupus.txt", header = TRUE)
rownames(Cheng) <- Cheng[,1]
Cheng <- Cheng[,-1]

MoSM <- read.csv("C:/Users/tt_kb/Documents/Bioinformatica/Crohn_Mo/GSE112057_series_matrix_2.csv", header = FALSE)
Mo <- read.table("C:/Users/tt_kb/Documents/Bioinformatica/Crohn_Mo/GSE112057_RawCounts_dataset.txt", header = TRUE)
orderMo <- arrange(Mo, desc(Mo[,2]))
uniqueMo <- orderMo[!duplicated(orderMo[,1]),]
rownames(uniqueMo) <- uniqueMo[,1]
Mo <- uniqueMo[,-1]

#Opening index 38.98 GTF file
index <- as.data.frame(rtracklayer::import('C:/Users/tt_kb/Documents/Bioinformatica/Index.gtf'))

#Matching datasets gene id with Index
Cheng_match <- match(rownames(Cheng), index$gene_id)
Cheng_match_clean <- Cheng_match[!is.na(Cheng_match)]

Shchet_match <- match(rownames(uniqueShchet), index$gene_id)
Shchet_match_clean <- Shchet_match[!is.na(Shchet_match)]

Mo_match <- match(rownames(Mo), index$gene_name)
Mo_match_clean <- Mo_match[!is.na(Mo_match)]

#Obtener nombres a partir de gene id
geneNamesCheng <- index$gene_name[Cheng_match_clean]
geneNamesShchet <- index$gene_name[Shchet_match_clean]
geneNamesMo <- index$gene_name[Mo_match_clean]
geneIDMo <- index$gene_id[Mo_match_clean]

#Eliminar filas sin match con index en dataset original
Cheng <- Cheng[!is.na(Cheng_match),]
uniqueShchet <- uniqueShchet[!is.na(Shchet_match),]
Mo <- Mo[!is.na(Mo_match),]
rownames(Mo) <- geneIDMo

#Separar dataset de Mo en grupos de enfermedades

Control_Mo <- Mo[grep("Control", MoSM[1,-1])]
colnames(Control_Mo) <- c("SM1", "SM2", "SM3", "SM4", "SM5", "SM6", "SM7", "SM8", "SM9", "SM10", "SM11", "SM12")

Crohn <- Mo[,grep("Crohn", MoSM[1,-1])]
colnames(Crohn) <- c("EC1", "EC2",	"EC3","EC4","EC5","EC6","EC7","EC8","EC9","EC10","EC11","EC12",	"EC13",	"EC14",	
                     "EC15","EC16","EC17","EC18",	"EC19","EC20","EC21","EC22","EC23","EC24","EC25","EC26","EC27",
                     "EC28","EC29","EC30","EC31","EC32","EC33","EC34","EC35","EC36","EC37","EC38","EC39","EC40",
                     "EC41",	"EC42","EC43","EC44","EC45","EC46","EC47","EC48", "EC49","EC50","EC51","EC52","EC53",
                     "EC54",	"EC55",	"EC56",	"EC57",	"EC58",	"EC59",	"EC60")


OligoJIA <- Mo[grep("Oligoarticular", MoSM[1,-1])]
colnames(OligoJIA) <- c("EO1",	"EO2",	"EO3",	"EO4",	"EO5",	"EO6",	"EO7",	"EO8",	"EO9",	"EO10",	"EO11",	"EO12",	
                        "EO13",	"EO14",	"EO15",	"EO16",	"EO17",	"EO18",	"EO19",	"EO20",	"EO21",	"EO22",	"EO23",	"EO24",	
                        "EO25",	"EO26",	"EO27",	"EO28",	"EO29",	"EO30",	"EO31",	"EO32",	"EO33",	"EO34",	"EO35",	"EO36",	
                        "EO37",	"EO38",	"EO39",	"EO40",	"EO41",	"EO42",	"EO43")

Ulc_Col <- Mo[grep("Colitis", MoSM[1,-1])]
colnames(Ulc_Col) <- c("EU1",	"EU2",	"EU3",	"EU4",	"EU5",	"EU6",	"EU7",	"EU8",	"EU9",	"EU10",	"EU11",	"EU12",	
                       "EU13",	"EU14",	"EU15")


PolyJIA<- Mo[grep("Polyarticular", MoSM[1,-1])]
colnames(PolyJIA) <- c("EP1",	"EP2",	"EP3",	"EP4",	"EP5",	"EP6",	"EP7",	"EP8",	"EP9",	"EP10",	"EP11",	"EP12",	
                       "EP13",	"EP14",	"EP15",	"EP16",	"EP17",	"EP18",	"EP19",	"EP20",	"EP21",	"EP22",	"EP23",	"EP24",	
                       "EP25",	"EP26",	"EP27",	"EP28",	"EP29",	"EP30",	"EP31",	"EP32",	"EP33",	"EP34",	"EP35",	"EP36",	
                       "EP37",	"EP38",	"EP39",	"EP40",	"EP41",	"EP42",	"EP43",	"EP44",	"EP45",	"EP46")

SystJIA<- Mo[grep("Systemic", MoSM[1,-1])]
colnames(SystJIA) <- c("ES1",	"ES2",	"ES3",	"ES4",	"ES5",	"ES6",	"ES7",	"ES8",	"ES9",	"ES10",	"ES11",	"ES12",	
                       "ES13",	"ES14",	"ES15",	"E2S16",	"ES17",	"ES18",	"ES19",	"ES20",	"ES21",	"ES22",	"ES23",	"ES24",	
                       "ES25",	"ES26")

#Agregar columna de gene names
Cheng <- cbind(geneNamesCheng,Cheng)
Shchetynsky <- cbind(geneNamesShchet, uniqueShchet)
Mo_Crohn <- cbind(geneNamesMo, Control_Mo, Crohn)
Mo_OligoJIA <- cbind(geneNamesMo, Control_Mo, OligoJIA)
Mo_PolyJIA <- cbind(geneNamesMo, Control_Mo, PolyJIA)
Mo_SystJIA <- cbind(geneNamesMo, Control_Mo, SystJIA)
Mo_Ulc_Col <- cbind(geneNamesMo, Control_Mo, Ulc_Col)

#Quitando ceros
emptyCheng <- Cheng[-which(apply(Cheng[,-1], 1, mean) == 0),]
emptyShchet <- Shchetynsky[-which(apply(Shchetynsky[,-1], 1, mean) == 0),]
emptyMo_Crohn <- Mo_Crohn[-which(apply(Mo_Crohn[,-1], 1, mean) == 0),]
emptyMo_OligoJIA <- Mo_OligoJIA[-which(apply(Mo_OligoJIA[,-1], 1, mean) == 0),]
emptyMo_PolyJIA <- Mo_PolyJIA[-which(apply(Mo_PolyJIA[,-1], 1, mean) == 0),]
emptyMo_SystJIA <- Mo_SystJIA[-which(apply(Mo_SystJIA[,-1], 1, mean) == 0),]
emptyMo_Ulc_Col <- Mo_Ulc_Col[-which(apply(Mo_Ulc_Col[,-1], 1, mean) == 0),]