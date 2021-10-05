#datasets
#-artritis Schetynsky
#-lupus Cheng
#-crohn Mo
#-myositis Parkes

#PREPARACION PREVIA DE DATASETS 
library("dplyr")


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

ParkesDM <- read.csv("C:/Users/tt_kb/Documents/Bioinformatica/Myositis_Parkes/Dermatomyositis.csv", header = TRUE)
reducedID_DM <- substring(ParkesDM[,1],1,15)
ParkesDM <- ParkesDM[,-1]
rownames(ParkesDM) <- reducedID_DM

ParkesIBM <- read.csv("C:/Users/tt_kb/Documents/Bioinformatica/Myositis_Parkes/Inclusion_Body_Myositis.csv", header =  TRUE)
reducedID_IBM <- substring(ParkesIBM[,1],1,15)
ParkesIBM <- ParkesIBM[,-1]
rownames(ParkesIBM) <- reducedID_IBM

ParkesPM <- read.csv("C:/Users/tt_kb/Documents/Bioinformatica/Myositis_Parkes/Polymyositis.csv", header =  TRUE)
reducedID_PM <- substring(ParkesPM[,1],1,15)
ParkesPM <- ParkesPM[,-1]
rownames(ParkesPM) <- reducedID_PM


#Opening index 38.98 GTF file
index <- as.data.frame(rtracklayer::import('C:/Users/tt_kb/Documents/Bioinformatica/Index.gtf'))

#Matching datasets gene id with Index
Cheng_match <- match(rownames(Cheng), index$gene_id)
Cheng_match_clean <- Cheng_match[!is.na(Cheng_match)]

Shchet_match <- match(rownames(uniqueShchet), index$gene_id)
Shchet_match_clean <- Shchet_match[!is.na(Shchet_match)]

Mo_match <- match(rownames(Mo), index$gene_name)
Mo_match_clean <- Mo_match[!is.na(Mo_match)]

ParkesDM_match <- match(rownames(ParkesDM), index$gene_id)
ParkesDM_match_clean <- ParkesDM_match[!is.na(ParkesDM_match)]

ParkesIBM_match <- match(rownames(ParkesIBM), index$gene_id)
ParkesIBM_match_clean <- ParkesIBM_match[!is.na(ParkesIBM_match)]

ParkesPM_match <- match(rownames(ParkesPM), index$gene_id)
ParkesPM_match_clean <- ParkesPM_match[!is.na(ParkesPM_match)]


#Obtener nombres a partir de gene id
geneNamesCheng <- index$gene_name[Cheng_match_clean]
geneNamesShchet <- index$gene_name[Shchet_match_clean]
geneNamesMo <- index$gene_name[Mo_match_clean]
geneIDMo <- index$gene_id[Mo_match_clean]
geneNamesParkesDM <- index$gene_name[ParkesDM_match_clean]
geneNamesParkesIBM <- index$gene_name[ParkesIBM_match_clean]
geneNamesParkesPM <- index$gene_name[ParkesPM_match_clean]

#Eliminar filas sin match con index en dataset original
Cheng <- Cheng[!is.na(Cheng_match),] 
uniqueShchet <- uniqueShchet[!is.na(Shchet_match),] 
Mo <- Mo[!is.na(Mo_match),]
rownames(Mo) <- geneIDMo
ParkesDM <- ParkesDM[!is.na(ParkesDM_match),]
ParkesIBM <- ParkesIBM[!is.na(ParkesIBM_match),]
ParkesPM <- ParkesPM[!is.na(ParkesPM_match),]

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
colnames(SystJIA) <- c("EZ1",	"EZ2",	"EZ3",	"EZ4",	"EZ5",	"EZ6",	"EZ7",	"EZ8",	"EZ9",	"EZ10",	"EZ11",	"EZ12",	
                       "EZ13",	"EZ14",	"EZ15",	"EZ16",	"EZ17",	"EZ18",	"EZ19",	"EZ20",	"EZ21",	"EZ22",	"EZ23",	"EZ24",	
                       "EZ25",	"EZ26")

#Agregar columna de gene names
Cheng <- cbind(geneNamesCheng,Cheng)
Shchetynsky <- cbind(geneNamesShchet, uniqueShchet)
#Mo <- cbind(geneNamesMo, Control_Mo, Crohn, OligoJIA, PolyJIA, SystJIA, Ulc_Col)
Mo <- cbind(geneNamesMo, Control_Mo, SystJIA, Ulc_Col)
Mo_Crohn <- cbind(geneNamesMo, Control_Mo, Crohn)
Mo_OligoJIA <- cbind(geneNamesMo, Control_Mo, OligoJIA)
Mo_PolyJIA <- cbind(geneNamesMo, Control_Mo, PolyJIA)
Mo_SystJIA <- cbind(geneNamesMo, Control_Mo[,1:6], SystJIA)
Mo_Ulc_Col <- cbind(geneNamesMo, Control_Mo[,7:12], Ulc_Col)
ParkesDM <- cbind(geneNamesParkesDM, ParkesDM)
ParkesIBM <- cbind(geneNamesParkesIBM, ParkesIBM)
ParkesPM <- cbind(geneNamesParkesPM, ParkesPM)


#Quitando ceros
emptyCheng <- Cheng[-which(apply(Cheng[,-1], 1, mean) == 0),]
emptyShchet <- Shchetynsky[-which(apply(Shchetynsky[,-1], 1, mean) == 0),]
emptyMo <- Mo[-which(apply(Mo[,-1], 1, mean) == 0),]
emptyMo_Crohn <- Mo_Crohn[-which(apply(Mo_Crohn[,-1], 1, mean) == 0),]
emptyMo_OligoJIA <- Mo_OligoJIA[-which(apply(Mo_OligoJIA[,-1], 1, mean) == 0),]
emptyMo_PolyJIA <- Mo_PolyJIA[-which(apply(Mo_PolyJIA[,-1], 1, mean) == 0),]
emptyMo_SystJIA <- Mo_SystJIA[-which(apply(Mo_SystJIA[,-1], 1, mean) == 0),]
emptyMo_Ulc_Col <- Mo_Ulc_Col[-which(apply(Mo_Ulc_Col[,-1], 1, mean) == 0),]
#emptyParkesDM <- ParkesDM[-which(apply(ParkesDM[,-1], 1, mean) == 0),]
