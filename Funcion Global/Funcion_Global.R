library("rtracklayer")
library("dplyr")
library("limma")
library("tidyverse")
library("sva")
library("readr")

#pipeline para analisis por separado
#-normalizar
#-log
#-plot densities
#-PCA

#datasets
#-artritis Schetynsky
#-lupus Cheng
#-Crohn Mo

#PREPARACION DE DATASETS FPKM

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

#PRUEBAS PARA FUNCION
  
  #Grafico de densidades Cheng y Shchet
  logCheng <- log2(emptyCheng[,-1] + 1)
  normalizedCheng <- normalizeQuantiles(logCheng)
  plotDensities(normalizedCheng, legend = "right", main = "Enfermos Lupus")
                
  normalizedShchet <- normalizeQuantiles(emptyShchet[,-1])
  logShchet <- log2(normalizedShchet + 1)
  plotDensities(logShchet, legend = "right", main = "Enfermos Artitris")
  
  
  #PCA Cheng: Diferenciacion de columnas por figura
  madCheng <- apply(normalizedCheng, 1, mad)
  pcaCheng <- prcomp(normalizedCheng[order(madCheng, decreasing = T),])
 
  pcaMat <- as.data.frame(t(pcaCheng$rotation)) 
  pcaMatSan <- pcaMat %>% select(contains("S"))
  pcaMatEnf <- pcaMat %>% select(contains("E"))
  rangePC1 <- range(pcaMat[,1])
  rangePC2 <- range(pcaMat[,2])
  DT_Names<-colnames(pcaMatSan[1])
  print(substr(DT_Names, 1, 2))
  print(substr(colnames(pcaMatSan[1]), 1, 2))
  colnames(model$data[[1]])
  
  plot(t(pcaMatSan)[,1], t(pcaMatSan)[,2],  main = "Cheng: PC1 vs PC2", 
       xlab= paste("PCA1: ", round(summary(pcaCheng)$importance[2,1]*100,1),"%", sep=""), 
       ylab=paste("PCA2: ", round(summary(pcaCheng)$importance[2,2]*100,1),"%",sep=""),
       pch = 1,  ylim =  rangePC2 <- range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]) ) 
  par(xpd=TRUE)
  legend("bottomright", 
         legend = c(substr(colnames(pcaMatSan[1]), 1, 2), substr(colnames(pcaMatEnf[1]), 1, 2)), 
         col = c("red", 
                 "red"), 
         pch = c(16,4), 
         bty = "n", 
         pt.cex = 2, 
         cex = 1.2, 
         text.col = "black", 
         horiz = F , 
         inset = c(0, 0))
  par(new=TRUE)
  plot(t(pcaMatEnf)[,1], t(pcaMatEnf)[,2],  main = "Cheng: PC1 vs PC2",
       xlab= paste("PCA1: ", round(summary(pcaCheng)$importance[2,1]*100,1),"%", sep=""), 
       ylab=paste("PCA2: ", round(summary(pcaCheng)$importance[2,2]*100,1),"%",sep=""),
       pch = 8, ylim =  rangePC2 <- range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]))

  
  #Unión de dos datasets para comparación Cheng y Shchet
  compareCxS <- merge(emptyCheng[-1], emptyShchet[-1], by=0)
  rownames(compareCxS) <- compareCxS[,1]
  compareCxS <- compareCxS[,-1]
  
  
  #ChengxShchet: Normalizacion, logaritmo y grafico
  logCxS <- log2(compareCxS + 1)
  normalizedCxS <- normalizeQuantiles(logCxS)
  plotDensities(normalizedCxS, legend = "right", main = "Enfermos Artritis & Lupus")
  
  
  #ChengxShchet: PCA diferenciado por figura y color
  
  color <- c("red", "green", "cyan", "blue", "purple", "magenta", "yellow")
  inicio <- 1
  DataSets <- list(emptyCheng[-1], emptyShchet[-1])
  NuDataSets <- length(DataSets)
  
  madCxS <- apply(normalizedCxS, 1, mad)
  pcaCxS <- prcomp(normalizedCxS[order(madCxS, decreasing = T)[1:1000],]) #Datasets invertidos
  pcaMatCxS <- as.data.frame(t(pcaCxS$rotation)) 
  
  
  for(i in 1:NuDataSets){
    print(color[i])
    NuCol <- ncol(DataSets[[i]])
    print(NuCol)
    final <- NuCol + inicio
    
    pcaMatCxS_San <- pcaMatCxS[inicio:final] %>% select(contains("S"))
    pcaMatCxS_Enf <- pcaMatCxS[inicio:final] %>% select(contains("E"))
    
    plot(t(pcaMatCxS_San)[,1], t(pcaMatCxS_San)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
         xlab= paste("PCA1: ", round(summary(pcaCxS)$importance[2,1]*100,1),"%", sep=""), 
         ylab=paste("PCA2: ", round(summary(pcaCxS)$importance[2,2]*100,1),"%",sep=""),
         pch = 1, col= c(color[i]), ylim = range(t(pcaMatCxS)[,2]), xlim = range(t(pcaMatCxS)[,1])) 
    par(new=TRUE)
    plot(t(pcaMatCxS_Enf)[,1], t(pcaMatCxS_Enf)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
         xlab= paste("PCA1: ", round(summary(pcaCxS)$importance[2,1]*100,1),"%", sep=""), 
         ylab=paste("PCA2: ", round(summary(pcaCxS)$importance[2,2]*100,1),"%",sep=""),
         pch = 8, col= c(color[i]),  ylim = range(t(pcaMatCxS)[,2]), xlim = range(t(pcaMatCxS)[,1])) 
    par(new=TRUE)
    
    inicio <- NuCol
  }
  
  
  
  
##FUNCIONES
  
  #Funcion para normalizar dependiendo de conteos o fpkm
  getNormalizedMatrix <- function(count_mat){
    if(all(apply(count_mat,2,is.integer), na.rm = FALSE) ){
      samples_classA = count_mat %>% select(contains("S"))
      samples_classB = count_mat %>% select(contains("E"))
      aux_data <- count_mat
      aux_desc <- data.frame(condition=c(rep("Sanos",length(samples_classA)),rep("Enfermos",length(samples_classB))), type=rep("paired-end",c(length(samples_classA)+length(samples_classB))))
      aux_dds <- DESeqDataSetFromMatrix(countData = aux_data, colData = aux_desc, design = ~condition)
      aux_dds <- DESeq(aux_dds)
      normalized_counts <- counts(aux_dds, normalized=T)
      log_counts <- log2(normalized_counts + 1)
      return(log_counts)
    }
    else{
      logX <- log2(count_mat + 1)
      normalizeQuantiles(logX)
    }
  }
  


  #Funcion global y
  RunPCAY <- function(x, compare = FALSE,...){
    if(compare == FALSE){
      
      par(mfrow = c(2,2))
      
      logX <- log2(x[,-1] + 1)
      normalizedX <- normalizeQuantiles(logX)
      plotDensities(normalizedX, legend = "right", main = paste("FPKM", deparse(substitute(x)), sep = " "))
      
      madX <- apply(normalizedX, 1, mad)
      pcaX <- prcomp(normalizedX[order(madX, decreasing = T),])
      
      pcaMat <- as.data.frame(t(pcaX$rotation)) 
      pcaMatSan <- pcaMat %>% select(contains("S"))
      pcaMatEnf <- pcaMat %>% select(contains("E"))
      
      plot(t(pcaMatSan)[,1], t(pcaMatSan)[,2],  main = paste(deparse(substitute(x)), ": PC1 vs PC2", sep = " "), 
           xlab= paste("PCA1: ", round(summary(pcaX)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaX)$importance[2,2]*100,1),"%",sep=""),
           pch = 1,  ylim = range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]) ) 
      par(new=TRUE)
      plot(t(pcaMatEnf)[,1], t(pcaMatEnf)[,2],  main = paste(deparse(substitute(x)), ": PC1 vs PC2", sep = " "),
           xlab= paste("PCA1: ", round(summary(pcaX)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaX)$importance[2,2]*100,1),"%",sep=""),
           pch = 8, ylim = range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1])) 
    }
    
    else if(compare == TRUE){
      NuDatasets <- length(list(x, ...))
      
      par(mfrow = c(2,2))
      
      DatasetVec <- vector()
      for (i in 1:NuDatasets){
        Datasets <- list(x, ...)
        NormDatasets <- list(normalizeQuantiles(Datasets[[i]][-1]))
        DatasetVec <- c(DatasetVec, NormDatasets[1])
      }
      
      compareData <- DatasetVec[1]
      for (i in 2:NuDatasets){
        compareData <- merge(compareData,DatasetVec[i], by=0)
        rownames(compareData) <- compareData[,1]
        compareData <- compareData[,-1]
      }
      
      logCompare <- log2(compareData + 1)
      normalizedCompare <- normalizeQuantiles(logCompare)
      plotDensities(normalizedCompare, legend = "right", main = paste("FPKM", deparse(substitute(x)), "&", deparse(substitute(...)), sep = " "))
      
      madCompare <- apply(normalizedCompare, 1, mad)
      pcaCompare <- prcomp(normalizedCompare[order(madCompare, decreasing = T),])
      
      pcaMatCom <- as.data.frame(t(pcaCompare$rotation))
      
      color <- c("red", "green", "cyan", "blue", "purple", "magenta", "yellow")
      inicio <- 1
      insety <- 0.4
      for(i in 1:NuDatasets){
        NuCol <- ncol(DatasetVec[[i]])
        final <- NuCol + inicio - 1
        
        pcaMatComSan <- pcaMatCom[inicio:final] %>% select(contains("S"))
        pcaMatComEnf <- pcaMatCom[inicio:final] %>% select(contains("E"))
        
        plot(t(pcaMatComSan)[,1], t(pcaMatComSan)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep=""),
             pch = 16, col= c(color[i]), ylim = range(t(pcaMatCom)[,2]), xlim = range(t(pcaMatCom)[,1]))
        par(new=TRUE)
        plot(t(pcaMatComEnf)[,1], t(pcaMatComEnf)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep=""),
             pch = 4, col= c(color[i]), ylim = range(t(pcaMatCom)[,2]), xlim = range(t(pcaMatCom)[,1]))
        
        legend("bottomright",
               legend = c("SL", "EL"), 
               col = c(color[i], 
                       color[i]), 
               pch = c(16,4), 
               bty = "n", 
               pt.cex = 1, 
               cex = 1, 
               text.col = "black", 
               horiz = F,
               inset = c(0.0, insety))
        insety <- insety - 0.4
        par(new=TRUE)
        
        inicio <- final + 1
      }
    }
  }

  
  #Funcion global z
  RunPCAZ <- function(x, compare = FALSE,...){
    if(compare == FALSE){
      
      par(mfrow = c(2,2))
      
      normalizedX <- getNormalizedMatrix(x[,-1])
      plotDensities(normalizedX, legend = "right", main = paste("FPKM", deparse(substitute(x)), sep = " "))
      
      madX <- apply(normalizedX, 1, mad)
      pcaX <- prcomp(normalizedX[order(madX, decreasing = T),])
      
      pcaMat <- as.data.frame(t(pcaX$rotation)) 
      pcaMatSan <- pcaMat %>% select(contains("S"))
      pcaMatEnf <- pcaMat %>% select(contains("E"))
      
      plot(t(pcaMatSan)[,1], t(pcaMatSan)[,2],  main = paste(deparse(substitute(x)), ": PC1 vs PC2", sep = " "), 
           xlab= paste("PCA1: ", round(summary(pcaX)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaX)$importance[2,2]*100,1),"%",sep=""),
           pch = 16,  ylim = range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]) ) 
      par(new=TRUE)
      plot(t(pcaMatEnf)[,1], t(pcaMatEnf)[,2],  main = paste(deparse(substitute(x)), ": PC1 vs PC2", sep = " "),
           xlab= paste("PCA1: ", round(summary(pcaX)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaX)$importance[2,2]*100,1),"%",sep=""),
           pch = 4, ylim = range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1])) 
    }
    
    else if(compare == TRUE){
      NuDatasets <- length(list(x, ...))
      
      par(mfrow = c(2,2))
      
      DatasetVec <- vector()
      for (i in 1:NuDatasets){
        Datasets <- list(x, ...)
        NormDatasets <- list(getNormalizedMatrix(Datasets[[i]][-1]))
        DatasetVec <- c(DatasetVec, NormDatasets[1])
      }
      
      compareData <- DatasetVec[1]
      for (i in 2:NuDatasets){
        compareData <- merge(compareData,DatasetVec[i], by=0)
        rownames(compareData) <- compareData[,1]
        compareData <- compareData[,-1]
      }
      
      normalizedCompare <- normalizeQuantiles(compareData)
      plotDensities(normalizedCompare, legend = "right", main = paste("FPKM", deparse(substitute(x)), "&", deparse(substitute(...)), sep = " "))
      
      madCompare <- apply(normalizedCompare, 1, mad)
      pcaCompare <- prcomp(normalizedCompare[order(madCompare, decreasing = T),])
      
      pcaMatCom <- as.data.frame(t(pcaCompare$rotation))
      
      ###
      inicio <- 1
      sc_mat_batch <- c()
      for (i in 1:NuDatasets){
        NuCol <- ncol(DatasetVec[[i]])
        final <- NuCol + inicio - 1
        
        cutData <- normalizedCompare[inicio:final]
        
        sc_mat_batch <- c(sc_mat_batch, rep(i, ncol(cutData)))
        
        inicio <- NuCol + 1
      }
      
      
      color <- c("red", "green", "cyan", "blue", "purple", "magenta", "yellow")
      inicio <- 1
      insety <- 0.7
      for(i in 1:NuDatasets){
        NuCol <- ncol(DatasetVec[[i]])
        final <- NuCol + inicio - 1
        
        pcaMatComSan <- pcaMatCom[inicio:final] %>% select(contains("S"))
        pcaMatComEnf <- pcaMatCom[inicio:final] %>% select(contains("E"))
        
        plot(t(pcaMatComSan)[,1], t(pcaMatComSan)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep=""),
             pch = 16, col= c(color[i]), ylim = range(t(pcaMatCom)[,2]), xlim = range(t(pcaMatCom)[,1]))
        par(new=TRUE)
        plot(t(pcaMatComEnf)[,1], t(pcaMatComEnf)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep=""),
             pch = 4, col= c(color[i]), ylim = range(t(pcaMatCom)[,2]), xlim = range(t(pcaMatCom)[,1]))
        par(xpd=TRUE)
        legend("bottomright",
               legend = c(substr(colnames(pcaMatComSan[1]), 1, 2),substr(colnames(pcaMatComEnf[1]), 1, 2)), 
               col = c(color[i], 
                       color[i]), 
               pch = c(16,4), 
               bty = "n", 
               pt.cex = 0.8, 
               cex = 0.8, 
               text.col = "black", 
               horiz = F,
               inset = c(-0.12, insety))
        insety <- insety - 0.3
        par(new=TRUE)
        
        inicio <- final + 1
      }
      
      
      
      #Se introduce vector que define bathches
      modcombat <- model.matrix(~1, data=data.frame(cmb=sc_mat_batch))
      sc_mat_combat <- ComBat(dat=normalizedCompare, batch=sc_mat_batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
      sc_mat_combat_mad <- apply(sc_mat_combat, 1, mad)
      sc_mat_combat_pca <- prcomp(sc_mat_combat[order(sc_mat_combat_mad, decreasing = TRUE)[1:1000],])
      sc_mat_combat_mat<- as.data.frame(t(sc_mat_combat_pca$rotation)) 
      ###
      
      par(new=FALSE)
      inicio <- 1
      insety <- 0.7
      for(i in 1:NuDatasets){
        NuCol <- ncol(DatasetVec[[i]])
        final <- NuCol + inicio - 1
        
        pcaMatComSan <- sc_mat_combat_mat[inicio:final] %>% select(contains("S"))
        pcaMatComEnf <- sc_mat_combat_mat[inicio:final] %>% select(contains("E"))
        
        plot(t(pcaMatComSan)[,1], t(pcaMatComSan)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(sc_mat_combat_pca)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(sc_mat_combat_pca)$importance[2,2]*100,1),"%",sep=""),
             pch = 16, col= c(color[i]), ylim = range(t(sc_mat_combat_mat)[,2]), xlim = range(t(sc_mat_combat_mat)[,1]))
        par(new=TRUE)
        plot(t(pcaMatComEnf)[,1], t(pcaMatComEnf)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(sc_mat_combat_pca)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(sc_mat_combat_pca)$importance[2,2]*100,1),"%",sep=""),
             pch = 4, col= c(color[i]), ylim = range(t(sc_mat_combat_mat)[,2]), xlim = range(t(sc_mat_combat_mat)[,1]))
        #par(xpd=TRUE)
        legend("bottomright",
               legend = c(substr(colnames(pcaMatComSan[1]), 1, 2),substr(colnames(pcaMatComEnf[1]), 1, 2)), 
               col = c(color[i], 
                       color[i]), 
               pch = c(16,4), 
               bty = "n", 
               pt.cex = 0.8, 
               cex = 0.8, 
               text.col = "black", 
               horiz = F,
               inset = c(-0.12, insety))
        insety <- insety - 0.3
        par(new=TRUE)
        
        inicio <- final + 1
      }
    }
  }
  