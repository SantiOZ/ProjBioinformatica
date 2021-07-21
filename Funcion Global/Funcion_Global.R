library("rtracklayer")
library("dplyr")
library("limma")
library("tidyverse")

#pipeline para analisis por separado
#-normalizar
#-log
#-plot densities
#-PCA

#datasets
#-artritis Schetynsky
#-lupus Cheng

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
  
  
  #Opening index 38.98 GTF file
  index <- as.data.frame(rtracklayer::import('C:/Users/tt_kb/Documents/Bioinformatica/Index.gtf'))
  
  #Matching datasets gene id with Index
  Cheng_match <- match(rownames(Cheng), index$gene_id)
  Cheng_match_clean <- Cheng_match[!is.na(Cheng_match)]
  
  Shchet_match <- match(rownames(uniqueShchet), index$gene_id)
  Shchet_match_clean <- Shchet_match[!is.na(Shchet_match)]
  
  #Obtener nombres a partir de gene id
  geneNamesCheng <- index$gene_name[Cheng_match_clean]
  geneNamesShchet <- index$gene_name[Shchet_match_clean]
  
  #Eliminar filas sin match con index en dataset original
  Cheng <- Cheng[!is.na(Cheng_match),]
  uniqueShchet <- uniqueShchet[!is.na(Shchet_match),]
  
  #Agregar columna de gene names
  Cheng <- cbind(geneNamesCheng,Cheng)
  Shchetynsky <- cbind(geneNamesShchet, uniqueShchet)
  
  #Quitando ceros
  emptyCheng <- Cheng[-which(apply(Cheng[,-1], 1, mean) == 0),]
  emptyShchet <- Shchetynsky[-which(apply(Cheng[,-1], 1, mean) == 0),]

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
  
  
  plot(t(pcaMatSan)[,1], t(pcaMatSan)[,2],  main = "Cheng: PC1 vs PC2", 
       xlab= paste("PCA1: ", round(summary(pcaCheng)$importance[2,1]*100,1),"%", sep=""), 
       ylab=paste("PCA2: ", round(summary(pcaCheng)$importance[2,2]*100,1),"%",sep=""),
       pch = 1,  ylim =  rangePC2 <- range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]) ) 
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
  pcaCxS <- prcomp(normalizedCxS[order(madCxS, decreasing = T),]) #Datasets invertidos
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
  
  #Funcion global X
  RunPCAX <- function(x, compare = FALSE,...){
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
      #NuDataSets <- length(list(x, ...))

      par(mfrow = c(2,2))
      
      compareData <- merge(x[-1], ...[-1], by=0)
      rownames(compareData) <- compareData[,1]
      compareData <- compareData[,-1]
      
      
      logCompare <- log2(compareData + 1)
      normalizedCompare <- normalizeQuantiles(logCompare)
      plotDensities(normalizedCompare, legend = "right", main = paste("FPKM", deparse(substitute(x)), "&", deparse(substitute(...)), sep = " "))
      
      madCompare <- apply(normalizedCompare, 1, mad)
      pcaCompare <- prcomp(normalizedCompare[order(madCompare, decreasing = T),])
      
      pcaMatCom <- as.data.frame(t(pcaCompare$rotation)) 
      pcaMatComSan <- pcaMatCom %>% select(contains("S"))
      pcaMatComEnf <- pcaMatCom %>% select(contains("E"))
      
      plot(t(pcaMatComSan)[,1], t(pcaMatComSan)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
           xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep=""),
           col= c("green"), ylim = range(t(pcaMatCom)[,2]), xlim = range(t(pcaMatCom)[,1])) 
      par(new=TRUE)
      plot(t(pcaMatComEnf)[,1], t(pcaMatComEnf)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
           xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep=""),
           col= c("red"), ylim = range(t(pcaMatCom)[,2]), xlim = range(t(pcaMatCom)[,1])) 
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
      for(i in 1:NuDatasets){
        NuCol <- ncol(DatasetVec[[i]])
        final <- NuCol + inicio - 1
        
        pcaMatComSan <- pcaMatCom[inicio:final] %>% select(contains("S"))
        pcaMatComEnf <- pcaMatCom[inicio:final] %>% select(contains("E"))
        
        plot(t(pcaMatComSan)[,1], t(pcaMatComSan)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep=""),
             pch = 1, col= c(color[i]), ylim = range(t(pcaMatCom)[,2]), xlim = range(t(pcaMatCom)[,1]))
        par(new=TRUE)
        plot(t(pcaMatComEnf)[,1], t(pcaMatComEnf)[,2],  main = paste(deparse(substitute(x)), "&", deparse(substitute(...)),": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep=""),
             pch = 8, col= c(color[i]), ylim = range(t(pcaMatCom)[,2]), xlim = range(t(pcaMatCom)[,1]))
        par(new=TRUE)
        
        inicio <- NuCol + 1
      }
    }
  }
