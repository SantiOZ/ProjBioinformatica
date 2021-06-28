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

#PRUEBAS PARA FUCNION

  #Normalizaci?n, logaritmo y grafico individual
  normalizedCheng <- normalizeQuantiles(emptyCheng[,-1])
  logCheng <- log2(normalizedCheng + 1)
  plotDensities(logCheng, legend = "right", main = "Enfermos Lupus")
                
                
  normalizedShchet <- normalizeQuantiles(emptyShchet[,-1])
  logShchet <- log2(normalizedShchet + 1)
  plotDensities(logShchet, legend = "right", main = "Enfermos Artitris")
  
  
  #PCA individual Cheng
  madCheng <- apply(normalizedCheng, 1, mad)
  pcaCheng <- prcomp(normalizedCheng[order(madCheng, decreasing = T),])
  plot(pcaCheng$rotation[,1], pcaCheng$rotation[,2],  main = "Cheng: PC1 vs PC2", 
       xlab= paste("PCA1: ", round(summary(pcaCheng)$importance[2,1]*100,1),"%", sep=""), 
       ylab=paste("PCA2: ", round(summary(pcaCheng)$importance[2,2]*100,1),"%",sep=""),
       col= c("green", "red")) 
  
  pcaMat <- as.data.frame(t(pcaCheng$rotation)) 
  pcaMatSan <- pcaMat %>% select(contains("S"))
  pcaMatEnf <- pcaMat %>% select(contains("E"))
  rangePC1 <- range(pcaMat[,1])
  rangePC2 <- range(pcaMat[,2])
  
  
  plot(t(pcaMatSan)[,1], t(pcaMatSan)[,2],  main = "Cheng: PC1 vs PC2", 
       xlab= paste("PCA1: ", round(summary(pcaCheng)$importance[2,1]*100,1),"%", sep=""), 
       ylab=paste("PCA2: ", round(summary(pcaCheng)$importance[2,2]*100,1),"%",sep=""),
       col= c("green"),  ylim =  rangePC2 <- range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]) ) 
  par(new=TRUE)
  plot(t(pcaMatEnf)[,1], t(pcaMatEnf)[,2],  main = "Cheng: PC1 vs PC2",
       xlab= paste("PCA1: ", round(summary(pcaCheng)$importance[2,1]*100,1),"%", sep=""), 
       ylab=paste("PCA2: ", round(summary(pcaCheng)$importance[2,2]*100,1),"%",sep=""),
       col= c("red"), ylim =  rangePC2 <- range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]))
  
  #PCA individual Shchet
  madShchet <- apply(normalizedShchet, 1, mad)
  pcaShchet <- prcomp(normalizedShchet[order(madShchet, decreasing = T),])
  plot(pcaShchet$rotation[,1], pcaShchet$rotation[,2],  main = "Shchet: PC1 vs PC2", 
       xlab= paste("PCA1: ", round(summary(pcaShchet)$importance[2,1]*100,1),"%", sep=""), 
       ylab=paste("PCA2: ", round(summary(pcaShchet)$importance[2,2]*100,1),"%",sep=""),
       col= c("green", "red")) 
  
  pcaMat <- as.data.frame(t(pcaShchet$rotation)) 
  pcaMatSan <- pcaMat %>% select(contains("S"))
  pcaMatEnf <- pcaMat %>% select(contains("E"))
  rangePC1 <- range(pcaMat[,1])
  rangePC2 <- range(pcaMat[,2])
  
  
  plot(t(pcaMatSan)[,1], t(pcaMatSan)[,2],  main = "Cheng: PC1 vs PC2", 
       xlab= paste("PCA1: ", round(summary(pcaCheng)$importance[2,1]*100,1),"%", sep=""), 
       ylab=paste("PCA2: ", round(summary(pcaCheng)$importance[2,2]*100,1),"%",sep=""),
       col= c("green"),  ylim =  rangePC2 <- range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]) ) 
  par(new=TRUE)
  plot(t(pcaMatEnf)[,1], t(pcaMatEnf)[,2],  main = "Cheng: PC1 vs PC2",
       xlab= paste("PCA1: ", round(summary(pcaCheng)$importance[2,1]*100,1),"%", sep=""), 
       ylab=paste("PCA2: ", round(summary(pcaCheng)$importance[2,2]*100,1),"%",sep=""),
       col= c("red"), ylim =  rangePC2 <- range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]))
  
  #Comparacion de dos datasets
  compareData <- merge(emptyCheng[-1], emptyShchet[-1], by=0)
  rownames(compareData) <- compareData[,1]
  compareData <- compareData[,-1]
  
  #Normalizaci?n, logaritmo y grafico comparado
  normalizedCompare <- normalizeQuantiles(compareData)
  logCompare <- log2(normalizedCompare + 1)
  plotDensities(logCompare, legend = "right", main = "Enfermos Artritis & Lupus")

  
  
##FUNCIONES
  
  #Funcion global
  RunPCA <- function(x, compare = FALSE){
    if(compare == FALSE){
      if (x == "cheng"){
        par(mfrow = c(2,2))
        
        normalizedCheng <- normalizeQuantiles(emptyCheng[,-1])
        logCheng <- log2(normalizedCheng + 1)
        plotDensities(logCheng, legend = "right", main = "Enfermos Lupus")
        
        madCheng <- apply(normalizedCheng, 1, mad)
        pcaCheng <- prcomp(normalizedCheng[order(madCheng, decreasing = T),])
        plot(pcaCheng$rotation[,1], pcaCheng$rotation[,2],  main = "Cheng: PC1 vs PC2", 
             xlab= paste("PCA1: ", round(summary(pcaCheng)$importance[2,1]*100,1),"%", sep=""), 
             ylab=paste("PCA2: ", round(summary(pcaCheng)$importance[2,2]*100,1),"%",sep="")) 
        
      }
      else if (x == "shchet"){
        par(mfrow = c(2,2))
        
        normalizedShchet <- normalizeQuantiles(emptyShchet[,-1] + 1)
        logShchet <- log2(normalizedShchet)
        plotDensities(logShchet, legend = "right", main = "Enfermos Artitris")
        
        madShchet <- apply(normalizedShchet, 1, mad)
        pcaShchet <- prcomp(normalizedShchet[order(madShchet, decreasing = T),])
        plot(pcaShchet$rotation[,1], pcaShchet$rotation[,2],  main = "Artritis: PC1 vs PC2", 
             xlab= paste("PCA1: ", round(summary(pcaShchet)$importance[2,1]*100,1),"%", sep=""), 
             ylab=paste("PCA2: ", round(summary(pcaShchet)$importance[2,2]*100,1),"%",sep="")) 
      } 
      else{
        print("return monke")
      }
    }
     else if(compare == TRUE){
       par(mfrow = c(2,2))
       
       compareData <- merge(emptyCheng[-1], emptyShchet[-1], by=0)
       rownames(compareData) <- compareData[,1]
       compareData <- compareData[,-1]
       
       normalizedCompare <- normalizeQuantiles(compareData + 1)
       logCompare <- log2(normalizedCompare)
       plotDensities(logCompare, legend = "right", main = "Enfermos Artritis & Lupus")
       
       madCompare <- apply(normalizedCompare, 1, mad)
       pcaCompare <- prcomp(normalizedCompare[order(madCompare, decreasing = T),])
       plot(pcaCompare$rotation[,1], pcaCompare$rotation[,2],  main = "Artritis & Lupus: PC1 vs PC2",
            xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""), 
            ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep="")) 
     }
    }
  
  
  
  #Funcion global X
  RunPCAX <- function(x, compare = FALSE,...){
    if(compare == FALSE){
      
      par(mfrow = c(2,2))
      
      normalizedX <- normalizeQuantiles(x[,-1])
      logX <- log2(normalizedX + 1)
      plotDensities(logX, legend = "right", main = paste("Enfermos", deparse(substitute(x)), sep = " "))
      
      madX <- apply(normalizedX, 1, mad)
      pcaX <- prcomp(normalizedX[order(madX, decreasing = T),])
      
      pcaMat <- as.data.frame(t(pcaX$rotation)) 
      pcaMatSan <- pcaMat %>% select(contains("S"))
      pcaMatEnf <- pcaMat %>% select(contains("E"))
      
      plot(t(pcaMatSan)[,1], t(pcaMatSan)[,2],  main = paste(deparse(substitute(x)), ": PC1 vs PC2", sep = " "), 
           xlab= paste("PCA1: ", round(summary(pcaX)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaX)$importance[2,2]*100,1),"%",sep=""),
           col= c("green"),  ylim = range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]) ) 
      par(new=TRUE)
      plot(t(pcaMatEnf)[,1], t(pcaMatEnf)[,2],  main = paste(deparse(substitute(x)), ": PC1 vs PC2", sep = " "),
           xlab= paste("PCA1: ", round(summary(pcaX)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaX)$importance[2,2]*100,1),"%",sep=""),
           col= c("red"), ylim = range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1])) 
      
    }
    else if(compare == TRUE){
      par(mfrow = c(2,2))
      
      compareData <- merge(x[-1], ...[-1], by=0)
      rownames(compareData) <- compareData[,1]
      compareData <- compareData[,-1]
      
      normalizedCompare <- normalizeQuantiles(compareData + 1)
      logCompare <- log2(normalizedCompare)
      plotDensities(logCompare, legend = "right", main = paste("Enfermos", deparse(substitute(x)), "&", deparse(substitute(...)), sep = " "))
      
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



