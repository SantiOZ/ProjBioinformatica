library("rtracklayer")
library("dplyr")
library("limma")
library("tidyverse")
library("sva")
library("readr")
library("DESeq2")

#pipeline para analisis por separado
#-normalizar
#-log
#-plot densities
#-PCA

#datasets
#-artritis Schetynsky
#-lupus Cheng
#-Crohn Mo


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
      pcaX <- prcomp(normalizedX[order(madX, decreasing = T)[1:1000],])
      
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
      
      DatasetVec <- list()
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
      pcaCompare <- prcomp(normalizedCompare[order(madCompare, decreasing = T)[1:1000],])
      
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
  
  #Funcion global w
  RunPCAW <- function(DatasetList){
    if(length(DatasetList) == 1){
      x <- ...[[1]]
      
      par(mfrow = c(2,2))
      
      normalizedX <- getNormalizedMatrix(x[,-1])
      plotDensities(normalizedX, legend = "right", main = paste("FPKM", deparse(substitute(...[1])), sep = " "))
      
      madX <- apply(normalizedX, 1, mad)
      pcaX <- prcomp(normalizedX[order(madX, decreasing = T)[1:1000],])
      
      pcaMat <- as.data.frame(t(pcaX$rotation)) 
      pcaMatSan <- pcaMat %>% select(contains("S"))
      pcaMatEnf <- pcaMat %>% select(contains("E"))
      
      plot(t(pcaMatSan)[,1], t(pcaMatSan)[,2],  main = paste(deparse(substitute(...[1])), ": PC1 vs PC2", sep = " "), 
           xlab= paste("PCA1: ", round(summary(pcaX)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaX)$importance[2,2]*100,1),"%",sep=""),
           pch = 16,  ylim = range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1]) ) 
      par(new=TRUE)
      plot(t(pcaMatEnf)[,1], t(pcaMatEnf)[,2],  main = paste(deparse(substitute(...[1])), ": PC1 vs PC2", sep = " "),
           xlab= paste("PCA1: ", round(summary(pcaX)$importance[2,1]*100,1),"%", sep=""), 
           ylab=paste("PCA2: ", round(summary(pcaX)$importance[2,2]*100,1),"%",sep=""),
           pch = 4, ylim = range(t(pcaMat)[,2]), xlim = range(t(pcaMat)[,1])) 
    }
    
    else if(length(DatasetList) > 1){
      NuDatasets <- length(DatasetList)
      
      par(mfrow = c(2,2))
      
      DatasetVec <- list()
      for (i in 1:NuDatasets){
        NormDatasets <- getNormalizedMatrix(DatasetList[[i]][-1])
        DatasetVec[[i]] <- NormDatasets
      }

      
      compareData <- DatasetVec[[1]]
      for (i in 2:NuDatasets){
        compareData <- merge(compareData, DatasetVec[i], by=0)
        rownames(compareData) <- compareData[,1]
        compareData <- compareData[,-1]
      }
      
      
      normalizedCompare <- normalizeQuantiles(compareData)
      plotDensities(normalizedCompare, legend = "right", main = paste("FPKM", deparse(substitute(DatasetList)), sep = " "))
      
      madCompare <- apply(normalizedCompare, 1, mad)
      pcaCompare <- prcomp(normalizedCompare[order(madCompare, decreasing = T)[1:1000],])
    
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
        
        plot(t(pcaMatComSan)[,1], t(pcaMatComSan)[,2],  main = paste(deparse(substitute(DatasetList)),": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(pcaCompare)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(pcaCompare)$importance[2,2]*100,1),"%",sep=""),
             pch = 16, col= c(color[i]), ylim = range(t(pcaMatCom)[,2]), xlim = range(t(pcaMatCom)[,1]))
        par(new=TRUE)
        plot(t(pcaMatComEnf)[,1], t(pcaMatComEnf)[,2],  main = paste(deparse(substitute(DatasetList)),": PC1 vs PC2",sep = " "),
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
        
        plot(t(pcaMatComSan)[,1], t(pcaMatComSan)[,2],  main = paste(deparse(substitute(DatasetList)), ": PC1 vs PC2",sep = " "),
             xlab= paste("PCA1: ", round(summary(sc_mat_combat_pca)$importance[2,1]*100,1),"%", sep=""),
             ylab=paste("PCA2: ", round(summary(sc_mat_combat_pca)$importance[2,2]*100,1),"%",sep=""),
             pch = 16, col= c(color[i]), ylim = range(t(sc_mat_combat_mat)[,2]), xlim = range(t(sc_mat_combat_mat)[,1]))
        par(new=TRUE)
        plot(t(pcaMatComEnf)[,1], t(pcaMatComEnf)[,2],  main = paste(deparse(substitute(DatasetList)),": PC1 vs PC2",sep = " "),
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
  
  