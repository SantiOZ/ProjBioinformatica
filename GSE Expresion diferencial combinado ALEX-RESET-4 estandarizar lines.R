# ver 0.2 : Se le agrego el analisis de pca

{if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
  BiocManager::install("tximport")
}

{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("DESeq2")
}

{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("rtracklayer")
}

{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("limma")
}

{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("sva")
}

install.packages("doBy")
install.packages("gplots")



writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
Sys.which("make")
"C:\\rtools40\\usr\\bin\\make.exe"
install.packages("jsonlite", type = "source")
library(jsonlite)


install.packages("C:/Users/aleja/Documents/R/refGenome", repos = NULL, type="source")


install.packages("RSQLite")

install.packages("gplots")
install.packages("ActivePathways")
install.packages("rtacklayer")
install.packages("dplyr")
install.packages("tibble")
install.packages("scatterplot3d")


{
  library(tximport)
  library(DESeq2)
  library(doBy)
  library(RSQLite)
  library(gplots)
  library(ActivePathways)
  library(rtracklayer)
  library(dplyr)
  library(tibble)
  library(scatterplot3d)
  library(limma)
  library(sva)
  library(refGenome)
}


setwd("D:/Users/Disco Duro Intero de Alex/Documents/Bioinformatics Team/2020 LUPUS/Analisis en R")#Se accede a un directorio de trabajo donde se encuentran los archivos que se utilizaran

archivos_tsv = list.files(pattern="*.tsv")#Lee todos los archivos que sigan el patron indicado
abundancias <- list()
for (i in 1:length(archivos_tsv)) {#Para cada archivo encontrado en abundancias, le aplicara la funcion read.table, ademas de cortar el nombre para quitar ".tsv"
  if (nchar(archivos_tsv[i])== 7) {
    abundancias[[i]] <- assign(substr(archivos_tsv[i],1,3),read.table(archivos_tsv[i], sep = "\t", header=TRUE))
    names(abundancias)[i] <- c(substr(archivos_tsv[i],1,3))
  } else {
    abundancias[[i]] <- assign(substr(archivos_tsv[i],1,4),read.table(archivos_tsv[i], sep = "\t", header=TRUE))
    names(abundancias)[i] <- c(substr(archivos_tsv[i],1,17))
  }
}

automatedDESeq2=function(count_mat=count_mat, comparisons=comparisons){
  # Si no hacemos esta l?nea, y hay genes con ceros, sale pval/qval como NA
  count_mat = round(count_mat +1)
  #Renglones -> Genes, Columnas -> Replicas
  resList=list()
  for (n in 1:length(comparisons)){
    controls = comparisons[[n]][1] #Define un control a partir de los objetos de la lista comparisons
    samples = comparisons[[n]][2] #Define las muestras a partir de los objetos de la lista comparisons
    samples_classA = colnames(count_mat)[grep(controls, colnames(count_mat))] #Toma todas las replicas correspondientes al control
    samples_classB = colnames(count_mat)[grep(samples, colnames(count_mat))] #Toma todas las replicas correspondientes al las muestras
    aux_data <- count_mat[,c(samples_classA,samples_classB)] #Concatena los estimated counts de control y muestra
    aux_desc <- data.frame(condition=c(rep("treated",length(samples_classA)),rep("untreated",length(samples_classB))), type=rep("paired-end",c(length(samples_classA)+length(samples_classB)))) #Clasifica controles con la condicion treated y muestras como untreated
    aux_dds <- DESeqDataSetFromMatrix(countData = aux_data, colData = aux_desc, design = ~condition) #Corre la funcion DESeqDataSetFromMatrix con los objetos aux_data, aux_desc y condition
    aux_dds <- DESeq(aux_dds) #Corre la funcion DESeq, la cual realiza los estadisticos
    aux_res <- as.data.frame(results(aux_dds)) #Guarda en forma de data frame el resultado de la funcion DESeq
    # ponemos validaciones para que no salgan NAs
    aux_res$log2FoldChange = sapply(aux_res$log2FoldChange, function(x) ifelse(is.na(x), 0, x)) #Si el valor de x es NA, lo convierte en 0, de lo contrario queda con el valor original de x
    aux_res$pvalue = sapply(aux_res$pvalue, function(x) ifelse(is.na(x), 1, x)) #Si el valor de x es NA, lo convierte en 1, de lo contrario queda con el valor original de x
    aux_res$padj = sapply(aux_res$padj, function(x) ifelse(is.na(x), 1, x)) #Si el valor de x es NA, lo convierte en 1, de lo contrario queda con el valor original de x
    
    resList[[paste(controls,"_",samples,sep="")]] = cbind(aux_res, counts(aux_dds, normalized=T))
    # colnames(resList[[paste(controls,"_",samples,sep="")]]) = c("baseMean","log2FoldChange","lfcSE","stat", "pvalue", "padj", paste(controls[n],"_1",sep=""), paste(controls[n],"_2",sep=""), paste(samples[n],"_1",sep=""), paste(samples[n],"_2",sep=""))	
    colnames(resList[[paste(controls,"_",samples,sep="")]]) = c("baseMean","log2FoldChange","lfcSE","stat", "pvalue", "padj", samples_classA, samples_classB)	
  }
  return(resList) #Devuelve el data frame final que contiene los estadisticos y los conteos de control y muestra
}

comparisons <- list(c("C", "L"))#Un argumento necesario para la funcion DESeq2, el cual indica entre que conjunto de datos realizar las comparaciones

count_mat <- c()#Inicializa un vector vacio para construir una matriz en el for loop
for (i in 1:length(abundancias)) {#Concatena las columnas llamadas est_counts con el valor anterior de count_mat
  count_mat <- cbind(count_mat, abundancias[[i]]$est_counts)  
}
rownames(count_mat) <- abundancias[[1]]$target_id
colnames(count_mat) <- names(abundancias)
#count_mat1 <- round(count_mat + 1)

resMat <- automatedDESeq2(count_mat, comparisons)
genes <- as.data.frame(rtracklayer::import('D:/Users/Disco Duro Intero de Alex/Documents/Bioinformatics Team/2020 LUPUS/Analisis en R/gencode.v34.chr_patch_hapl_scaff.annotation.gtf.gz'))
ens <- ensemblGenome()

aux_tpm_S_E <- c()
for (i in 1:length(abundancias)) {#Genera variables que contienen las columnas de tpm 
  aux_tpm_S_E <- cbind(aux_tpm_S_E, abundancias[[i]]$tpm)
}
rownames(aux_tpm_S_E) <- abundancias$LN_TBR_SRR6730152$target_id
colnames(aux_tpm_S_E) <- names(abundancias)

aux_counts_S_E <- c()
for (i in 1:length(abundancias)) {#Genera variables que contienen las columnas de tpm 
  aux_counts_S_E <- cbind(aux_counts_S_E, abundancias[[i]]$est_counts)
}
rownames(aux_counts_S_E) <- abundancias$LN_TBR_SRR6730152$target_id
colnames(aux_counts_S_E) <- names(abundancias)


#####################################################################################################################1
#####################################################################################################################1
#####################################################################################################################1
#----------------------------------------------------------------PCA y comBat de los datos de los autores------------------------------------------------------

{
gz_paper <- file('C:/Users/aleja/Documents/Bioinformatics team/Analisis en R/GSE110685_PCarlucci_RNA-Seq_tc_rpkm_combined_submitted.txt', 'rt')
ftp_paper <- read.table(gz_paper,header=T, sep = '\t')
libs <- c('S01', 'S03', 'S14', 'S15', 'S16', 'S17', 'S19', 'S21', 'S34', 'S35', 'S36', 'S37', 'S46', 'S47', 'S51', 'S52', 'S53', 'S02', 'S06', 'S08', 'S10', 'S11', 'S13', 'S20', 'S23', 'S24')
etiq <- c('SAT1', 'SAT2', 'SAT3', 'SAT4', 'SAT5', 'SAT6', 'SAT7', 'SAT8', 'SAD1', 'SAD2', 'SAD3', 'SAD4', 'SAD5', 'SAD6', 'SAD7', 'SAD8', 'SAD9', 'EL1', 'EL2', 'EL3', 'EL4', 'EL5', 'EL6', 'EL7', 'EL8', 'EL9')
paper_df <- ftp_paper[,libs]
colnames(paper_df) <- etiq
paper_df <- log2(empty_row(paper_df + 1,2))
#Hacer normalizaci�n por cuantiles
paper_df<- normalizeQuantiles(paper_df)
plotDensities(paper_df, legend = "topright")
plot.densities(paper_df)
#plotDensities(paper_df)

paper_mad <- apply(paper_df, 1, mad)
paper_pca <- prcomp(paper_df[order(paper_mad, decreasing = T)[0:100],])


#jpeg(file="Paper PreCombat 1 v.s. 2.jpg",width = 500, height = 500)
plot(paper_pca$rotation[,1], paper_pca$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC2", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)
#dev.off()

#jpeg(file="Paper PreCombat 1 v.s. 3.jpg",width = 500, height = 500)
plot(paper_pca$rotation[,1], paper_pca$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)
#dev.off()

#jpeg(file="Paper PreCombat 2 v.s. 3.jpg",width = 500, height = 500)
plot(paper_pca$rotation[,2], paper_pca$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC2 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 3)
#dev.off()

paper_mat_batch <- c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1)
paper_modcombat <- model.matrix(~1, data=data.frame(cmb=paper_mat_batch))
paper_postComBat <- ComBat(dat= paper_df, batch=paper_mat_batch, mod=paper_modcombat, par.prior=TRUE, prior.plots=FALSE)

paper_mad_postcombat <- apply(paper_postComBat, 1, mad)
paper_pca_postcombat <- prcomp(paper_postComBat[order(paper_mad_postcombat, decreasing = T)[0:1000],])


max(paper_postComBat)
min(paper_postComBat)

#jpeg(file="Paper PostCombat 1 v.s. 2.jpg",width = 500, height = 500)
plot(paper_pca_postcombat$rotation[,1], paper_pca_postcombat$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Post ComBat: PC1 vs PC2", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)
#dev.off()

#jpeg(file="Paper PostCombat 1 v.s. 3.jpg",width = 500, height = 500)
plot(paper_pca_postcombat$rotation[,1], paper_pca_postcombat$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Post ComBat: PC1 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)
#dev.off()

#jpeg(file="Paper PostCombat 2 v.s. 3.jpg",width = 500, height = 500)
plot(paper_pca_postcombat$rotation[,2], paper_pca_postcombat$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Post ComBat: PC2 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 4)
#dev.off()
}

#--------------------------------------------------------------------------------------------------------------------------------------------------
#####################################################################################################################2
#####################################################################################################################2
#####################################################################################################################2
#----------------------------------------------------------------PCA y comBat de KALLISTO------------------------------------------------------

# Recolecta los TSV�s para tener los data sets  
{
{
  setwd("C:/Users/aleja/Documents/Bioinformatics team/Analisis en R/")
  
  
  Ct1<-read.table(file="CTRL_SRR6730133_abundance.tsv", sep = "\t", header=TRUE)
  Ct2<-read.table(file="CTRL_SRR6730135_abundance.tsv", sep = "\t", header=TRUE)
  Ct3<-read.table(file="CTRL_SRR6730146_abundance.tsv", sep = "\t", header=TRUE)
  Ct4<-read.table(file="CTRL_SRR6730147_abundance.tsv", sep = "\t", header=TRUE)
  Ct5<-read.table(file="CTRL_SRR6730148_abundance.tsv", sep = "\t", header=TRUE)
  Ct6<-read.table(file="CTRL_SRR6730149_abundance.tsv", sep = "\t", header=TRUE)
  Ct7<-read.table(file="CTRL_SRR6730151_abundance.tsv", sep = "\t", header=TRUE)
  Ct8<-read.table(file="CTRL_SRR6730153_abundance.tsv", sep = "\t", header=TRUE)
  Ct9<-read.table(file="CTRL_SRR6730166_abundance.tsv", sep = "\t", header=TRUE)
  Ct10<-read.table(file="CTRL_SRR6730167_abundance.tsv", sep = "\t", header=TRUE)
  Ct11<-read.table(file="CTRL_SRR6730168_abundance.tsv", sep = "\t", header=TRUE)
  Ct12<-read.table(file="CTRL_SRR6730169_abundance.tsv", sep = "\t", header=TRUE)
  Ct13<-read.table(file="CTRL_SRR6730178_abundance.tsv", sep = "\t", header=TRUE)
  Ct14<-read.table(file="CTRL_SRR6730179_abundance.tsv", sep = "\t", header=TRUE)
  Ct15<-read.table(file="CTRL_SRR6730183_abundance.tsv", sep = "\t", header=TRUE)
  Ct16<-read.table(file="CTRL_SRR6730184_abundance.tsv", sep = "\t", header=TRUE)
  Ct17<-read.table(file="CTRL_SRR6730185_abundance.tsv", sep = "\t", header=TRUE)
  
  
  Ln1<-read.table(file="LN_TBR_SRR6730134_abundance.tsv", sep = "\t", header=TRUE)
  Ln2<-read.table(file="LN_TBR_SRR6730138_abundance.tsv", sep = "\t", header=TRUE)
  Ln3<-read.table(file="LN_TBR_SRR6730140_abundance.tsv", sep = "\t", header=TRUE)
  Ln4<-read.table(file="LN_TBR_SRR6730142_abundance.tsv", sep = "\t", header=TRUE)
  Ln5<-read.table(file="LN_TBR_SRR6730143_abundance.tsv", sep = "\t", header=TRUE)
  Ln6<-read.table(file="LN_TBR_SRR6730145_abundance.tsv", sep = "\t", header=TRUE)
  Ln7<-read.table(file="LN_TBR_SRR6730152_abundance.tsv", sep = "\t", header=TRUE)
  Ln8<-read.table(file="LN_TBR_SRR6730155_abundance.tsv", sep = "\t", header=TRUE)
  Ln9<-read.table(file="LN_TBR_SRR6730156_abundance.tsv", sep = "\t", header=TRUE)
  
  
  n= 1
}

# Recolecta la base del gtf
# Aparte crea valores de los tpm de cada muestra
{
  genes <- as.data.frame(rtracklayer::import('gencode.v30.chr_patch_hapl_scaff.annotation.gtf.gz'))
  tpm_Ct1 <- c(Ct1$tpm)
  tpm_Ct2 <- c(Ct2$tpm)
  tpm_Ct3 <- c(Ct3$tpm)
  tpm_Ct4 <- c(Ct4$tpm)
  tpm_Ct5 <- c(Ct5$tpm)
  tpm_Ct6 <- c(Ct6$tpm)
  tpm_Ct7 <- c(Ct7$tpm)
  tpm_Ct8 <- c(Ct8$tpm)
  tpm_Ct9 <- c(Ct9$tpm)
  tpm_Ct10 <- c(Ct10$tpm)
  tpm_Ct11 <- c(Ct11$tpm)
  tpm_Ct12 <- c(Ct12$tpm)
  tpm_Ct13 <- c(Ct13$tpm)
  tpm_Ct14 <- c(Ct14$tpm)
  tpm_Ct15 <- c(Ct15$tpm)
  tpm_Ct16 <- c(Ct16$tpm)
  tpm_Ct17 <- c(Ct17$tpm)
  
  
  tpm_Ln1 <- c(Ln1$tpm)
  tpm_Ln2 <- c(Ln2$tpm)
  tpm_Ln3 <- c(Ln3$tpm)
  tpm_Ln4 <- c(Ln4$tpm)
  tpm_Ln5 <- c(Ln5$tpm)
  tpm_Ln6 <- c(Ln6$tpm)
  tpm_Ln7 <- c(Ln7$tpm)
  tpm_Ln8 <- c(Ln8$tpm)
  tpm_Ln9 <- c(Ln9$tpm)

  targid <- c(Ct1$target_id)
  
  }



ftp_kallisto <- data.frame(targid ,tpm_Ct1, tpm_Ct2, tpm_Ct3, tpm_Ct4, tpm_Ct5, tpm_Ct6, tpm_Ct7, tpm_Ct8, tpm_Ct9, tpm_Ct10, tpm_Ct11, tpm_Ct12, tpm_Ct13, tpm_Ct14, tpm_Ct15, tpm_Ct16, tpm_Ct17, tpm_Ln1, tpm_Ln2, tpm_Ln3, tpm_Ln4, tpm_Ln5, tpm_Ln6, tpm_Ln7, tpm_Ln8, tpm_Ln9, 'rt')


Klibs <- c('tpm_Ct1', 'tpm_Ct2', 'tpm_Ct3', 'tpm_Ct4', 'tpm_Ct5', 'tpm_Ct6', 'tpm_Ct7', 'tpm_Ct8', 'tpm_Ct9', 'tpm_Ct10', 'tpm_Ct11', 'tpm_Ct12', 'tpm_Ct13', 'tpm_Ct14', 'tpm_Ct15', 'tpm_Ct16', 'tpm_Ct17', 'tpm_Ln1', 'tpm_Ln2', 'tpm_Ln3', 'tpm_Ln4', 'tpm_Ln5', 'tpm_Ln6', 'tpm_Ln7', 'tpm_Ln8', 'tpm_Ln9')
Ketiq <- c('KSAT1', 'KSAT2', 'KSAT3', 'KSAT4', 'KSAT5', 'KSAT6', 'KSAT7', 'KSAT8', 'KSAD1', 'KSAD2', 'KSAD3', 'KSAD4', 'KSAD5', 'KSAD6', 'KSAD7', 'KSAD8', 'KSAD9', 'KEL1', 'KEL2', 'KEL3', 'KEL4', 'KEL5', 'KEL6', 'KEL7', 'KEL8', 'KEL9')
#as.double(kallisto_df <- ftp_kallisto[,Klibs])

kallisto_df <- ftp_kallisto[,Klibs]
colnames(kallisto_df) <- Ketiq
kallisto_df <- log2(empty_row(kallisto_df + 1,2))
#Hacer normalizacion por cuantiles
kallisto_df<-normalizeQuantiles(kallisto_df)
plot.Densities(kallisto_df, legend = "topright")
plot.densities(kallisto_df)
#plotDensities(paper_df)



kallisto_mad <- apply(kallisto_df, 1, mad)
kallisto_pca <- prcomp(kallisto_df[order(kallisto_mad, decreasing = T)[0:1000],])

{
#jpeg(file="Kallisto PreCombat 1 v.s. 2.jpg",width = 500, height = 500)
plot(kallisto_pca$rotation[,1], kallisto_pca$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC2", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F,ncol = 2)
#dev.off()

#jpeg(file="Kallisto PreCombat 1 v.s. 3.jpg",width = 500, height = 500)
plot(kallisto_pca$rotation[,1], kallisto_pca$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)
#dev.off()

#jpeg(file="Kallisto PreCombat 2 v.s. 3.jpg",width = 500, height = 500)
plot(kallisto_pca$rotation[,2], kallisto_pca$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC2 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 3)
#dev.off()

kallisto_mat_batch <- c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1)
kallisto_modcombat <- model.matrix(~1, data=data.frame(cmb=kallisto_mat_batch))
kallisto_postComBat <- ComBat(dat= kallisto_df, batch=kallisto_mat_batch, mod=kallisto_modcombat, par.prior=TRUE, prior.plots=FALSE)

kallisto_mad_postcombat <- apply(kallisto_postComBat, 1, mad)
kallisto_pca_postcombat <- prcomp(kallisto_postComBat[order(kallisto_mad_postcombat, decreasing = T)[0:1000],])

#jpeg(file="Kallisto PostCombat 1 v.s. 2.jpg",width = 500, height = 500)
plot(kallisto_pca_postcombat$rotation[,1], kallisto_pca_postcombat$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Kallisto Post ComBat: PC1 vs PC2", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)
#dev.off()

#jpeg(file="Kallisto PostCombat 1 v.s. 3.jpg",width = 500, height = 500)
plot(kallisto_pca_postcombat$rotation[,1], kallisto_pca_postcombat$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Kallisto Post ComBat: PC1 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)
#dev.off()


#jpeg(file="Kallisto PostCombat 2 v.s. 3.jpg",width = 500, height = 500)
plot(kallisto_pca_postcombat$rotation[,2], kallisto_pca_postcombat$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Kallisto Post ComBat: PC2 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 4)
#dev.off()
}

#--------------------------------------------------------------------------------------------------------------------------------------------------
#JUNTAR LOS PLOTS

#Paper(pre,post)Kallisto(pre,post)


#1 vs 2
jpeg(file="1 v.s. 2.jpg",width = 1280, height = 720)
{
par(mfrow=c(2,2))

plot(paper_pca$rotation[,1], paper_pca$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC2", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)

plot(paper_pca_postcombat$rotation[,1], paper_pca_postcombat$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Post ComBat: PC1 vs PC2", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)

plot(kallisto_pca$rotation[,1], kallisto_pca$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC2", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F,ncol = 2)

plot(kallisto_pca_postcombat$rotation[,1], kallisto_pca_postcombat$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Kallisto Post ComBat: PC1 vs PC2", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)
}
dev.off()


#1 vs 3
jpeg(file="1 v.s. 3.jpg",width = 1280, height = 720)
{
par(mfrow=c(2,2))

plot(paper_pca$rotation[,1], paper_pca$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)

plot(paper_pca_postcombat$rotation[,1], paper_pca_postcombat$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Post ComBat: PC1 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)

plot(kallisto_pca$rotation[,2], kallisto_pca$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC2 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 3)

plot(kallisto_pca_postcombat$rotation[,2], kallisto_pca_postcombat$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Kallisto Post ComBat: PC2 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 4)
}
dev.off()


#2 vs 3
jpeg(file="2 v.s. 3.jpg",width = 1280, height = 720)
{
par(mfrow=c(2,2))

plot(paper_pca$rotation[,2], paper_pca$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC2 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 3)

plot(paper_pca_postcombat$rotation[,2], paper_pca_postcombat$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Post ComBat: PC2 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(paper_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 4)

plot(kallisto_pca$rotation[,1], kallisto_pca$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)

plot(kallisto_pca_postcombat$rotation[,1], kallisto_pca_postcombat$rotation[,3], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26), pch=c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Kallisto Post ComBat: PC1 vs PC3", xlim=range(-0.4,0.4), ylim=range(-0.6,0.8))
legend("topleft", legend=rownames(kallisto_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26),pch = c(15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)
}
dev.off()












#--------------------------------------------------------------------------------------------------------------------------------------------------

#####################################################################################################################3
#####################################################################################################################3
#####################################################################################################################3
#----------------------------------------------------------------PCA y comBat de los datos de los autores-----------------------------------------------------{-
{
gz_paper <- file('D:/Users/Disco Duro Intero de Alex/Documents/Bioinformatics Team/2020 LUPUS/Analisis en R/GSE110685_PCarlucci_RNA-Seq_tc_rpkm_combined_submitted.txt', 'rt')
ftp_paper <- read.table(gz_paper,header=T, sep = '\t')
libs <- c('S01', 'S03', 'S14', 'S15', 'S16', 'S17', 'S19', 'S21', 'S02', 'S06', 'S08', 'S10', 'S11', 'S13', 'S20', 'S23', 'S24')
etiq <- c('SAT1', 'SAT2', 'SAT3', 'SAT4', 'SAT5', 'SAT6', 'SAT7', 'SAT8', 'EL1', 'EL2', 'EL3', 'EL4', 'EL5', 'EL6', 'EL7', 'EL8', 'EL9')
paper_df <- ftp_paper[,libs]
colnames(paper_df) <- etiq
paper_df <- log2(empty_row(paper_df + 1,2))
plotDensities(paper_df, legend = "topright")
#plot.densities(paper_df)
plotDensities(paper_df)

paper_mad <- apply(paper_df, 1, mad)
paper_pca <- prcomp(paper_df[order(paper_mad, decreasing = T)[0:200],])
}

plot(paper_pca$rotation[,1], paper_pca$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), pch=c(15,15,15,15,15,15,15,15,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC2")
legend("topleft", legend=rownames(paper_pca$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),pch = c(15,15,15,15,15,15,15,15,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F)

{

paper_mat_batch <- c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2)
paper_modcombat <- model.matrix(~1, data=data.frame(cmb=paper_mat_batch))
paper_postComBat <- ComBat(dat= paper_df, batch=paper_mat_batch, mod=paper_modcombat, par.prior=TRUE, prior.plots=FALSE)

paper_mad_postcombat <- apply(paper_postComBat, 1, mad)
paper_pca_postcombat <- prcomp(paper_postComBat[order(paper_mad_postcombat, decreasing = T)[0:200],])

}

plot(paper_pca_postcombat$rotation[,1], paper_pca_postcombat$rotation[,2], col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), pch=c(15,15,15,15,15,15,15,15,19,19,19,19,19,19,19,19,19), cex = 1, xlab="PC1",ylab="PC2", main = "Post ComBat: PC1 vs PC2")
legend("topright", legend=rownames(paper_pca_postcombat$rotation), col=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),pch = c(15,15,15,15,15,15,15,15,19,19,19,19,19,19,19,19,19), pt.cex = 1, cex=0.75, horiz = F, ncol = 2)

#--------------------------------------------------------------------------------------------------------------------------------------------------





getNormalizedCounts <- function(count_mat){
  samples_classA = 1:round(ncol(count_mat)/2)
  samples_classB = (1:ncol(count_mat))[-samples_classA]
  aux_data <- count_mat
  aux_desc <- data.frame(condition=c(rep("treated",length(samples_classA)),rep("untreated",length(samples_classB))), type=rep("paired-end",c(length(samples_classA)+length(samples_classB))))
  aux_dds <- DESeqDataSetFromMatrix(countData = aux_data, colData = aux_desc, design = ~condition)
  aux_dds <- DESeq(aux_dds)
  normalized_counts <- counts(aux_dds, normalized=T)
  list(dds=aux_dds, norm_counts=normalized_counts)
}

#t.test(exp_data[gen,”controles”], exp_data[gen,”enfermos”])


tpm_adj = aux_tpm_S_E + 1

expressed_tpm_log2 = log2(empty_row(tpm_adj,2))

#Dend <- as.dendrogram(hclust(dist(t(expressed_tpm_log2)),method = "ward.D2"))
#plot(Dend)

Artritis <- colcat(a ='A', x = expressed_tpm_log2, pos = 2)
Crohn <- colcat(a= 'C', x= expressed_tpm_log2, pos=2)
Diabetes <- colcat(a='D', b='J', x=expressed_tpm_log2, pos=2)
Esclerosis <- colcat(a='E', b='J', x=expressed_tpm_log2, pos=2)
Lupus <- colcat(a= 'L', x=expressed_tpm_log2, pos=2)

Norm_Artritis <- normalizeQuantiles(Artritis, ties = TRUE)
Norm_Crohn <- normalizeQuantiles(Crohn, ties = TRUE)
Norm_Diabetes <- normalizeQuantiles(Diabetes, ties = TRUE)
Norm_Esclerosis <- normalizeQuantiles(Esclerosis, ties = TRUE)
Norm_Lupus <- normalizeQuantiles(Lupus, ties = TRUE)

merged_Norm <- empty_row(cbind(Norm_Artritis, Norm_Crohn, Norm_Diabetes, Norm_Esclerosis, Norm_Lupus), 0, equal = TRUE)
merged_preComBat <- normalizeQuantiles(merged_Norm, ties = TRUE)

plotDensities(Artritis, main = 'Artritis pre qn', legend =  'topright')
plotDensities(Norm_Artritis, main = 'Artritis post qn1', legend = 'topright')
plotDensities(merged_preComBat[,1:15], main = 'Artritis post qn2', legend = 'topright')

plotDensities(Crohn, main = 'Crohn pre qn', legend= 'topright')
plot.densities(Norm_Crohn, main = 'Crohn post qn1', legend = 'topright')
plotDensities(merged_preComBat[,16:35], main = 'Crohn post qn2', legend = 'topright')

plotDensities(Diabetes, main = 'Diabetes pre qn', legend = 'topright')
plotDensities(Norm_Diabetes, main = 'Diabetes post qn1', legend = 'topright')
plotDensities(merged_preComBat[,36:43], main = 'Diabetes post qn2', legend = 'topright')

plotDensities(Esclerosis, main = 'Esclerosis pre qn', legend = 'topright')
plotDensities(Norm_Esclerosis, main = 'Esclerosis post qn1', legend = 'topright')
plotDensities(merged_preComBat[,44:51], main = 'Esclerosis post qn2', legend = 'topright')

plotDensities(Lupus, main = 'Lupus pre qn', legend = 'topright')
plotDensities(Norm_Lupus, main = 'Lupus post qn1', legend = 'topright')
plotDensities(merged_preComBat[,52:77], main = 'Lupus post qn2', legend = 'topright')

mad_filt <- apply(merged_preComBat, 1, mad) >= 2

pca_preComBat <- prcomp(merged_preComBat[mad_filt,])

color_vector = c(rep(1, times = length(grep("EA", rownames(pca_preComBat$rotation)))),
                 rep(1, times = length(grep("SA", rownames(pca_preComBat$rotation)))),
                 rep(2, times = length(grep("EC", rownames(pca_preComBat$rotation)))),
                 rep(2, times = length(grep("SC", rownames(pca_preComBat$rotation)))),
                 rep(3, times = length(grep("ED", rownames(pca_preComBat$rotation)))),
                 rep(3, times = length(grep("SJ.", rownames(pca_preComBat$rotation)))),
                 rep(4, times = length(grep("EE", rownames(pca_preComBat$rotation)))),
                 rep(4, times = length(grep("SJ", rownames(pca_preComBat$rotation)))),
                 rep(5, times = length(grep("EL", rownames(pca_preComBat$rotation)))),
                 rep(5, times = length(grep("SL", rownames(pca_preComBat$rotation))))
)

pch_vector =   c(rep(17, times = length(grep("EA", rownames(pca_preComBat$rotation)))),
                 rep(19, times = length(grep("SA", rownames(pca_preComBat$rotation)))),
                 rep(17, times = length(grep("EC", rownames(pca_preComBat$rotation)))),
                 rep(19, times = length(grep("SC", rownames(pca_preComBat$rotation)))),
                 rep(17, times = length(grep("ED", rownames(pca_preComBat$rotation)))),
                 rep(19, times = length(grep("SJ.", rownames(pca_preComBat$rotation)))),
                 rep(17, times = length(grep("EE", rownames(pca_preComBat$rotation)))),
                 rep(19, times = length(grep("SJ", rownames(pca_preComBat$rotation)))),
                 rep(17, times = length(grep("EL", rownames(pca_preComBat$rotation)))),
                 rep(19, times = length(grep("SL", rownames(pca_preComBat$rotation))))
)

#store_lgnd = c()

# #for (i in 1:length(rownames(data_pca_sample$rotation))){
#   if (nchar(rownames(data_pca_sample$rotation)[i]) < 5){
#     store_lgnd[i] = substr(rownames(data_pca_sample$rotation)[i], 1,2)
#   } else {
#     store_lgnd[i] = substr(rownames(data_pca_sample$rotation)[i], 1,4)
#   }
# }

#pca_lgnd <- unique(store_lgnd)
pca_lgnd <- c('EA', 'SA', 'EC', 'SC', 'ED', 'SJ.D', 'EE', 'SJ.E', 'EL', 'SL')

plot(pca_preComBat$rotation[,1], pca_preComBat$rotation[,2], col=color_vector, pch=pch_vector, cex = 1, xlab="PC1",ylab="PC2", main = "Pre ComBat: PC1 vs PC2")
legend("topright", legend=pca_lgnd, col=c(1,1,2,2,3,3,4,4,5,5),pch = c(17,19,17,19,17,19,17,19,17,19), pt.cex = 1, cex=0.75, horiz = F)

plot(pca_preComBat$rotation[,1], pca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC1",ylab="PC3", main = "Pre ComBat: PC1 vs PC3")
legend("topright", legend=pca_lgnd,  col=c(1,1,2,2,3,3,4,4,5,5),pch = c(17,19,17,19,17,19,17,19,17,19), pt.cex = 1, cex=0.75)

plot(pca_preComBat$rotation[,2], pca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC2",ylab="PC3", main = "Pre ComBat: PC2 vs PC3")
legend("topright", legend=pca_lgnd, col=c(1,1,2,2,3,3,4,4,5,5),pch = c(17,19,17,19,17,19,17,19,17,19), pt.cex=1, cex=0.75, horiz = F)

s3d = scatterplot3d(pca_preComBat$rotation[,"PC1"], pca_preComBat$rotation[,"PC2"], pca_preComBat$rotation[,"PC3"], xlab='Comp.1', ylab='Comp.2', zlab='Comp.3', color=color_vector, pch = pch_vector)
s3d.coords = s3d$xyz.convert(pca_preComBat$rotation[,"PC1"], pca_preComBat$rotation[,"PC2"], pca_preComBat$rotation[,"PC3"])
text(s3d.coords$x, s3d.coords$y, labels=rownames(pca_preComBat$rotation), col=color_vector, cex=.8, pos=4)

#Combat
mat_batch <- c(rep(1, ncol(Artritis)), rep(2, ncol(Crohn)), rep(3, ncol(Diabetes)), rep(4, ncol(Esclerosis)), rep(5, ncol(Lupus)))
modcombat <- model.matrix(~1, data=data.frame(cmb=mat_batch))
merged_postComBat <- ComBat(dat= merged_preComBat, batch=mat_batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

mad_filt2 <- apply(merged_postComBat, 1, mad) >= 2.5

mad_postComBat <- merged_postComBat[mad_filt2,]

pca_postComBat <- prcomp(mad_postComBat)

plot(x=pca_postComBat$rotation[,1], y=pca_postComBat$rotation[,2], pch=pch_vector, col=color_vector, xlab="PCA1", ylab="PCA2", main = "Post ComBat: PC1 vs PC2")
legend(x="bottomright", legend=pca_lgnd, col=c(1,1,2,2,3,3,4,4,5,5),pch = c(17,19,17,19,17,19,17,19,17,19))

plot(x=pca_postComBat$rotation[,1], y=pca_postComBat$rotation[,3], pch=pch_vector, col=color_vector, xlab="PCA1", ylab="PCA3", main = "Post ComBat: PC1 vs PC3")
legend(x="bottomright", legend=pca_lgnd, col=c(1,1,2,2,3,3,4,4,5,5),pch = c(17,19,17,19,17,19,17,19,17,19))

plot(x=pca_postComBat$rotation[,2], y=pca_postComBat$rotation[,3], pch=pch_vector, col=color_vector, xlab="PCA2", ylab="PCA3", main = "Post ComBat: PC2 vs PC3")
legend(x="bottomright", legend=pca_lgnd, col=c(1,1,2,2,3,3,4,4,5,5),pch = c(17,19,17,19,17,19,17,19,17,19))

names <- row.names.data.frame(resMat[[1]])
gene_pos <- match(names, genes$transcript_id)
gene_name <- genes$gene_name
gene_type <- genes$gene_type
transcript_type <- genes$transcript_type

resMat[[1]] <- cbind(gene_name[gene_pos], resMat[[1]], gene_type[gene_pos], transcript_type[gene_pos])
estadisticos <- cbind(resMat[[1]]$`gene_name[gene_pos]`, resMat[[1]]$padj, resMat[[1]]$log2FoldChange)

expresion <- estadisticos
rownames(expresion) <- rownames(resMat[[1]])

nombres_expresion <- c("gene_name[gene_pos]", "padj", "log2FoldChange")

for (i in 1:length(archivos_tsv)) {#Genera variables que contienen las columnas de tpm 
  if (nchar(archivos_tsv[i])== 7) {
    nombres_expresion <- c(nombres_expresion, paste0("tpm_",substr(archivos_tsv[i],1,3)))
  } else {
    nombres_expresion <-  c(nombres_expresion, paste0("tpm_",substr(archivos_tsv[i],1,4)))
  }
}

for (i in 1:length(abundancias)) {
  expresion <- cbind(expresion, abundancias[[i]]$tpm)
}
expresion <- as.data.frame(expresion)


colnames(expresion) <- nombres_expresion
expresion[, 2:ncol(expresion)] <- lapply(expresion[, 4:ncol(expresion)], function(x) as.numeric(as.character(x)))

seleccion <- expresion %>%
  rownames_to_column("nombres") %>%
  filter(padj < .01) %>%
  filter(abs(log2FoldChange)>= 30) %>%
  column_to_rownames("nombres")
sel_idx <- (rowSums(seleccion[,4:ncol(seleccion)] > 1) > 1)
expS_E <- seleccion[sel_idx,]


scale_aux_tpm_S_E <- t(apply(aux_tpm_S_E[rownames(expS_E),],1,function(x) scale(x)))
colnames(scale_aux_tpm_S_E) <- colnames(aux_tpm_S_E)
rownames(scale_aux_tpm_S_E) <-gene_name[match(row.names(scale_aux_tpm_S_E), genes$transcript_id)] 
quant<- quantile(scale_aux_tpm_S_E)

#colors_S_E <- c(seq(quant[1],quant[2],length=100),seq(quant[2]+.001,quant[4],length=100),seq(quant[4]+.001,quant[5],length=100))

#my_palette <- colorRampPalette(c("green", "white", "red"))(n = 299)

#HM_S_E <- heatmap.2(as.matrix(scale_aux_tpm_S_E), col= my_palette, 
#          breaks= colors_S_E, density.info="none", trace="none", 
#         dendrogram=c("column"), symm=F,symkey=F,symbreaks=T, scale="none")

HM_S_E <- heatmap.2(as.matrix(scale_aux_tpm_S_E), density.info="none", trace="none", 
                    dendrogram=c("column"), symm=F,symkey=F,symbreaks=T, scale="none")
# 
# HC_MS_name <- as.character(resMat[[1]]$`gene_name[gene_pos]`)
# HC_MS_pvalue_logfold <- -log10(resMat[[1]]$pvalue)*sign(resMat[[1]]$log2FoldChange)
# NF_HC_MS <- data.frame(HC_MS_name, HC_MS_pvalue_logfold, stringsAsFactors = FALSE)#Non filtered
# HC_MS_filter <- match(unique(HC_MS_name), HC_MS_name)
# Filtered_HC_MS <- NF_HC_MS[HC_MS_filter,]
# HC_MS_order <- order(Filtered_HC_MS$HC_MS_pvalue_logfold, decreasing = TRUE)
# GSEA_HC_MS <- Filtered_HC_MS[HC_MS_order,]
# colnames(GSEA_HC_MS) <- c("geneName", "Log10(pvalue)*-sign(log2foldchange)")
# 
# GSEA_list <- list("GSEA_HC_MS" = GSEA_HC_MS, "GSEA_HC_DM" = GSEA_HC_DM, "GSEA_MS_DM" = GSEA_MS_DM)



