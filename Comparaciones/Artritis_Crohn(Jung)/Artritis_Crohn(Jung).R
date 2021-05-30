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
library(edgeR)
library(PanNETassigner)


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

empty_row <- function(mat, ref, equal = FALSE){
  if (equal==FALSE) {
    filt = apply(mat, 1, function(x) ifelse(all(x < ref), 0, 1))
    filtered = mat[which(filt == 1),]
  } else if (equal == TRUE) {
    filt = apply(mat, 1, function(x) ifelse(all(x == ref), 0, 1))
    filtered = mat[which(filt == 1),]
  }
  return(filtered)
}
#-----------------------------------Preparando datos de Crohn----------------------------------------------------------------------

setwd("C:/Users/tocht/Downloads")

load("C:/Users/tocht/Downloads/rse_gene_blood.Rdata")
row <- rowData(rse_gene) #used to view the range information
col <- colData(rse_gene) #Sample meta-data describing the samples as dataframe

rpkm_control <- rpkm.SummarizedExperiment(rse_gene, gene.length = row$bp_length, normalized.lib.sizes = TRUE)
control_genematch <- match(rownames(rpkm_control), row$gene_id)
control_names <- unlist(row$symbol)
rpkm_control <- rpkm_control[!is.na(control_names[control_genematch]),]
rownames(rpkm_control) <- control_names[control_genematch][!is.na(control_names[control_genematch])]
rpkm_control <- rpkm_control[!duplicated(rownames(rpkm_control)),]

celltype_match <- match(colnames(rpkm_control), col$run)
celltype_filter <- col$smtsd[celltype_match]
celltype_filter <- celltype_filter[!is.na(celltype_filter)] == "Whole Blood"

rpkm_control <- rpkm_control[, celltype_filter]

control_labs <- c()
for (i in 1:length(colnames(rpkm_control))) {
  control_labs <- c(control_labs, paste("SG", i, sep = ""))
}
colnames(rpkm_control) <- control_labs
rpkm_control <- as.data.frame(empty_row(rpkm_control, ref = 0, equal = T))

plotDensities(log2(rpkm_control+1), legend = "right", main = "Controles Gtex no normalizados")
rpkm_control <- log2(normalizeQuantiles(rpkm_control)+1)
plotDensities(rpkm_control, legend = "right", main = "Controles Gtex normalizados")

data_Crohn <- read.csv("totalFile.csv", header = TRUE, stringsAsFactors = FALSE )
aux_Crohn <- data_Crohn[,-1:-3]
prom_C <- apply(aux_Crohn, 1, mean)

data_Crohn <- data_Crohn[order(prom_C, decreasing = T),]
data_Crohn <- data_Crohn[!duplicated(data_Crohn$Gene_name),]
rownames(data_Crohn) <- data_Crohn$Gene_name
data_Crohn <- as.data.frame(empty_row(data_Crohn, ref = 0 , equal = T))

plotDensities(log2(data_Crohn[,-1:-3] + 1), legend = "right", main = "Enfermos Crohn no normalizados")
data_Crohn <- log2(normalizeQuantiles(data_Crohn[,-1:-3]) + 1)
plotDensities(data_Crohn, legend = "right", main = "Enfermos Crohn no normalizados")

match_crohn_control <- match(rownames(rpkm_control), rownames(data_Crohn))
match_crohn_control <- match_crohn_control[!is.na(match_crohn_control)]
crohn_sbst <- data_Crohn[match_crohn_control,]
control_sbst <- rpkm_control[match(rownames(crohn_sbst), rownames(rpkm_control)),]

# crohn_mad_nmf <- apply(crohn_sbst, 1, mad)
# crohn_nmf <- crohn_sbst[order(crohn_mad_nmf, decreasing = T),][1:1500,]
# setwd("C://Users/tocht/Documents/Graficas NMF/Crohn_1500/")
# 
# nmfconsensus(
#   input.ds =           "",
#   k.init =             2,
#   k.final =            5,
#   num.clusterings =    100,
#   maxniter =           1000,
#   error.function =     "euclidean",
#   rseed =              1234,
#   stopconv =           40,
#   stopfreq =           10,
#   non.interactive.run = T,
#   doc.string =         "crohn_mat",
#   data = crohn_nmf)

crohn_clusters <- read.table("C://Users/tocht/Documents/Graficas NMF/Crohn_1500/crohn_matconsensus.k.2.gct", header = T, sep = "\t")

crohn_clust1_samples <- crohn_clusters[crohn_clusters[,2] == 1,1]
crohn_clust2_samples <- crohn_clusters[crohn_clusters[,2] == 2,1]

mad_Crohn <- apply(data_Crohn, 1, mad)
pca_Crohn <- prcomp(data_Crohn[order(mad_Crohn, decreasing = T)[0:2000], c(crohn_clust1_samples,crohn_clust2_samples)])

crohn_col_vector <- c(rep(1, times = length(crohn_clust1_samples)),
                  rep(2, times = length(crohn_clust2_samples)))

plot(pca_Crohn$rotation[,1], pca_Crohn$rotation[,2], col = crohn_col_vector, pch = 15 ,xlab= paste("PCA1: ", round(summary(pca_Crohn)$importance[2,1]*100,1),"%", sep=""), ylab=paste("PCA2: ", round(summary(pca_Crohn)$importance[2,2]*100,1),"%",sep=""),  cex = 1.1,  main = "Autores Crohn: PC1 vs PC2")
legend("topright", legend=c("Cluster 1", "Cluster 2"), col=c(1,2), pch = c(15,15), pt.cex = 1, cex=0.6, horiz = F, ncol = 1, inset= c(0.05,0.1))
#legend("top", legend=c("Control", "Enfermedad"), pch = c(2,0), pt.cex = 1, cex=0.6, horiz = F, ncol = 1, inset=c(0.1,0.35))
#text(pca_Crohn$rotation[,1], pca_Crohn$rotation[,2], labels = rownames(pca_Crohn$rotation))

# control_mad_nmf <- apply(control_sbst, 1, mad)
# control_nmf <- control_sbst[order(control_mad_nmf, decreasing = T),][1:1500,]
# setwd("C://Users/tocht/Documents/Graficas NMF/Control_1500_2/")


# timestamp()
# nmfconsensus(
#   input.ds =           "",
#   k.init =             2,
#   k.final =            5,
#   num.clusterings =    100,
#   maxniter =           1000,
#   error.function =     "euclidean",
#   rseed =              1234,
#   stopconv =           40,
#   stopfreq =           10,
#   non.interactive.run = T,
#   doc.string =         "control_mat",
#   data = control_nmf)
# timestamp()

control_clusters <- read.table("C://Users/tocht/Documents/Graficas NMF/Control_1500/control_matconsensus.k.2.gct", header = T, sep = "\t")


control_clusters1_samples <- control_clusters[control_clusters[,2] == 1,1]
control_clusters2_samples <- control_clusters[control_clusters[,2] == 2,1]

data_crohn_control <- log2(normalizeQuantiles(data.frame(crohn_sbst[,crohn_clust1_samples], control_sbst[,control_clusters1_samples])) + 1)
prep_mad_filt <- apply(data_crohn_control, 1, mad)
crohnprep_pca <- prcomp(data_crohn_control[order(prep_mad_filt, decreasing = T)[1:2000],])


color_vector <- c(rep(1, times = length(grep("EC", rownames(crohnprep_pca$rotation)))),
                  rep(2, times = length(grep("SG", rownames(crohnprep_pca$rotation))))
)
pch_vector <- c(rep(15, times = length(grep("EC", rownames(crohnprep_pca$rotation)))),
                rep(17, times = length(grep("EA", rownames(crohnprep_pca$rotation))))
)

plot(crohnprep_pca$rotation[,1], crohnprep_pca$rotation[,2], pch = pch_vector, col= color_vector, xlab= paste("PCA1: ", round(summary(crohnprep_pca)$importance[2,1]*100,1),"%", sep=""), ylab=paste("PCA2: ", round(summary(crohnprep_pca)$importance[2,2]*100,1),"%",sep=""),  cex = 1.1,  main = "Autores Crohn preparacion con controles GTEX: PC1 vs PC2")
#legend("topleft", legend=rownames(pca_Crohn$rotation), pt.cex = 1, cex=0.6, horiz = F, ncol = 5, inset=c(0.1,0.25))
#text(crohnprep_pca$rotation[,1], crohnprep_pca$rotation[,2], labels = rownames(crohnprep_pca$rotation))

#-----------------------------------Carga de la anotacion del genoma---------------------------------------------------------------
genes <- as.data.frame(rtracklayer::import('D:/Archivos_de_trabajo/Analisis_combinado/Homo_sapiens.GRCh38.96.chr_patch_hapl_scaff.gtf.gz'))


#-----------------------------------Preparando controles de Artirtis para juntar con enfermos de Crohn-----------------------------
paper_RA <- gzfile('D:/Archivos_de_trabajo/Analisis_combinado/GSE90081_5nRA_12HC_7tRA_fpkm_table.txt.gz', 'rt')
data_Artritis <- read.table(paper_RA,header=T, sep = '\t')

#Orden por promedios
prom_A <- apply(data_Artritis[-1:-2], 1, mean)
data_Artritis <- data_Artritis[order(prom_A, decreasing = T),]

#Eliminacion de Na y repetidos
data_Artritis <- data_Artritis[!duplicated(data_Artritis$Gene_ID),]
data_Artritis <- data_Artritis[!is.na(data_Artritis$Gene_ID),]

#Matriz final de controles de artritis con log2
rownames(data_Artritis) <- data_Artritis$Gene_ID
data_Artritis <- normalizeQuantiles(log2(data_Artritis[,-1:-2] + 1))
samples_A <- c('ntRA_01', 'ntRA_02', 'ntRA_03', 'ntRA_04', 'ntRA_05', 'HC_01', 'HC_02', 'HC_03', 'HC_04', 'HC_05', 'HC_06', 'HC_07', 'HC_08', 'HC_09', 'HC_10' )
etiq_A <- c('EA1', 'EA2', 'EA3', 'EA4', 'EA5', 'SA1', 'SA2', 'SA3', 'SA4', 'SA5', 'SA6', 'SA7', 'SA8', 'SA9', 'SA10')

data_Artritis <- data_Artritis[,samples_A]
colnames(data_Artritis) <- etiq_A
data_Artritis <- as.data.frame(empty_row(data_Artritis, ref = 0, equal = T))

#Match de genes para juntarse en una sola matriz
gene_match_A <- match(rownames(data_Artritis), genes$gene_name)
gene_match_A <- gene_match_A[!is.na(gene_match_A)]
gene_match_C <- match(rownames(data_Artritis), rownames(data_crohn_control))
gene_match_C <- gene_match_C[!is.na(gene_match_C)]
Crohn_Artritis<- normalizeQuantiles(data.frame(data_crohn_control[gene_match_C,], data_Artritis[match(rownames(data_crohn_control[gene_match_C,]), rownames(data_Artritis)),]))

#Filtro con MAD y PCA
mad_Crohn_artritis <- apply(Crohn_Artritis, 1, mad)
pca_Crohn_artritis <- prcomp(Crohn_Artritis[order(mad_Crohn_artritis, decreasing = T)[0:2000],])

#Grafica de PCA
color_vector1 <- c(rep(1, times = length(grep("EC", rownames(pca_Crohn_artritis$rotation)))),
                  rep(2, times = length(grep("SG", rownames(pca_Crohn_artritis$rotation)))),
                  rep(3, times = length(grep("EA", rownames(pca_Crohn_artritis$rotation)))),
                  rep(4, times = length(grep("SA", rownames(pca_Crohn_artritis$rotation))))
)
pch_vector1 =   c(rep(15, times = length(grep("EC", rownames(pca_Crohn_artritis$rotation)))),
                 rep(17, times = length(grep("SG", rownames(pca_Crohn_artritis$rotation)))),
                 rep(15, times = length(grep("EA", rownames(pca_Crohn_artritis$rotation)))),
                 rep(17, times = length(grep("SA", rownames(pca_Crohn_artritis$rotation))))
)

plot(pca_Crohn_artritis$rotation[,1], pca_Crohn_artritis$rotation[,2], pch = pch_vector1, col = color_vector1, xlab= paste("PCA1: ", round(summary(pca_Crohn_artritis)$importance[2,1]*100,1),"%", sep=""), ylab=paste("PCA2: ", round(summary(pca_Crohn_artritis)$importance[2,2]*100,1),"%",sep=""),  cex = 1.1,  main = "Autores Crohn vs Artritis: PC1 vs PC2")



#Correccion del batch effect con ComBat
batches <- c(rep(1, times = length(grep("EC", rownames(pca_Crohn_artritis$rotation)))),
                              rep(1, times = length(grep("SG", rownames(pca_Crohn_artritis$rotation)))),
                              rep(2, times = length(grep("EA", rownames(pca_Crohn_artritis$rotation)))),
                              rep(2, times = length(grep("SA", rownames(pca_Crohn_artritis$rotation))))
)
modcombat <- model.matrix(~1, data = data.frame(cmb=batches))
Crohn_combat <- ComBat(dat= Crohn_Artritis, batch = batches, mod = modcombat, par.prior = T, prior.plots = F)
Crohn_combat <- Crohn_combat - range(Crohn_combat)[1]

#Filtro con MAD y PCA despues de ComBat
mad_Crohn_PostCombat <- apply(Crohn_combat, 1, mad)
pca_Crohn_PostCombat <- prcomp(Crohn_combat[order(mad_Crohn_PostCombat, decreasing = T)[0:2000],])

#Grafica de PCA de datos sin batch effect
plot(pca_Crohn_PostCombat$rotation[,1], pca_Crohn_PostCombat$rotation[,2], pch = pch_vector1, col = color_vector1, xlab= paste("PCA1: ", round(summary(pca_Crohn_PostCombat)$importance[2,1]*100,1),"%", sep=""), ylab=paste("PCA2: ", round(summary(pca_Crohn_PostCombat)$importance[2,2]*100,1),"%",sep=""),  cex = 1.1,  main = "Autores Crohn con controles (PostComBat): PC1 vs PC2")
