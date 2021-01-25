library (DESeq2)
library(refGenome)
library(tximport)
library(stringr)
library(gplots)
library(ActivePathways)
library(rtracklayer)
library(dplyr)
library(tibble)
library(scatterplot3d)
library(RColorBrewer)
library(BiocManager)
library(limma)
library(dplyr)
setwd("C:/Users/angel/Documents/Bioinformatica/AR/tsvs")
AR1<-read.table(file="RA1ab.tsv", sep = "\t", header=TRUE)
AR2<-read.table(file="RA2ab.tsv", sep = "\t", header=TRUE)
AR3<-read.table(file="RA3ab.tsv", sep = "\t", header=TRUE)
AR4<-read.table(file="RA4ab.tsv", sep = "\t", header=TRUE)
AR5<-read.table(file="RA5ab.tsv", sep = "\t", header=TRUE)
PC1<-read.table(file="SA1.tsv", sep = "\t", header=TRUE)
PC2<-read.table(file="SA2.tsv", sep = "\t", header=TRUE)
PC3<-read.table(file="SA3.tsv", sep = "\t", header=TRUE)
PC4<-read.table(file="SA4.tsv", sep = "\t", header=TRUE)
PC5<-read.table(file="SA5.tsv", sep = "\t", header=TRUE)
PC6<-read.table(file="SA6.tsv", sep = "\t", header=TRUE)
PC7<-read.table(file="SA7.tsv", sep = "\t", header=TRUE)
PC8<-read.table(file="SA8.tsv", sep = "\t", header=TRUE)
PC9<-read.table(file="SA9.tsv", sep = "\t", header=TRUE)
PC10<-read.table(file="SA10.tsv", sep = "\t", header=TRUE)




automatedDESeq2=function(count_mat=count_mat, comparisons=comparisons){
  # Si no hacemos esta línea, y hay genes con ceros, sale pval/qval como NA
  # count_mat = count_mat +1
  
  #Renglones -> Genes, Columnas -> Replicas
  resList=list()
  for (n in 1:length(comparisons)){
    controls = comparisons[[n]][1]
    samples = comparisons[[n]][2]
    samples_classA = colnames(count_mat)[grep(controls, colnames(count_mat))]
    samples_classB = colnames(count_mat)[grep(samples, colnames(count_mat))]
    aux_data <- count_mat[,c(samples_classA,samples_classB)]
    aux_desc <- data.frame(condition=c(rep("treated",length(samples_classA)),rep("untreated",length(samples_classB))), type=rep("paired-end",c(length(samples_classA)+length(samples_classB))))
    aux_dds <- DESeqDataSetFromMatrix(countData = aux_data, colData = aux_desc, design = ~condition)
    aux_dds <- DESeq(aux_dds)
    aux_res <- as.data.frame(results(aux_dds))
    # ponemos validaciones para que no salgan NAs
    aux_res$log2FoldChange = sapply(aux_res$log2FoldChange, function(x) ifelse(is.na(x), 0, x))
    aux_res$pvalue = sapply(aux_res$pvalue, function(x) ifelse(is.na(x), 1, x))
    aux_res$padj = sapply(aux_res$padj, function(x) ifelse(is.na(x), 1, x))
    
    resList[[paste(controls,"_",samples,sep="")]] = cbind(aux_res, counts(aux_dds, normalized=T))
    # colnames(resList[[paste(controls,"_",samples,sep="")]]) = c("baseMean","log2FoldChange","lfcSE","stat", "pvalue", "padj", paste(controls[n],"_1",sep=""), paste(controls[n],"_2",sep=""), paste(samples[n],"_1",sep=""), paste(samples[n],"_2",sep=""))	
    colnames(resList[[paste(controls,"_",samples,sep="")]]) = c("baseMean","log2FoldChange","lfcSE","stat", "pvalue", "padj", samples_classA, samples_classB)	
  }
  return(resList)
}

count_mat <- cbind(AR1$est_counts, AR2$est_counts, AR3$est_counts, AR4$est_counts, AR5$est_counts,  PC1$est_counts, PC2$est_counts, PC3$est_counts, PC4$est_counts, PC5$est_counts, PC6$est_counts, PC7$est_counts, PC8$est_counts, PC9$est_counts,PC10$est_counts )
rownames(count_mat) <- AR1$target_id
colnames(count_mat) <- c("AR_1", "AR_2","AR_3", "AR_4","AR_5","PC_1","PC_2","PC_3","PC_4","PC_5","PC_6","PC_7","PC_8","PC_9","PC_10")
count_mat1 <- count_mat + 1
count_mat1 <- round(count_mat1)
seleccion <- count_mat1<=0
seleccion
count_mat1[seleccion]
a<- count_mat1[seleccion]
comparisons <- list(c("AR","PC"))
comparisons
resMat <- automatedDESeq2(count_mat1, comparisons)
resMat

resFrame<-resMat[[1]]

setwd("C:/Users/angel/Documents/Bioinformatica/AR")

##PCA y plot densities



tpm_aux<-cbind(AR1$tpm, AR2$tpm, AR3$tpm, AR4$tpm, AR5$tpm,  PC1$tpm, PC2$tpm, PC3$tpm, PC4$tpm, PC5$tpm, PC6$tpm, PC7$tpm, PC8$tpm, PC9$tpm,PC10$tpm)
colnames(tpm_aux) <- c("AR_1", "AR_2","AR_3", "AR_4","AR_5","PC_1","PC_2","PC_3","PC_4","PC_5","PC_6","PC_7","PC_8","PC_9","PC_10")
rownames(tpm_aux)<-AR1$target_id
tpm_adj<-tpm_aux+1
expressed_tpm_log2 = log2(empty_file(tpm_adj,2))
tpmnorm<-normalizeQuantiles(expressed_tpm_log2)
plot.densities(expressed_tpm_log2, main = 'Artritis pre qn (tpm)', legend.cex = .8, legend.cols = 4, legend.pos = 'bottomleft')
plot.densities(tpmnorm, main = 'Artritis post qn (tpm)', legend.cex = .8, legend.cols = 4, legend.pos = 'bottomleft')



fpkm<-read.table(gzfile("FPKM.txt.gz"),sep="\t",header=TRUE)
fpkmdf<-as.data.frame(fpkm[3:17])
colnames(fpkmdf) <- c("AR_1", "AR_2","AR_3", "AR_4","AR_5","PC_1","PC_2","PC_3","PC_4","PC_5","PC_6","PC_7","PC_8","PC_9","PC_10")
fpkm_adj<-fpkmdf+1
expressed_fpkm_log2 = log2(empty_file(fpkm_adj,2))
fpkmnorm<- normalizeQuantiles(expressed_fpkm_log2)
plot.densities(expressed_fpkm_log2, main = 'Artritis pre qn (fpkm)', legend.cex = .8, legend.cols = 4, legend.pos = 'bottomleft')
plot.densities(fpkmnorm, main = 'Artritis post qn (fpkm)', legend.cex = .8, legend.cols = 4, legend.pos = 'bottomleft')


#PCA with no name matching or ordering
# tmad_filt <- apply(tpmnorm, 1, mad) 
# 
# tpca_preComBat <- prcomp(tpmnorm[order(tmad_filt, decreasing = T)[0:1000],])
# 
# color_vector = c(rep(1, times = length(grep("AR", rownames(tpca_preComBat$rotation)))),
#                  rep(1, times = length(grep("PC", rownames(tpca_preComBat$rotation)))))
# pch_vector =   c(rep(17, times = length(grep("AR", rownames(tpca_preComBat$rotation)))),
#                  rep(19, times = length(grep("PC", rownames(tpca_preComBat$rotation)))))
# pca_lgnd <- c('EA', 'SA')
# 
# plot(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,2], col=color_vector, pch=pch_vector, cex = 1, xlab="PC1",ylab="PC2", main = "tpm Pre ComBat: PC1 vs PC2")
# legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75, horiz = F)
# text(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,2], labels = row.names(tpca_preComBat$rotation), pos = 4,cex=0.5)
# 
# 
# plot(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC1",ylab="PC3", main = "tpm Pre ComBat: PC1 vs PC3")
# legend("topright", legend=pca_lgnd,  col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75)
# text(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,3], labels = row.names(tpca_preComBat$rotation), pos = 4,cex=0.5)
# 
# 
# plot(tpca_preComBat$rotation[,2], tpca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC2",ylab="PC3", main = "tpm Pre ComBat: PC2 vs PC3")
# legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex=1, cex=0.75, horiz = F)
# text(tpca_preComBat$rotation[,2], tpca_preComBat$rotation[,3], labels = row.names(tpca_preComBat$rotation), pos = 4,cex=0.5)
# 
# 
# fmad_filt <- apply(fpkmnorm, 1, mad) 
# 
# pca_preComBat <- prcomp(fpkmnorm[order(fmad_filt, decreasing = T)[0:1000],])
# 
# color_vector = c(rep(1, times = length(grep("AR", rownames(pca_preComBat$rotation)))),
#                  rep(1, times = length(grep("PC", rownames(pca_preComBat$rotation)))))
# pch_vector =   c(rep(17, times = length(grep("AR", rownames(pca_preComBat$rotation)))),
#                  rep(19, times = length(grep("PC", rownames(pca_preComBat$rotation)))))
# pca_lgnd <- c('EA', 'SA')
# 
# plot(pca_preComBat$rotation[,1], pca_preComBat$rotation[,2], col=color_vector, pch=pch_vector, cex = 1, xlab="PC1",ylab="PC2", main = "fpkm Pre ComBat: PC1 vs PC2")
# legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75, horiz = F)
# text(pca_preComBat$rotation[,1], pca_preComBat$rotation[,2], labels = row.names(pca_preComBat$rotation), pos = 4,cex=0.5)
# 
# 
# plot(pca_preComBat$rotation[,1], pca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC1",ylab="PC3", main = "fpkm Pre ComBat: PC1 vs PC3")
# legend("topright", legend=pca_lgnd,  col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75)
# text(pca_preComBat$rotation[,1], pca_preComBat$rotation[,3], labels = row.names(pca_preComBat$rotation), pos = 4,cex=0.5)

# 
# plot(pca_preComBat$rotation[,2], pca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC2",ylab="PC3", main = "fpkm Pre ComBat: PC2 vs PC3")
# legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex=1, cex=0.75, horiz = F)
# text(pca_preComBat$rotation[,2], pca_preComBat$rotation[,3], labels = row.names(pca_preComBat$rotation), pos = 4,cex=0.5)
# 


#finding genenames of tpmnorm
human_index<- rtracklayer::import('gencode.v30.chr_patch_hapl_scaff.annotation.gtf.gz')
human_df=as.data.frame(human_index)
gene_pos<-match(rownames(tpmnorm),human_df$transcript_id)
gene_name<-human_df$gene_name[gene_pos]
tpmnorm_name<-cbind(tpmnorm)
tpmnorm_name<-as.data.frame(tpmnorm_name)
tpmnorm_name$gene_name=gene_name
fpk_pos<-match(rownames(fpkmnorm),rownames(fpkm))
fpkm_name<-fpkm$Gene_ID[fpk_pos]
fpkmnorm_name<-cbind(fpkm_name,fpkmnorm)

#eliminar genes repetidos con el valor menor de expresion
tpm_avg<-rowMeans(tpmnorm_name[,1:15])
order_tpm<-tpmnorm_name[order(tpm_avg,decreasing=TRUE),]
order_tpm<-distinct(order_tpm,order_tpm$gene_name,.keep_all=TRUE)
#46625 genes en database original solo 12663 genes no repetidos
#ordenar por dispersion
tmad_filt <- apply(order_tpm[,1:15], 1, mad) 
fin_tpm<-order_tpm[order(tmad_filt, decreasing = T),]
matching_tpm<-fin_tpm

#hacer match para graficar fpkm dependiente de tpm
fpkm_avg<-rowMeans(fpkmnorm_name[,2:16])
order_fpkm<-fpkmnorm_name[order(fpkm_avg,decreasing=TRUE),]
order_fpkm<-distinct(order_fpkm,order_fpkm$fpkm_name,.keep_all=TRUE)

fpkm_pos<-match(matching_tpm$gene_name,order_fpkm$fpkm_name)
#166 genes identificados en el top 1000 de fpkm no existen en los gene names de tpm 
matched_fpkm<-order_fpkm[fpkm_pos,2:16]

final_fpkm<-matched_fpkm[rowSums(is.na(matched_fpkm)) != ncol(matched_fpkm),1:15]

final_tpm<-matching_tpm[rowSums(is.na(matched_fpkm)) != ncol(matched_fpkm),1:15]

final_fpkm<-final_fpkm[1:2000,]
final_tpm<-final_tpm[1:2000,]


tpca_preComBat <- prcomp(final_tpm)

color_vector = c(rep(1, times = length(grep("AR", rownames(tpca_preComBat$rotation)))),
                 rep(1, times = length(grep("PC", rownames(tpca_preComBat$rotation)))))
pch_vector =   c(rep(17, times = length(grep("AR", rownames(tpca_preComBat$rotation)))),
                 rep(19, times = length(grep("PC", rownames(tpca_preComBat$rotation)))))
pca_lgnd <- c('EA', 'SA')

plot(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,2], col=color_vector, pch=pch_vector, cex = 1, xlab=paste("PCA2: ", round(summary(tpca_preComBat)$importance[2,2]*100,1),"%", sep=""),ylab=paste("PCA1: ", round(summary(tpca_preComBat)$importance[2,1]*100,1),"%",sep=""), main = "tpm Selected Genes: PC1 vs PC2")
legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75, horiz = F)
text(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,2], labels = row.names(tpca_preComBat$rotation), pos = 4,cex=0.5)


# plot(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC1",ylab="PC3", main = "tpm selected genes: PC1 vs PC3")
# legend("topright", legend=pca_lgnd,  col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75)
# text(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,3], labels = row.names(tpca_preComBat$rotation), pos = 4,cex=0.5)
# 
# 
# plot(tpca_preComBat$rotation[,2], tpca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC2",ylab="PC3", main = "tpm selected genes: PC2 vs PC3")
# legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex=1, cex=0.75, horiz = F)
# text(tpca_preComBat$rotation[,2], tpca_preComBat$rotation[,3], labels = row.names(tpca_preComBat$rotation), pos = 4,cex=0.5)


pca_preComBat <- prcomp(final_fpkm)

color_vector = c(rep(1, times = length(grep("AR", rownames(pca_preComBat$rotation)))),
                 rep(1, times = length(grep("PC", rownames(pca_preComBat$rotation)))))
pch_vector =   c(rep(17, times = length(grep("AR", rownames(pca_preComBat$rotation)))),
                 rep(19, times = length(grep("PC", rownames(pca_preComBat$rotation)))))
pca_lgnd <- c('EA', 'SA')

plot(pca_preComBat$rotation[,1], pca_preComBat$rotation[,2], col=color_vector, pch=pch_vector, cex = 1, xlab=paste("PCA2: ", round(summary(pca_preComBat)$importance[2,2]*100,1),"%", sep=""),ylab=paste("PCA1: ", round(summary(pca_preComBat)$importance[2,1]*100,1),"%",sep=""), main = "fpkm Selected Genes: PC1 vs PC2")
legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75, horiz = F)
text(pca_preComBat$rotation[,1], pca_preComBat$rotation[,2], labels = row.names(pca_preComBat$rotation), pos = 4,cex=0.5)

# 
# plot(pca_preComBat$rotation[,1], pca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC1",ylab="PC3", main = "fpkm selected genes: PC1 vs PC3")
# legend("topright", legend=pca_lgnd,  col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75)
# text(pca_preComBat$rotation[,1], pca_preComBat$rotation[,3], labels = row.names(pca_preComBat$rotation), pos = 4,cex=0.5)
# 
# 
# plot(pca_preComBat$rotation[,2], pca_preComBat$rotation[,3], col=color_vector, pch=pch_vector, cex = 1, xlab="PC2",ylab="PC3", main = "fpkm selected genes: PC2 vs PC3")
# legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex=1, cex=0.75, horiz = F)
# text(pca_preComBat$rotation[,2], pca_preComBat$rotation[,3], labels = row.names(pca_preComBat$rotation), pos = 4,cex=0.5)

# graficar independientmente fpkm y tpm

cor(final_fpkm,final_tpm,method="spearman")
cor(final_fpkm,final_tpm,method="pearson")


ind_tpm<-matching_tpm[1:2000,1:15]
fmad_filt <- apply(order_fpkm[,2:16], 1, mad) 
fin_fpkm<-order_fpkm[order(fmad_filt, decreasing = T),]
ind_fpkm<-fin_fpkm[1:2000,2:16]

tpca_preComBat <- prcomp(ind_tpm)

color_vector = c(rep(1, times = length(grep("AR", rownames(tpca_preComBat$rotation)))),
                 rep(1, times = length(grep("PC", rownames(tpca_preComBat$rotation)))))
pch_vector =   c(rep(17, times = length(grep("AR", rownames(tpca_preComBat$rotation)))),
                 rep(19, times = length(grep("PC", rownames(tpca_preComBat$rotation)))))
pca_lgnd <- c('EA', 'SA')

plot(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,2], col=color_vector, pch=pch_vector, cex = 1, xlab=paste("PCA2: ", round(summary(tpca_preComBat)$importance[2,2]*100,1),"%", sep=""),ylab=paste("PCA1: ", round(summary(tpca_preComBat)$importance[2,1]*100,1),"%",sep=""), main = "Kallisto: Pre Combat PC1 vs PC2")
legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75, horiz = F)
text(tpca_preComBat$rotation[,1], tpca_preComBat$rotation[,2], labels = row.names(tpca_preComBat$rotation), pos = 4,cex=0.5)


pca_preComBat <- prcomp(ind_fpkm)

color_vector = c(rep(1, times = length(grep("AR", rownames(pca_preComBat$rotation)))),
                 rep(1, times = length(grep("PC", rownames(pca_preComBat$rotation)))))
pch_vector =   c(rep(17, times = length(grep("AR", rownames(pca_preComBat$rotation)))),
                 rep(19, times = length(grep("PC", rownames(pca_preComBat$rotation)))))
pca_lgnd <- c('EA', 'SA')

plot(pca_preComBat$rotation[,1], pca_preComBat$rotation[,2], col=color_vector, pch=pch_vector, cex = 1, xlab=paste("PCA2: ", round(summary(pca_preComBat)$importance[2,2]*100,1),"%", sep=""),ylab=paste("PCA1: ", round(summary(pca_preComBat)$importance[2,1]*100,1),"%",sep=""), main = "Autores: Pre Combat PC1 vs PC2")
legend("topright", legend=pca_lgnd, col=c(1,1),pch = c(17,19), pt.cex = 1, cex=0.75, horiz = F)
text(pca_preComBat$rotation[,1], pca_preComBat$rotation[,2], labels = row.names(pca_preComBat$rotation), pos = 4,cex=0.5)


cor(ind_fpkm,ind_tpm,method="spearman")
cor(ind_fpkm,ind_tpm,method="pearson")



########################################################################################################



filter<-AR1[5]>2 | AR2[5]>2 |AR3[5]>2 |AR4[5]>2 |AR5[5]>2 | PC1[5]>2|PC2[5]>2|PC3[5]>2|PC4[5]>2|PC5[5]>2|PC6[5]>2|PC7[5]>2|PC8[5]>2|PC9[5]>2|PC10[5]>2 
resFrame<-resFrame[filter,]
filterres<-resFrame[2]>=2 | resFrame[2]<=-2
resFrame<-resFrame[filterres,]
filterp<-resFrame[6]<=0.05
resFrame<-resFrame[filterp,]
name<-row.names.data.frame(resFrame)
#name<-substring(name,1,15)
resFrame<-cbind(name,resFrame)

#matching con el indice para encontrar los genes
human_index<- rtracklayer::import('gencode.v30.chr_patch_hapl_scaff.annotation.gtf.gz')
human_df=as.data.frame(human_index)
name<-
gene_pos<-match(name,human_df$transcript_id)
gene_name<-human_df$gene_name[gene_pos]
gene_type<-human_df$gene_type[gene_pos]
resFrame<-cbind(gene_type,resFrame)
resFrame <- cbind(gene_name, resFrame)


#tpm frame
tpm_pos<-match(name,PC1$target_id)
tpm_aux<-cbind(AR1$tpm[tpm_pos], AR2$tpm[tpm_pos], AR3$tpm[tpm_pos], AR4$tpm[tpm_pos], AR5$tpm[tpm_pos],  PC1$tpm[tpm_pos], PC2$tpm[tpm_pos], PC3$tpm[tpm_pos], PC4$tpm[tpm_pos], PC5$tpm[tpm_pos], PC6$tpm[tpm_pos], PC7$tpm[tpm_pos], PC8$tpm[tpm_pos], PC9$tpm[tpm_pos],PC10$tpm[tpm_pos])
colnames(tpm_aux) <- c("AR_1", "AR_2","AR_3", "AR_4","AR_5","PC_1","PC_2","PC_3","PC_4","PC_5","PC_6","PC_7","PC_8","PC_9","PC_10")
rownames(tpm_aux)<-gene_name
tpm_frame<-as.data.frame(tpm_aux)
scale_tpm_frame<-t(apply(tpm_frame,1,function(x) scale(x)))
colnames(scale_tpm_frame) <- c("AR1", "AR2","AR3", "AR4","AR5","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
#funcion encontrada de internet

colors = c(seq(-2,-0.4,length=100),seq(-0.39,1.5,length=100),seq(1.51,4,length=100))
my_palette <- colorRampPalette(c("green", "white", "red"))(n = 299)


dev.off()
heatmap.2(scale_tpm_frame, dendrogram=c("column"), col=my_palette, breaks=colors, trace="none",scale="none",symbreaks=T, symm=F, symkey=F, density.info="none" )

write.table(resFrame$gene_name, sep="\t", quote=F, row.names=F, col.names=F)

#GSEA
full_gene_name<-toupper(row.names(resMat[[1]]))
pvals<-resMat[[1]]$pvalue
signs<-sign(resMat[[1]]$log2FoldChange)
pvals<- -log10(pvals)*signs
GSEO_frame<-cbind(resMat[[1]]$pvalue,full_gene_name,pvals)
GSEO_frame<-as.data.frame(GSEO_frame)
colnames(GSEO_frame)<-c("pvalue","Transcript_id","rankingvals")
gene_pos_GSEO<-match(full_gene_name,human_df$transcript_id)
gene_name_GSEO<-human_df$gene_name[gene_pos_GSEO]
gene_type_GSEO<-human_df$gene_type[gene_pos_GSEO]
GSEO_frame<-cbind(gene_type_GSEO,GSEO_frame)
GSEO_frame <- cbind(gene_name_GSEO, GSEO_frame)

top_GSEO<-GSEO_frame[order(GSEO_frame$rankingvals,decreasing = TRUE),c(1,5)]
top_GSEO<-top_GSEO[!duplicated(top_GSEO$gene_name_GSEO),]

mainDir<-"C:/Users/angel/Documents/Bioinformatica/AR/"
setwd(mainDir)
#gene_set<-"C:/Users/angel/Documents/Bioinformatica/GMT/c2.cp.v6.2.symbols.gmt,C:/Users/angel/Documents/Bioinformatica/GMT/h.all.v6.2.symbols.gmt,C:/Users/angel/Documents/Bioinformatica/GMT/c2.cp.reactome.v6.2.symbols.gmt,C:/Users/angel/Documents/Bioinformatica/GMT/c5.all.v6.2.symbols.gmt,C:/Users/angel/Documents/Bioinformatica/GMT/c2.cp.biocarta.v6.2.symbols.gmt,C:/Users/angel/Documents/Bioinformatica/GMT/c2.cp.kegg.v6.2.symbols.gmt"
gene_set<-"C:/Users/angel/Documents/Bioinformatica/GMT/c2.cp.v6.2.symbols.gmt"
outputDir<- "gsea_reports_c2_v6"
system(paste("mkdir", outputDir))

gsea_pos <- c()
gsea_neg <- c()
comparison_name <- "AR"
rnk_file <- paste(comparison_name, ".rnk", sep="")

write.table(top_GSEO, sep="\t", quote=F, row.names=F, col.names=F, file=rnk_file)
java_command <- "java -cp C:/Users/angel/Documents/Bioinformatica/GMT/gsea-3.0.jar -Xmx2000m xtools.gsea.GseaPreranked"
java_command <- paste(java_command, "-gmx", gene_set, "-collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk")
java_command <- paste(java_command, rnk_file, "-scoring_scheme weighted -rpt_label", comparison_name, "-include_only_symbols true")
java_command <- paste(java_command, "-make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 5 -zip_report false -out")
java_command <- paste(java_command, outputDir, "-gui false")
system(java_command)
setwd(outputDir)
aux_dir <- dir(pattern=comparison_name)
setwd(aux_dir)
aux_files <- dir()
aux_pos_file <- read.table(file=aux_files[grep("gsea_report_for_na_pos", aux_files)[2]], header=T, sep="\t", as.is=TRUE)
aux_neg_file <- read.table(file=aux_files[grep("gsea_report_for_na_neg", aux_files)[2]], header=T, sep="\t", as.is=TRUE)
if (nrow(aux_pos_file)>0)
  gsea_pos <- rbind(gsea_pos, cbind(aux_pos_file[,c("NAME","NES","NOM.p.val","FDR.q.val")]))
if (nrow(aux_neg_file)>0)
  gsea_neg <- rbind(gsea_neg, cbind(aux_neg_file[,c("NAME","NES","NOM.p.val","FDR.q.val")]))
setwd(paste(mainDir, outputDir, sep=""))
system(paste("rm -rf", aux_dir))
setwd(mainDir)
system(paste("rm", rnk_file))
rm(top_GSEO)
rm(aux_pos_file)
rm(aux_neg_file)	

colnames(gsea_pos)[ncol(gsea_pos)] = "Comparison"
colnames(gsea_neg)[ncol(gsea_neg)] = "Comparison"
gsea_pos[,ncol(gsea_pos)]=as.character(gsea_pos[,ncol(gsea_pos)])
gsea_neg[,ncol(gsea_neg)]=as.character(gsea_neg[,ncol(gsea_neg)])

save(gsea_pos, gsea_neg, file="gsea_res.Rdata")



#hipergeometrica
geneNames<-toupper(GSEO_frame$gene_name_GSEO )  # 188753 
gene_sets=readLines("C:/Users/angel/Documents/Bioinformatica/GMT/c5.all.v6.2.symbols.gmt")            # 5917



geneSet_genes=list()
geneSet_names=c()

for (n in 1:length(gene_sets)){
  geneSet_names=c(geneSet_names,strsplit(strsplit(gene_sets[n], "\t")[[1]][1], "\\(")[[1]][1])
  aux_genes=strsplit(gene_sets[n], "\t")[[1]][3:length(strsplit(gene_sets[n], "\t")[[1]])]
  geneSet_genes[[n]]= intersect(aux_genes, geneNames)
}
geneSet_size=unlist(lapply(geneSet_genes, length))
geneSet_genes= geneSet_genes[-which(geneSet_size < 5)]
length(geneSet_genes)          # 5911
geneSet_names= geneSet_names[-which(geneSet_size < 5)]
names(geneSet_genes) = geneSet_names
all_genes = toupper(GSEO_frame$gene_name_GSEO[which(GSEO_frame$gene_type_GSEO=="protein_coding")])   # 166735



up_degs = resFrame$gene_name[which(resFrame$log2FoldChange > 0)]
down_degs = resFrame$gene_name[which(resFrame$log2FoldChange < 0)]



gene_set_res=list()
#significant_proteins = toupper(up_degs)
significant_proteins = toupper(down_degs)
geneSets= c()
genesInSet = list()
for(g in 1:length(geneSet_genes)) {
  geneSets= rbind(geneSets, p.overlap.gene_set(all_genes, geneSet_genes[[g]], significant_proteins))
  genesInSet[[names(geneSet_genes)[g]]] = significant_proteins[which(significant_proteins %in% geneSet_genes[[g]])]
}
rownames(geneSets) = names(geneSet_genes)
geneSets = cbind(geneSets, p.adjust(geneSets[,1]))
colnames(geneSets) = c("p.value", "Num", "FDR")
aux_genesInSet=unlist(lapply(genesInSet, function(x) ifelse(length(x) > 0, paste(x, collapse=","), 0)))
geneSets = data.frame(geneSets, genes=aux_genesInSet)
geneSets= geneSets[order(geneSets[,1]),]
geneSets= geneSets[order(geneSets[,2], decreasing=T),]
geneSets = geneSets[,c(1,3,2,4)]
#gene_set_res[["up"]]= geneSets
gene_set_res[["down"]]= geneSets
p.overlap.gene_set <- function(total_genes, geneset_genes, significant_genes, lower.tail=FALSE) {
  total_genes_in_set <- length(which(total_genes %in% geneset_genes))
  total_genes_not_in_set <- length(total_genes) - total_genes_in_set
  significant_genes_in_set <- length(which(significant_genes %in% geneset_genes))
  c(phyper(significant_genes_in_set, total_genes_in_set, total_genes_not_in_set, length(significant_genes), lower.tail=lower.tail), significant_genes_in_set)
}


newnames<-as.character(resFrame$gene_name)
newpos<-match(newnames,fpkm$Gene_ID)
newfpkm<-fpkm[newpos,1:17]
newfpkm<-data.frame(row.names=newnames)
#de esta forma no existen algunos de los genes encontrados en la otra manera. ya que provienen de ENSG y no de ENST
