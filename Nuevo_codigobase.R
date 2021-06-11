library(tximport)
library(rtracklayer)
library(DESeq2)

loadabundancedata = function(dir, filetype) { #dir = directorio donde se encuentran tus archivos, filetype = "extension del archivo"
  files = list.files(path = dir, pattern=filetype)
  temp_list <- list()
  for (i in 1:length(files)) {#Guarda las tablas de datos en una lista, ademas de cortar el nombre para quitar ".tsv"
    temp_list[[i]] <- assign(substr(files[i], 1,nchar(files[i])-4),read.table(files[i], sep = "\t", header=TRUE))
    names(temp_list)[i] <- c(substr(files[i], 1,nchar(files[i])-4))
  }
  return(temp_list)
}

dir <- "D:/Archivos_de_trabajo/Diabetes_Esclerosis"
setwd(dir)
abundancias_ED <- loadabundancedata(dir, ".tsv")

genes <- as.data.frame(rtracklayer::import('D:/Archivos_de_trabajo/Analisis_combinado/Homo_sapiens.GRCh38.96.chr_patch_hapl_scaff.gtf.gz'))

tx2gene_transcriptid <- substr(abundancias_ED$ED1$target_id,1,15)
tx2gene_geneid <- genes$gene_id[match(tx2gene_transcriptid, genes$transcript_id)]

tx2gene_df <- data.frame(tx2gene_transcriptid, tx2gene_geneid)

files <- list.files(dir, pattern = ".tsv")
txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = tx2gene_df, ignoreTxVersion = TRUE)

colnames(txi.kallisto.tsv$counts) <- substr(files, 1,nchar(files)-4)

sampleTable <- data.frame(condition = factor(substr(files, 1,1)))
rownames(sampleTable) <- colnames(txi.kallisto.tsv$counts)
dds <- DESeqDataSetFromTximport(txi.kallisto.tsv, sampleTable, ~condition)

dds <- DESeq(dds)
dds$condition <- relevel(dds$condition, ref = "E")
dds <- DESeq(dds)
res <- as.data.frame(results(dds))
res$log2FoldChange = sapply(res$log2FoldChange, function(x) ifelse(is.na(x), 0, x))
res$pvalue = sapply(res$pvalue, function(x) ifelse(is.na(x), 1, x))
res$padj = sapply(res$padj, function(x) ifelse(is.na(x), 1, x))

genename_match <- match(rownames(res), genes$gene_id)
res <- data.frame(gene_name = genes$gene_name[genename_match], res)

res_filt <- res[abs(res$log2FoldChange) > 1,]
res_filt <- res_filt[res_filt$padj < .05,]
