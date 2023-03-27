

# limma ------------------------------------------------------------------------
library(limma)
group <- factor(list_ID$type,levels = c("COVID","HD"))

limma_design<-model.matrix(~-1+group)

rownames(limma_design) <- list_ID$ID

limma_expr <- view_bulkRNAseq

limma_expr <- log2(limma_expr+1)

contrast.matrix<-makeContrasts(contrasts = "groupCOVID-groupHD",
                               levels = limma_design)  

fit <- lmFit(limma_expr,limma_design)   

fit1 <- contrasts.fit(fit, contrast.matrix)  

fit2 <- eBayes(fit1,0.01)    

tempOutput = topTable(fit2, 
                      coef="groupCOVID-groupHD", 
                      n=nrow(fit2),
                      lfc=log2(1),
                      adjust="fdr")    

dif<-tempOutput[tempOutput[,"adj.P.Val"]<0.05,]
sd(tempOutput$logFC)

dif_logFC <- dif[abs(dif[,"logFC"])>sd(tempOutput$logFC)*3,]
# VN ----------------------------------------------------------------------
library(VennDiagram)
venn.diagram(list(RNAseq = rownames(dif_logFC_bulk),
                  Protein = rownames(dif_protein)),
             height = 5000,
             width = 5000,
             resolution = 300, 
             imagetype = "tiff", 
             alpha=c(0.5,0.5),
             fill=c("red","blue"), 
             cat.fontface="bold",
             fontfamily=5,
             cat.cex = 3,
             main="Venn diagram of RNAseq and protein",
             main.cex = 2, 
             main.fontface = 2, 
             main.fontfamily = 3,
             cex = 3,
             filename = "result/dif/VennDiagram2.tif")

list_gene <- merge(Panepitranscriptome_gene_list,
                   list,
                   by = 'Row.names')
rownames(list_gene) <- list_gene$Row.names
dif_logFC_bulk$id <- rownames(dif_logFC_bulk)
dif_protein$id <- rownames(dif_protein)

merge_dif <- merge(dif_logFC_bulk,dif_protein,by = 'id')
merge_dif

list_merge_bulkRNAseq_id <- merge(merge_dif,
                                  raw_bulkRNAseq,
                                  by.x = 'id',
                                  by.y = 'Associated_Gene_Name')
list_merge_protein_id <- merge(merge_dif,
                               raw_protein,
                               by.x = 'id',
                               by.y = 'Gene_name')
merge_dif <- c(rownames(dif_logFC_bulk),rownames(dif_protein))
merge_dif <- unique(merge_dif)
merge_dif <- data.frame(id =merge_dif)

