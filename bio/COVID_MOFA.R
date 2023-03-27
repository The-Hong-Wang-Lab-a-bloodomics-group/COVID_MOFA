
# 1 package  -------------------------------------------------------------------
library(data.table)
library(MOFAdata)
library(MOFA2)
library(grid)
library(gplots)
library(ggplot2)
library(readxl)
library(readr)
# 2 read data ------------------------------------------------------------------
# read and clean data 
## 2.1 Platelet transcriptome --------------------------------------------------
raw_bulkRNAseq <- read_excel("data/bulkRNAseq_protein_coding.xlsx")
data_bulkRNAseq <- aggregate(raw_bulkRNAseq,
                             by = list(raw_bulkRNAseq$Associated_Gene_Name),
                             mean)
rownames(data_bulkRNAseq) <- data_bulkRNAseq$Group.1
data_bulkRNAseq <- subset(data_bulkRNAseq,
                          select = -c(`...1`,
                                      `...96`,
                                      `Group.1`,
                                      Ensembl_Gene_ID))
data_bulkRNAseq <- subset(data_bulkRNAseq,
                          select = -c(`avr`,
                                      `Associated_Gene_Name`,
                                      `Gene_Biotype`,
                                      Position,
                                      Gene_Description))

## 2.2 Platelet proteome -------------------------------------------------------
raw_protein <- read_csv("data/protein.csv")
data_protein <- aggregate(raw_protein,
                          by = list(raw_protein$Gene_name),
                          mean)
rownames(data_protein) <- data_protein$Group.1
data_protein <- subset(data_protein,
                       select = -c(`Group.1`,
                                   Uniprot_ID,
                                   Gene_name))

## 2.3 Plasma proteomics -------------------------------------------------------
raw_Plasma_proteomics <- read_excel("data/Proteomics.xlsx")
View(raw_Plasma_proteomics)
data_Plasma_proteomics <- aggregate(raw_Plasma_proteomics,
                                    by = list(raw_Plasma_proteomics$`Gene Symbol`),
                                    mean)
rownames(data_Plasma_proteomics) <- data_Plasma_proteomics$Group.1
data_Plasma_proteomics <- subset(data_Plasma_proteomics,
                                 select = -c(`Group.1`,
                                             Proteins,
                                             `Gene Symbol`))

## 2.4 Plasma metabolomic ------------------------------------------------------
raw_Metabolomic <- read_excel("data/Metabolomic.xlsx")

data_Metabolomic <- aggregate(raw_Metabolomic,
                              by = list(raw_Metabolomic$Name),
                              mean)
# replace unrecognized characters
raw_Metabolomic$Name <- gsub('¡À',' ',raw_Metabolomic$Name)
raw_Metabolomic$Name <- gsub('¦Á','Alpha',raw_Metabolomic$Name)
raw_Metabolomic$Name<- gsub('¦Ä','Delta',raw_Metabolomic$Name)
raw_Metabolomic$Name <- gsub('¦Â','Beta',raw_Metabolomic$Name)
raw_Metabolomic$Name <- gsub('¦Ø','Omega',raw_Metabolomic$Name)
rownames(data_Metabolomic) <- data_Metabolomic$Group.1
data_Metabolomic <- subset(data_Metabolomic,
                           select = -c(`Group.1`,
                                       Compound_ID,
                                       Name,
                                       HMDB_ID,
                                       `SuperClass(HMDB)`))

# 3 ID conversion --------------------------------------------------------------

# 3.1 merge sampleID -----------------------------------------------------------

list_metadata <- read_excel("../../data/list_metadata.xlsx")
list_ID <- subset(list_metadata,
                  select = c(`ÐÕÃû`,
                             bulk_ID,
                             protein_ID,
                             ID))
list_protein_ID <- colnames(data_protein)
list_protein_ID <- as.data.frame(list_protein_ID)
colnames(list_protein_ID)[1] <- c('protein_ID')
list_protein_ID <- merge(list_ID,list_protein_ID,by = 'protein_ID')
list_bulk_ID <- colnames(data_bulkRNAseq)
list_bulk_ID <- as.data.frame(list_bulk_ID)
colnames(list_bulk_ID) <- c('bulk_ID')
list_bulk_ID <- merge(list_ID,list_bulk_ID,by = 'bulk_ID')
list_ID <- merge(list_bulk_ID,list_protein_ID,by ='ID',all = T)
## 3.2 Platelet transcriptome --------------------------------------------------

data_bulkRNAseq <- as.data.frame(t(data_bulkRNAseq))
data_bulkRNAseq$bulk_ID <- rownames(data_bulkRNAseq)
view_bulkRNAseq <- merge(list_ID,
                         data_bulkRNAseq,
                         by = 'bulk_ID',
                         all.x =T)
data_bulkRNAseq <- subset(data_bulkRNAseq,
                          select = -c(bulk_ID))
data_bulkRNAseq <- as.data.frame(t(data_bulkRNAseq))
rownames(view_bulkRNAseq) <- view_bulkRNAseq$ID
view_bulkRNAseq <- view_bulkRNAseq[,-1:-4]
view_bulkRNAseq <- as.data.frame(t(view_bulkRNAseq))

## 3.3 Platelet proteome -------------------------------------------------------
data_protein <- as.data.frame(t(data_protein))
data_protein$protein_ID <- rownames(data_protein)
view_protein <- merge(list_ID,
                      data_protein,
                      by = 'protein_ID',
                      all.x =T)
rownames(view_protein) <- view_protein$ID
data_protein <- subset(data_protein,
                       select = -c(protein_ID))
data_protein <- as.data.frame(t(data_protein))
view_protein <- view_protein[,-1:-4]
view_protein <- as.data.frame(t(view_protein))

## 3.4 Plasma proteomics -------------------------------------------------------
data_Plasma_proteomics <- as.data.frame(t(data_Plasma_proteomics))
data_Plasma_proteomics$protein_ID <- gsub("HD_",
                                          "HD-",
                                          rownames(data_Plasma_proteomics))
view_Plasma_proteomics <- merge(list_ID,
                                data_Plasma_proteomics,
                                by = 'protein_ID',
                                all.x =T)
data_Plasma_proteomics <- subset(data_Plasma_proteomics,
                                 select = -c(protein_ID))
data_Plasma_proteomics <- as.data.frame(t(data_Plasma_proteomics))
rownames(view_Plasma_proteomics) <- view_Plasma_proteomics$ID
view_Plasma_proteomics <- view_Plasma_proteomics[,-1:-4]
view_Plasma_proteomics <- as.data.frame(t(view_Plasma_proteomics))

## 3.5 Plasma metabolomic ------------------------------------------------------
data_Metabolomic <- as.data.frame(t(data_Metabolomic))
data_Metabolomic$protein_ID <- gsub("HD_","HD-",rownames(data_Metabolomic))
view_Metabolomic <- merge(list_ID,
                          data_Metabolomic,
                          by = 'protein_ID',
                          all.x =T)
data_Metabolomic <- subset(data_Metabolomic,
                           select = -c(protein_ID))
data_Metabolomic <- as.data.frame(t(data_Metabolomic))
rownames(view_Metabolomic) <- view_Metabolomic$ID
view_Metabolomic <- view_Metabolomic[,-1:-4]
view_Metabolomic <- as.data.frame(t(view_Metabolomic))

# 4 data filter ----------------------------------------------------------------
# filter the result of bio/MOFA_limma.R 
view_bulkRNAseq_filter <- as.data.frame(t(subset(as.data.frame(t(view_bulkRNAseq)),
                                                 select = list_merge_bulkRNAseq_id$id)))
view_protein_filter <- as.data.frame(t(subset(as.data.frame(t(view_protein)),
                                              select = list_merge_protein_id$id)))
# sort the sample ID
view_bulkRNAseq_filter <- subset(view_bulkRNAseq_filter,
                                 select = list_metadata$ID)

view_protein_filter <- subset(view_protein_filter,
                              select = list_metadata$ID)

view_Plasma_proteomics_filter <- subset(view_Plasma_proteomics,
                                        select = list_metadata$ID)

view_Metabolomic_filter <- subset(view_Metabolomic,
                                  select = list_metadata$ID)

# 5 data standardization -------------------------------------------------------
view_bulkRNAseq_filter_scale <- scale(view_bulkRNAseq_filter)
view_protein_filter_scale <- scale(view_protein_filter)
view_Plasma_proteomics_filter_scale <- scale(view_Plasma_proteomics_filter)
view_Metabolomic_filter_scale <- scale(view_Metabolomic_filter)

# 6 Ensembl ID data(GSEA need) -------------------------------------------------

## 6.1 Platelet transcriptome --------------------------------------------------
view_bulkRNAseq_filter_scale_Ensembl <- view_bulkRNAseq_filter_scale
view_bulkRNAseq_filter_scale_Ensembl <- as.data.frame(view_bulkRNAseq_filter_scale_Ensembl)
view_bulkRNAseq_filter_scale_Ensembl$Associated_Gene_Name <- rownames(view_bulkRNAseq_filter_scale_Ensembl)

view_bulkRNAseq_filter_scale_Ensembl <- merge(view_bulkRNAseq_filter_scale_Ensembl,
                                           list_RNAID_mapping,
                                           by = 'Associated_Gene_Name',
                                           all.x = T)
write.csv(view_bulkRNAseq_filter_scale_Ensembl,
          file = 'data_input/view_bulkRNAseq_filter_scale_Ensembl.csv')
# Some ids are not matched. need to manually match them
view_bulkRNAseq_filter_scale_Ensembl <- read_csv("data_input/view_bulkRNAseq_filter_scale_Ensembl.csv")
view_bulkRNAseq_filter_scale_Ensembl <- as.data.frame(view_bulkRNAseq_filter_scale_Ensembl)
rownames(view_bulkRNAseq_filter_scale_Ensembl) <- view_bulkRNAseq_filter_scale_Ensembl$Ensembl
view_bulkRNAseq_filter_scale_Ensembl <- subset(view_bulkRNAseq_filter_scale_Ensembl,
                                               select = -c(`...1`, 
                                                           Associated_Gene_Name, 
                                                           Ensembl))

## 6.2  Platelet proteome ------------------------------------------------------
list_RNAID_mapping <- subset(raw_bulkRNAseq,
                             select = c(Associated_Gene_Name,
                                        Ensembl))
view_protein_filter_scale_Ensembl <- view_protein_filter_scale
view_protein_filter_scale_Ensembl <- as.data.frame(view_protein_filter_scale_Ensembl)
view_protein_filter_scale_Ensembl$Associated_Gene_Name <- rownames(view_protein_filter_scale_Ensembl)

view_protein_filter_scale_Ensembl <- merge(view_protein_filter_scale_Ensembl,
                                           list_RNAID_mapping,
                                           by = 'Associated_Gene_Name',
                                           all.x = T)
# Some ids are not matched. need to manually match them
write.csv(view_protein_filter_scale_Ensembl,file = 'data_input/view_protein_filter_scale_Ensembl.csv')

view_protein_filter_scale_Ensembl <- read_csv("data_input/view_protein_filter_scale_Ensembl.csv")
view_protein_filter_scale_Ensembl <- as.data.frame(view_protein_filter_scale_Ensembl)
rownames(view_protein_filter_scale_Ensembl) <- view_protein_filter_scale_Ensembl$Ensembl
view_protein_filter_scale_Ensembl <- subset(view_protein_filter_scale_Ensembl,
                                            select = -c(`...1`,
                                                        Associated_Gene_Name,
                                                        Ensembl))

## 6.3 Plasma_proteomics -------------------------------------------------------
view_Plasma_proteomics_filter_scale_Ensembl <- view_Plasma_proteomics_filter_scale
view_Plasma_proteomics_filter_scale_Ensembl <- as.data.frame(view_Plasma_proteomics_filter_scale_Ensembl)
view_Plasma_proteomics_filter_scale_Ensembl$Associated_Gene_Name <- rownames(view_Plasma_proteomics_filter_scale_Ensembl)
view_Plasma_proteomics_filter_scale_Ensembl <- merge(view_Plasma_proteomics_filter_scale_Ensembl,
                                           list_RNAID_mapping,
                                           by = 'Associated_Gene_Name',
                                           all.x = T)
write.csv(view_Plasma_proteomics_filter_scale_Ensembl,
          file = 'data_input/view_Plasma_proteomics_filter_scale_Ensembl.csv')
# Some ids are not matched. need to manually match them
view_Plasma_proteomics_filter_scale_Ensembl <- read_csv("data_input/view_Plasma_proteomics_filter_scale_Ensembl.csv")
view_Plasma_proteomics_filter_scale_Ensembl <- as.data.frame(view_Plasma_proteomics_filter_scale_Ensembl)
rownames(view_Plasma_proteomics_filter_scale_Ensembl) <- view_Plasma_proteomics_filter_scale_Ensembl$Ensembl
view_Plasma_proteomics_filter_scale_Ensembl <- subset(view_Plasma_proteomics_filter_scale_Ensembl,
                                                      select = -c(`...1`,
                                                                  Associated_Gene_Name,
                                                                  Ensembl))


# 7 data input -----------------------------------------------------------------
data <- list(as.matrix(view_bulkRNAseq_filter_scale),
             as.matrix(view_protein_filter_scale),
             as.matrix(view_Plasma_proteomics_filter_scale),
             as.matrix(view_Metabolomic_filter_scale))

names(data) <- c('Platelet transcriptome',
                 'Platelet proteome',
                 'Plasma proteome',
                 'Plasma metabolomic')
lapply(data,dim)

# data heatmap
data_heatmap <- view_bulkRNAseq_filter_scale
ComplexHeatmap::Heatmap(data_heatmap,
                        show_heatmap_legend = F,
                        na_col = "gray",
                        name = "value",
                        show_column_names = F,
                        show_row_names = F,
                        cluster_rows = F,
                        cluster_columns = F,
                        row_names_gp = grid::gpar(ngle=45,fontsize = 12),
                        column_names_gp = grid::gpar(fontsize = 12),
                        col = colorRamp2(c(quantile(data_heatmap,na.rm = T)[2],
                                           quantile(data_heatmap,na.rm = T)[3]), 
                                         c("#3F66E1FF","#ED1E1AFF")))

# 8 Model initialization -------------------------------------------------------
MOFAobject <- create_mofa(data)
plot_data_overview(MOFAobject)

## 8.1 data summary ------------------------------------------------------------
plot_data_overview(MOFAobject)
print(MOFAobject)

## 8.2 set options -------------------------------------------------------------

### 8.2.1 data options ---------------------------------------------------------
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

### 8.2.2 model options --------------------------------------------------------
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15
head(model_opts)

### 8.2.3 training options -----------------------------------------------------
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
head(train_opts)

## 8.3 prepare for MOFA model --------------------------------------------------
MOFAobject_data <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

## 8.4 loading python environment ----------------------------------------------

#library(reticulate)
#use_python()

# 9 run MOFA -------------------------------------------------------------------
MOFAobject.trained <- run_mofa(MOFAobject_data,
                               outfile = file.path(getwd(),
                                                   "data/model_TPM_symbol_102.hdf5"))

# 10 model summary -------------------------------------------------------------

# Overview of the trained MOFA model 
slotNames(MOFAobject.trained)
names(MOFAobject.trained@data)
dim(MOFAobject.trained@data$RNAseq)
names(MOFAobject.trained@expectations)

## 10.1 Add sample metadata to the model ---------------------------------------

# Check whether the sample ID is corresponding
# the type of metadata must be data.frame
stopifnot(all(sort(list_metadata$ID) == sort(unlist(samples_names(MOFAobject.trained)))))
list_metadata$`Clinical spectrum` <- factor(list_metadata$`Clinical spectrum`,
                                            levels = c('Non-omicron',
                                                       'Moderate',
                                                       'Mild'))
list_metadata$`¸´Ñô£¨1Ñô0Òõ£©CT£¼40` <- factor(list_metadata$`¸´Ñô£¨1Ñô0Òõ£©CT£¼40`,
                                        levels = c('Yes',
                                                   'No'))
list_metadata$`¸´Ñô£¨1Ñô0Òõ£©CT£¼35` <- factor(list_metadata$`¸´Ñô£¨1Ñô0Òõ£©CT£¼35`,
                                        levels = c('Yes',
                                                   'No'))
samples_metadata(MOFAobject.trained) <- list_metadata

## 10.2 Correlation between factors --------------------------------------------

# the good sanity check is no correlation by factors.
# MOFA does not have orthogonal constraints like PCA.
# if there is a strong correlation between factors, 
# it means that the data may be incorrectly normalized 
# or too many factors are used.

plot_factor_cor(MOFAobject.trained)

## 10.3 Plot variance decomposition --------------------------------------------
# MOFAobject.trained@cache[["variance_explained"]]
# variance explained--per factor
plot_variance_explained(MOFAobject.trained,max_r2=10,plot_total =T)[[1]]

# variance explained--total
plot_variance_explained(MOFAobject.trained,max_r2=10,plot_total =T)[[2]]

# this shows whether the factors can provide a good fit to the data
# if the data has a strong nonlineary noisy 
# the result will only have a small amount of variance explained(<10%) 

# 11 Characterisation of Factor  -----------------------------------------------

## 11.1 Association analysis ---------------------------------------------------
# the value of covariates must be in colnames(metadata)
# the function visualize it by corrplot and pheatmap.
# the correlation coefficient is pearson's coefficient.
# if the logical value of return_data if TURE, 
# the function will return the result.
# after that, you can visualize it however you want.
# if you need the analysis of correlation, you can use MOFA_factor_metadata_correlation.R in bio/
covariates <- c("Gender",
                "Age",
                "Omicron phases",
                "Clinical spectrum",
                "Ct value",
                "PLT num",
                "Vaccine type",
                "Underling medical conditions")
metadata_cor <- correlate_factors_with_covariates(MOFAobject.trained, 
                                  covariates = covariates,
                                  tl.cex = 0.7,
                                  plot=c("r"),
                                  return_data = T)
metadata_log_pval <- correlate_factors_with_covariates(MOFAobject.trained, 
                                                  covariates = covariates,
                                                  tl.cex = 0.7,
                                                  plot=c("log_pval"),
                                                  return_data = T)
correlate_factors_with_covariates(MOFAobject.trained,
                                  covariates = covariates, 
                                  tl.cex = 0.7, 
                                  plot=c("log_pval"))


## 11.2 the analysis of correlation between factors and metadata ---------------

### 11.2.1 as a whole ----------------------------------------------------------
# if you want to use pearson's coefficient, 
# you can use metadata_log_pval and metadata_cor
# if you want to use spearman's coefficient, 
# you need use bio/MOFA_factor_metadata_correlation.R to get cor_p and cor_r
ComplexHeatmap::Heatmap(cor_p,
                        name = "log10 pvalue",
                        cluster_rows = F,
                        cluster_columns = F,
                        row_names_gp = grid::gpar(ngle=45,fontsize = 12),
                        column_names_gp = grid::gpar(fontsize = 12),
                        col = colorpanel(50,low="white", high="red"),
                        rect_gp = gpar(col = "black", lwd = 1))
ComplexHeatmap::Heatmap(cor_r,
                        name = "cor value",
                        cluster_rows = F,
                        cluster_columns = F,
                        row_names_gp = grid::gpar(ngle=45,fontsize = 12),
                        column_names_gp = grid::gpar(fontsize = 12),
                        col = colorRamp2(c(-1, 0, 1), 
                                         c("#2727F6", "white", "red")),
                        rect_gp = gpar(col = "black", lwd = 1),)
heatmap(metadata_log_pval,
        Colv = NA,
        Rowv = NA)
pheatmap::pheatmap(metadata_log_pval)

## Warning in correlate_factors_with_covariates(MOFAobject, covariates =
## c("Gender", : There are non-numeric values in the covariates data.frame,
## converting to numeric...

### 11.2.2 a single factor -----------------------------------------------------
# use bio/MOFA_factor_Omicron_phases_scatter.R

## 11.3 Plot factor values -----------------------------------------------------
for(i in 1:15){
  plot_list[[i]] <- plot_factor(MOFAobject.trained,dot_size = 1.5,
                                factors = i, 
                                color_by = paste0("Factor",i)
                                )
}
prow <- plot_grid(plot_list[[1]]+theme(legend.position="none"),
                  plot_list[[2]]+theme(legend.position="none"),
                  plot_list[[3]]+theme(legend.position="none"),
                  plot_list[[4]]+theme(legend.position="none"),
                  plot_list[[5]]+theme(legend.position="none"),
                  plot_list[[6]]+theme(legend.position="none"),
                  plot_list[[7]]+theme(legend.position="none"),
                  plot_list[[8]]+theme(legend.position="none"),
                  plot_list[[9]]+theme(legend.position="none"),
                  plot_list[[10]]+theme(legend.position="none"),
                  plot_list[[11]]+theme(legend.position="none"),
                  plot_list[[12]]+theme(legend.position="none"),
                  plot_list[[13]]+theme(legend.position="none"),
                  plot_list[[14]]+theme(legend.position="none"),
                  plot_list[[15]]+theme(legend.position="none"),
                  ncol = 5,
                  nrow = 3,
                  align = "v")
prow

## 11.4 Plot feature weights ---------------------------------------------------

plot_weights(MOFAobject.trained,
             view = "Platelet transcriptome",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)
plot_weights(MOFAobject.trained,
             view = "protein",
             factor = 1,
             nfeatures = 10,     
             scale = T           
)
plot_weights(MOFAobject.trained,
             view = "Plasma proteomics",
             factor = 1,
             nfeatures = 10,     
             scale = T           
)
plot_weights(MOFAobject.trained,
             view = "Plasma metabolomic",
             factor = 1,
             nfeatures = 10,     
             scale = T           
)
plot_top_weights(MOFAobject.trained,
                 view = "Platelet transcriptome",
                 factor = 1:5,
                 nfeatures = 10,
                 scale = T
)
plot_top_weights(MOFAobject.trained,
                 view = "protein",
                 factor = 11:15,
                 nfeatures = 10, 
                 scale = T
)
plot_top_weights(MOFAobject.trained,
                 view = "Plasma proteomics",
                 factor = 1:5,
                 nfeatures = 10,
                 scale = T
)
plot_top_weights(MOFAobject.trained,
                 view = "Plasma metabolomic",
                 factor = 1:5,
                 nfeatures = 10,
                 scale = T
)
# 12 other analysis ------------------------------------------------------------
# Plot molecular signatures in the input data

## 12.1 scatter plot -----------------------------------------------------------
# check whether there is a linear relationship beteen factors and metadata.
plot_data_scatter(MOFAobject.trained, 
                  view = "Platelet transcriptome",
                  factor = 1,  
                  features = 5,
                  sign = "positive",
                  color_by = "Omicron phases")
plot_data_scatter(MOFAobject.trained, 
                  view = "Platelet transcriptome",
                  factor = 1,  
                  features = 5,
                  sign = "negative",
                  color_by = "Omicron phases")
## 12.2 heatmap ----------------------------------------------------------------
# check the distribution of factor values
plot_data_heatmap(MOFAobject.trained, 
                  view = "Platelet transcriptome",
                  factor = 1,  
                  features = 50,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
                  )

## 12.3 Inspection of combinations of Factors ----------------------------------
# instance
p <- plot_factors(MOFAobject.trained, 
                  factors = c(1,3), 
                  dot_size = 2.5,
                  show_missing = T)
p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")
print(p)

# 13 GSEA ----------------------------------------------------------------------
# All GSEA results are obtained in MOFA_GSEA
# the resulting visualization is performed by bio/GSEA_dotplot.R
library(PCGSE)
utils::data(reactomeGS)

head(colnames(reactomeGS))

# GSEA on positive weights, with default options
res.positive_Plasma_proteomics <- run_enrichment(MOFAobject.trained, 
                               feature.sets = reactomeGS, 
                               view = "Plasma proteomics",
                               sign = "positive"
)
res.negative_Plasma_proteomics <- run_enrichment(MOFAobject.trained, 
                               feature.sets = reactomeGS, 
                               view = "Plasma proteomics",
                               sign = "negative"
)
plot_enrichment_heatmap(res.positive_protein)
plot_enrichment_heatmap(res.negative_protein)
## 13.1 Plasma_proteomics ------------------------------------------------------
for(i in 1:15){
  tryCatch({
    a <- plot_enrichment(res.positive_Plasma_proteomics,
                         factor = i,
                         max.pathways = 10)
    ggsave(filename = paste0("GSEA_Plasma_proteomics_factor",i,"_positive.pdf"),
           a,
           width = 5,
           height = 3,
           dpi = 300)
    rm(a)
  },
  error=function(e){cat("ERROR :",i,"positive",conditionMessage(e), "\n")})
  tryCatch({
    b <- plot_enrichment(res.negative_Plasma_proteomics, 
                         factor = i,
                         max.pathways = 10)
    ggsave(filename = paste0("GSEA_Plasma_proteomics_factor", i, "_negative.pdf"),
           b,
           width = 5,
           height = 3,
           dpi = 300)
    rm(b)
  }, 
  error=function(e){cat("ERROR :",i,"negative",conditionMessage(e), "\n")})
}
res.positive_protein <- run_enrichment(MOFAobject.trained, 
                                       feature.sets = reactomeGS, 
                                       view = "protein",
                                       sign = "positive")
res.negative_protein <- run_enrichment(MOFAobject.trained, 
                                       feature.sets = reactomeGS, 
                                       view = "protein",
                                       sign = "negative")

## 13.2 Platelet proteome ------------------------------------------------------
for(i in 1:15){
  tryCatch({
    a <- plot_enrichment(res.positive_protein, 
                         factor = i, 
                         max.pathways = 10)
    ggsave(filename = paste0("GSEA_protein_factor",i,"_positive.pdf"),
           a,
           width = 5,
           height = 3,
           dpi = 300)
    rm(a)
  }, error=function(e){cat("ERROR :",i,"positive",conditionMessage(e), "\n")})
  tryCatch({
    b <- plot_enrichment(res.negative_protein, 
                         factor = i, 
                         max.pathways = 10)
    ggsave(filename = paste0("GSEA_protein_factor",i,"_negative.pdf"),
           b,
           width = 5,
           height = 3,
           dpi = 300)
    rm(b)
  }, error=function(e){cat("ERROR :",i,"negative",conditionMessage(e), "\n")})
}
# RNAseq
res.positive_RNAseq <- run_enrichment(MOFAobject.trained, 
                                      feature.sets = reactomeGS, 
                                      view = "RNAseq",
                                      sign = "positive")
res.negative_RNAseq <- run_enrichment(MOFAobject.trained, 
                                      feature.sets = reactomeGS, 
                                      view = "RNAseq",
                                      sign = "negative")
## 13.3 Plasma proteome --------------------------------------------------------
for(i in 1:15){
  tryCatch({
    a <- plot_enrichment(res.positive_Plasma_proteomics, 
                         factor = i, 
                         max.pathways = 10)
    ggsave(filename = paste0("GSEA_Plasma_RNAseq",i,"_positive.pdf"),
           a,
           width = 5,
           height = 3,
           dpi = 300)
    rm(a)
  }, 
  error=function(e){cat("ERROR :",i,"positive",conditionMessage(e), "\n")})
  tryCatch({
    b <- plot_enrichment(res.negative_RNAseq, 
                         factor = i, 
                         max.pathways = 10)
    ggsave(filename = paste0("GSEA_Plasma_RNAseq",i,"_negative.pdf"),
           b,
           width = 5,
           height = 3,
           dpi = 300)
    rm(b)
  }, error=function(e){cat("ERROR :",i,"negative",conditionMessage(e), "\n")})
}
# 14 export MOFA results -------------------------------------------------------
names(MOFAobject.trained@data)
MOFAobject_factors <- get_factors(MOFAobject.trained, 
                       views = "all", 
                       factors = "all", 
                       as.data.frame = TRUE)
MOFAobject_weights <- get_weights(MOFAobject.trained, 
                       views = "Plasma metabolomic", 
                       factors = "all", 
                       as.data.frame = TRUE)
MOFAobject_weights_Metabolomic <- MOFAobject_weights
MOFAobject_weights_Metabolomic$Name <- MOFAobject_weights_Metabolomic$feature
list_Metaboiomic_mapping <- raw_Metabolomic[,1:4]
MOFAobject_weights_Metabolomic <- merge(list_Metaboiomic_mapping,
                                        MOFAobject_weights_Metabolomic,
                                        by = "Name",
                                        all.y = T)
head(MOFAobject_weights)
MOFAobject_factors <- get_factors(MOFAobject.trained,
                                  as.data.frame = TRUE)
weights <- MOFAobject_weights
weights$feature <- gsub('_RNAseq','',weights$feature)
table(weights$factor)
weights_1 <- subset(weights,weights$factor == "Factor1")
weights_2 <- subset(weights,weights$factor == "Factor2")
weights_3 <- subset(weights,weights$factor == "Factor3")
weights_4 <- subset(weights,weights$factor == "Factor4")
weights_5 <- subset(weights,weights$factor == "Factor5")
weights_6 <- subset(weights,weights$factor == "Factor6")
weights_7 <- subset(weights,weights$factor == "Factor7")
weights_8 <- subset(weights,weights$factor == "Factor8")
weights_9 <- subset(weights,weights$factor == "Factor9")
weights_10 <- subset(weights,weights$factor == "Factor10")
weights_11 <- subset(weights,weights$factor == "Factor11")
weights_12 <- subset(weights,weights$factor == "Factor12")
weights_13 <- subset(weights,weights$factor == "Factor13")
weights_14 <- subset(weights,weights$factor == "Factor14")
weights_15 <- subset(weights,weights$factor == "Factor15")
weights$factor <- factor(weights$factor)
levels(weights$factor)
weights_select <- data.frame()
for (i in 1:15){
  weights_one <- subset(weights,weights$factor == levels(weights$factor)[i])
  weights_one <- weights_one[order(weights_one$value),]
  weights_one <- weights_one[c(1:10,(length(weights_one$value)-9):(length(weights_one$value))),]
  weights_one <- weights_one[order(abs(weights_one$value)),]
  weights_one$x <- seq(-2,1.8,0.2)
  weights_select <- rbind(weights_select,weights_one)
}
