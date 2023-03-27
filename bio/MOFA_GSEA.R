library(data.table)
library(MOFAdata)
library(MOFA2)
data <- list(as.matrix(view_bulkRNAseq_filter_scale),
             as.matrix(view_protein_filter_scale_Ensembl),
             as.matrix(view_Plasma_proteomics_filter_scale),
             as.matrix(view_Metabolomic_filter_scale))
names(data) <- c('Platelet transcriptome',
                 'Platelet proteome',
                 'Plasma proteome',
                 'Plasma metabolomic')
lapply(data,dim)
MOFAobject <- create_mofa(data)
# run MOFA ----------------------------------------------------------------
plot_data_overview(MOFAobject)
print(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 15
head(model_opts)
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
head(train_opts)
MOFAobject_data <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

## run MOFA --------------------------------------------------------------------

MOFAobject.trained <- run_mofa(MOFAobject_data,
                               outfile = file.path(getwd(),"model_TPM_t_Ensembl.hdf5"))

# GSEA --------------------------------------------------------------------

res.positive_protein <- run_enrichment(MOFAobject.trained, 
                                                 feature.sets = reactomeGS, 
                                                 view = "Platelet proteome",
                                                 sign = "positive"
)
res.negative_protein <- run_enrichment(MOFAobject.trained, 
                                                 feature.sets = reactomeGS, 
                                                 view = "Platelet proteome",
                                                 sign = "negative"
)
save(res.positive_protein,file = 'data_GSEA/res.positive_protein.rdata')
save(res.negative_protein,file = 'data_GSEA/res.negative_protein.rdata')
save(res.positive,file = 'data_GSEA/res.positive_RNA.rdata')
save(res.negative,file = 'data_GSEA/res.negative_RNA.rdata')
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
save(res.positive_Plasma_proteomics,file = 'data_GSEA/res.positive_Plasma_proteomics.rdata')
save(res.negative_Plasma_proteomics,file = 'data_GSEA/res.negative_Plasma_proteomics.rdata')

plot_enrichment(res.positive_protein, factor = 11, max.pathways = 20)
plot_enrichment(res.negative_protein, factor = 11, max.pathways = 20)

