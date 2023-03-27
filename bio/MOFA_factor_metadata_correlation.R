library(Hmisc)
z <- MOFAobject.trained@expectations[["Z"]][["group1"]]
z <- as.data.frame(z)
z$ID <- rownames(z)
merge_z <- merge(list_metadata,z,by = 'ID')
# The original function of MOFA -------------------------------------------
# covariates <- subset(merge_z,select = c(`Omicron phases`,
#                                         `Underling medical conditions`,
#                                         `Ct value`,
#                                         Gender,
#                                         Age,
#                                         `PLT num`,
#                                         `Clinical spectrum`,
#                                         `Vaccine type`))
# covariates <- metadata[match(rownames(covariates), metadata$sample),]
# cols <- which(sapply(covariates, is.character))
# if (length(cols >= 1)) {
#   covariates[cols] <- lapply(covariates[cols], as.factor)
# }
# cols <- which(!sapply(covariates, class) %in% c("numeric",
#                                                 "integer"))
# if (length(cols >= 1)) {
#   cols.factor <- which(sapply(covariates, class) == "factor")
#   covariates[cols] <- lapply(covariates[cols], as.numeric)
#   warning("There are non-numeric values in the covariates data.frame, converting to numeric...")
#   covariates[cols] <- lapply(covariates[cols], as.numeric)
# }
# stopifnot(all(sapply(covariates, class) %in% c("numeric",
#                                                "integer")))
x_cor <- merge_z
rownames(x_cor) <- x_cor$ID
# y_cor = factors
y_cor <- x_cor[,21:35]
# x_cor = metadata
x_cor <- subset(x_cor,select = c(`Omicron phases`,
                                 `Underling medical conditions`,
                                 `Ct value`,
                                 Gender,
                                 Age,
                                 `PLT num`,
                                 `Clinical spectrum`,
                                 `Vaccine type`))

cortest <- rcorr(x=as.matrix(x_cor),
                 y=as.matrix(y_cor),
                 type = "spearman")# spearman or pearson
View(cortest[["P"]][9:23,1:8])
# cortest <- corr.test(x_cor, y_cor, method = "pearson", 
#                      adjust = "BH")
cor_r <- t(cortest$r)
cor_p <- t(cortest$p)
cor_p[cor_p > 0.05] <- 1
if (all(cor_p == 1)) 
  stop("All p-values are 1.0, nothing to plot")
cor_p <- -log10(cor_p)
cor_p[is.infinite(cor_p)] <- 1000