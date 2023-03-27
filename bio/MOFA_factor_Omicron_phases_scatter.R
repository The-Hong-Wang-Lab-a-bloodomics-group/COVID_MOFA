library(tidyr)
library(vioplot)
library(ggsci)
library(ggpubr)
library(data.table)
library(reshape2)
library(cowplot)

# data input --------------------------------------------------------------
## The original function of MOFA -------------------------------------------
# plot_factor(MOFAobject.trained, 
#             factors = 1:9, 
#             color_by = "Omicron phases",
#             dodge = TRUE,
#             add_violin = TRUE,
#             show_missing =F
# )
# object <- MOFAobject.trained
# factors = 1:9
# color_by = "Omicron phases"
# dodge = TRUE
# add_violin = TRUE
# show_missing =F
# groups <- "group"
# a <- get_weights(MOFAobject.trained)


z <- MOFAobject.trained@expectations[["Z"]][["group1"]]
z <- as.data.frame(z)
z$ID <- rownames(z)
merge_z <- merge(list_metadata,z,by = 'ID')

# plot --------------------------------------------------------------------
data_ggplot <- merge_z

data_ggplot <- subset(data_ggplot,data_ggplot$type == 'COVID')

data_ggplot$`Omicron phases` <- factor(data_ggplot$`Omicron phases`,
                                       levels = c("Non-omicron",
                                                  "Acute omicron",
                                                  "Post-acute (Inpatient)"))
colnames(data_ggplot)
data_ggplot <- subset(data_ggplot,
                      select = -c(`...1`,
                                  `����`,
                                  bulk_ID,
                                  protein_ID,
                                  Gender,Age,
                                  `ȷ������`,
                                  `Clinical spectrum`))
data_ggplot <- subset(data_ggplot,
                      select = -c(`PLT num`,
                                  `positive CT <40`,
                                  `Ct value`,
                                  `re positive CT <40`,
                                  `re positive CT <35`,
                                  `Vaccine type`,
                                  `������������`,
                                  `Underling medical conditions`,
                                  sample,
                                  type))
data_ggplot$x <- as.numeric(factor(data_ggplot$`Omicron phases`,
                                   levels = c("Non-omicron",
                                              "Acute omicron",
                                              "Post-acute (Inpatient)")))
# comparisons <- list(c("Yes","No"))
comparisons <- list(c("Acute omicron","Post-acute (Inpatient)"),
                    c("Acute omicron","Non-omicron"),
                    c("Post-acute (Inpatient)","Non-omicron"))
# comparisons <- list(c("Moderate","Mild"),
#                     c("Moderate","Non-omicron"),
#                     c("Mild","Non-omicron"))
value_colour <- c("Non-omicron"="#00A087FF",
                  "Acute omicron"="#E64B35FF",
                  "Post-acute (Inpatient)"="#F39B7FFF")

data_ggplot_melt <- melt(data_ggplot,id.vars = c("ID","Omicron phases"))
data_ggplot_melt$x <- as.numeric(factor(data_ggplot_melt$`Omicron phases`,
                                        levels = c("Non-omicron",
                                                   "Acute omicron",
                                                   "Post-acute (Inpatient)"),
                                        labels = c(1,2,3)))
plot_list = list()


 
for( i in 1:15){
  data_ggplot_melt_one <- subset(data_ggplot_melt,
                                 data_ggplot_melt$variable == levels(data_ggplot_melt$variable)[i])
  a <- ggplot(data_ggplot_melt_one,aes(x =`Omicron phases` ,y = value))+
    geom_violin(aes(fill = `Omicron phases`))+
    geom_boxplot(outlier.shape = NA,fill = 'white',width = 0.3)+
    geom_jitter(mapping = aes(fill=`Omicron phases`),
                width = 0.4,
                shape = 21,
                size=0.5,stroke = 0.5)+
    geom_smooth(mapping = aes(x = x,y = value,colour = 'Purple'),
                alpha = 0.85,
                size = 1,
                show.legend = FALSE)+
    scale_color_manual(values = c("Purple"="#7A51A1FF"))+
    scale_y_continuous(limits = c(min(data_ggplot_melt_one$value,na.rm = T),
                                  max(data_ggplot_melt_one$value,na.rm = T)+
                                    0.38*(max(data_ggplot_melt_one$value,na.rm = T)-
                                            min(data_ggplot_melt_one$value,na.rm = T))))+
    stat_compare_means(comparisons = comparisons,
                       method = 'wilcox.test',
                       label = "p.format",size =3.5)+# format
    scale_fill_manual(values=value_colour)+
    
    theme_classic()+
    theme(axis.text.x = element_blank(),
          #axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Calibri",size=16),
          axis.text.y=element_text(size=10,face="plain",colour="black"),
          axis.title.y=element_text(family="Arial",size = 10,face="plain"),
          axis.title.x=element_text(family="Arial",size = 10,face="plain"), 
          axis.title=element_text(family="Arial",size = 10,face="plain"), 
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5),
          #plot.margin=unit(c(0.5,1.5,0.5,0.5),'lines'),
          axis.line = element_line(colour = "black",size=0.4)
          )+
    labs(title=paste0("Factor",i), x="", y="")
  plot_list[[i]] <-  a
}

prow <- plot_grid(plot_list[[1]]+theme(legend.position="none"),
                  plot_list[[13]]+theme(legend.position="none"),
                  plot_list[[8]]+theme(legend.position="none"),
                  plot_list[[9]]+theme(legend.position="none"),
                  plot_list[[3]]+theme(legend.position="none"),
                  plot_list[[5]]+theme(legend.position="none"),
                  plot_list[[11]]+theme(legend.position="none"),
                  plot_list[[12]]+theme(legend.position="none"),
                  plot_list[[15]]+theme(legend.position="none"),
                  plot_list[[4]]+theme(legend.position="none"),
                  plot_list[[2]]+theme(legend.position="none"),
                  plot_list[[6]]+theme(legend.position="none"),
                  ncol = 4,
                  nrow = 3,
                  align = "v")
prow
legend <- get_legend(plot_list[[2]] + theme(legend.position="bottom"))# bottom
p <- plot_grid( prow, legend, nrow = 2, rel_heights  = c(1, .1),rel_widths = c(1, .1))
p
