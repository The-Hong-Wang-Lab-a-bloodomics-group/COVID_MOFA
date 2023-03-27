# package -----------------------------------------------------------------
library(ggplot2)
library(tidyr)
library(vioplot)
library(ggsci)
library(ggpubr)
library(data.table)
library(reshape2)
library(cowplot)
library(readxl)

# data input --------------------------------------------------------------
data_ggplot <- as.data.frame(t(view_Metabolomic_filter_scale))
data_ggplot$sample <- rownames(data_ggplot)
data_ggplot <- merge(list_metadata_type,
                     data_ggplot,
                     by = 'sample')

data_ggplot$`Omicron phases` <- factor(data_ggplot$`Omicron phases,
                                       levels = c("Non-omicron",
                                                  "Acute omicron",
                                                  "Post-acute (Inpatient)"))
colnames(data_ggplot)
Top_gene_to_draw_figures <- read_excel("data/Top gene to draw figures.xlsx")
## Plasma_metabolome ------------------------------------------------------
data_ggplot <- subset(data_ggplot,select = c(`sample`,
                                             `Omicron phases`,
                                             `Uric acid`,
                                             `DL-Tryptophan`,
                                             `DL-Carnitine`,
                                             `L-Phenylalanine`,
                                             `2-Hydroxycinnamic acid`,
                                             `L-Glutamine`,
                                             `Creatinine`))
colnames(data_ggplot) <- c("sample",
                           "Omicron phases",
                           Top_gene_to_draw_figures$`Plasma metabolome`)
# others ------------------------------------------------------------------
# change the column of Top_gene_to_draw_figures
# Plasma proteome
data_ggplot <- subset(data_ggplot,
                      select = c(`sample`,
                                 `Omicron phases`,
                                 Top_gene_to_draw_figures$`Plasma proteome`))
# Platelet proteome
data_ggplot <- subset(data_ggplot,
                      select = c(`sample`,
                                 `Omicron phases`,
                                 Top_gene_to_draw_figures$`Platelet proteome`))
# Platelet transcriptome
data_ggplot <- subset(data_ggplot,
                      select = c(`sample`,
                                 `Omicron phases`,
                                 Top_gene_to_draw_figures$`Platelet transcriptome`))

# ggplot parameters -------------------------------------------------------
comparisons <- list(c("Acute omicron","Post-acute (Inpatient)"),
                    c("Acute omicron","Non-omicron"),
                    c("Post-acute (Inpatient)","Non-omicron"))
value_colour <- c("Non-omicron"="#00A087FF",
                  "Acute omicron"="#E64B35FF",
                  "Post-acute (Inpatient)"="#F39B7FFF")
plot_list = list()
data_ggplot_melt <- reshape::melt(data_ggplot,
                                  id.vars = c("sample","Omicron phases"))
data_ggplot_melt$x <- as.numeric(factor(data_ggplot_melt$`Omicron phases`,
                                        levels = c("Non-omicron",
                                                   "Acute omicron",
                                                   "Post-acute (Inpatient)")))
plot_list = list()
# data_ggplot_melt_Plasma_metabolome <- data_ggplot_melt
# data_ggplot_melt <- data_ggplot_melt_Platelet_transcriptome
for( i in 1:8){
  data_ggplot_melt_one <- subset(data_ggplot_melt,
                                 data_ggplot_melt$variable == levels(data_ggplot_melt$variable)[i])
  plot_list[[i]] <- ggplot(data_ggplot_melt_one,aes(x =`Omicron phases` ,y = value))+
    geom_violin(aes(fill = `Omicron phases`))+
    geom_boxplot(outlier.shape = NA,fill = 'white',width = 0.2)+
    geom_jitter(mapping = aes(fill=`Omicron phases`),
                width = 0.4,
                shape = 21,
                size=0.5,stroke = 0.5)+
    geom_smooth(mapping = aes(x = x,y = value),colour = '#7A51A1',
                alpha = 0.85,
                size = 1,
                show.legend = FALSE)+
    scale_y_continuous(limits = c(min(data_ggplot_melt_one$value,na.rm = T),
                                  max(data_ggplot_melt_one$value,na.rm = T)+
                                    0.38*(max(data_ggplot_melt_one$value,na.rm = T)-
                                            min(data_ggplot_melt_one$value,na.rm = T))))+
    stat_compare_means(comparisons = comparisons,
                       method = 'wilcox.test',
                       label = "p.format",size =2.5)+# format
    scale_fill_manual(values=value_colour)+
    scale_color_manual(values = c("Purple"="#7A51A1FF"))+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          #axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Arial",size=16), 
          axis.text.y=element_text(family="Arial",size=9,face="plain",colour="black"),
          axis.title.y=element_text(family="Arial",size = 9,face="plain"),
          axis.title.x=element_text(family="Arial",size = 9,face="plain"), 
          axis.title=element_text(family="Arial",size = 9,face="plain"), 
          panel.border = element_blank(),
          #plot.title = element_text(hjust = 0.5),
          #plot.margin=unit(c(0.5,1.5,0.5,0.5),'lines'),
          axis.line = element_line(colour = "black",size=0.4)
    )+
    labs(title=paste0(levels(data_ggplot_melt$variable)[i]), x="", y="")
  plot_list[[i]] <-  a
}

prow <- plot_grid(plot_list[[1]]+theme(legend.position="none"),
                  plot_list[[2]]+theme(legend.position="none"),
                  plot_list[[3]]+theme(legend.position="none"),
                  plot_list[[4]]+theme(legend.position="none"),
                  plot_list[[5]]+theme(legend.position="none"),
                  plot_list[[6]]+theme(legend.position="none"),
                  plot_list[[7]]+theme(legend.position="none"),
                  plot_list[[8]]+theme(legend.position="none"),
                  ncol = 4,
                  nrow = 2,
                  align = "v")
prow
legend <- get_legend(plot_list[[2]] + theme(legend.position="bottom"))# bottom
p <- plot_grid( prow, 
                legend, 
                nrow = 2, 
                rel_heights = c(1, .1),
                rel_widths = c(1, .1))
p
