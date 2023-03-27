library(ggplot2)
library(readxl)

# data input --------------------------------------------------------------
#GSEA_all <- read_excel("scale_TPM/102/COVID_MOFA_GSEA_table_selectedv02/GSEA_all.xlsx")
GSEA_all_select <- read_excel("scale_TPM/102/COVID_MOFA_GSEA_table_selectedv02/GSEA_all_select.xlsx")
# GSEA Metabolomic result select has been processed
GSEA_Metabolomic_result_select <- read_excel("scale_TPM/102/Metabolomic/GSEA_Metabolomic_result_select.xlsx")

# data conversion ---------------------------------------------------------
data_ggplot <- reshape2::melt(GSEA_all_select,
                              id = c("Description","View","sign"))
data_ggplot$logp <- -log10(data_ggplot$value + 1e-100)
data_ggplot <- rbind(data_ggplot,
                     GSEA_Metabolomic_result_select)
data_ggplot$View <- factor(data_ggplot$View,
                           levels = c("Platelet transcriptome",
                                      "Platelet proteome",
                                      "Plasma proteome","Metabolomic"))
data_ggplot$sign <- factor(data_ggplot$sign,
                           levels = c("positive",
                                      "negative"))
data_ggplot$variable <- factor(data_ggplot$variable,
                               levels = c("Factor1",
                                          "Factor13",
                                          "Factor8",
                                          "Factor9",
                                          "Factor3",
                                          "Factor5",
                                          "Factor11",
                                          "Factor12",
                                          "Factor15",
                                          "Factor4",
                                          "Factor2",
                                          "Factor6",
                                          "Factor7",
                                          "Factor10",
                                          "Factor14"))
data_ggplot <- subset(data_ggplot,data_ggplot$variable != "Factor7")
data_ggplot <- subset(data_ggplot,data_ggplot$variable != "Factor10")
data_ggplot <- subset(data_ggplot,data_ggplot$variable != "Factor14")
for (i in 1:nrow(data_ggplot)) {
  if(data_ggplot$value[i] > 0.05){
    data_ggplot$logp[i] <- NA
  }
}

# data visualization ------------------------------------------------------
value_colour <- c("positive" = "#ED1E1AFF",
                  "negative" = "#3F66E1FF")
# basic setting
ggplot(data_ggplot,aes(x = variable, y = Description))+
  # point figure
  geom_point(aes(color=sign, size=logp))+
  # theme setting
  theme_classic()+
  # color settings
  scale_color_manual(values = value_colour)+
  # details setting
  theme(#axis.text.x = element_blank(),
        axis.text.x=element_text(angle=45,
                                 hjust = 1,
                                 colour="black",
                                 family="Arial",
                                 size=8), 
        axis.text.y=element_text(size=10,
                                 family = "Arial",
                                 face="plain",
                                 colour="black"),
        axis.title.y=element_text(family="Arial",
                                  size = 8,
                                  face="plain"),
        axis.title.x=element_text(family="Arial",
                                  size = 8,
                                  face="plain"), 
        axis.title=element_text(family="Arial",
                                size = 8,
                                face="plain"), 
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour = "black",
                                 size=0.4),
        legend.position = "none",
        plot.margin=unit(c(0.5,1.5,0.5,0.5),
                         'lines'))+
  # Faceted setting
    facet_grid(View~., scale='free',space = 'free')+
  labs(title="", x="", y="")
