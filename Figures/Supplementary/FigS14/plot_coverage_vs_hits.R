library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggrastr)
source("../../color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))

df_sub_sp_cover<-read.delim("Targeting_level/uvig_hq_cover-vs-n_spacer_subsampled.tsv",stringsAsFactors = T)
summary(df_sub_sp_cover)
## Make categories to get percentage
tmp <- df_sub_sp_cover %>%
  mutate(cat_sp = cut(n_unique_spacer,c(0,10,max(n_unique_spacer)+1), right=F)) %>%
  mutate(cat_cover = cut(covered,c(0,200,max(covered)+1), right=F)) %>%
  group_by(cat_sp, cat_cover) %>%
  summarise(count=n()) %>%
  ungroup() %>%
  mutate(pcent=count/sum(count)*100)
# placement of labels
tmp$x<-c(5,5,5000,5000)
tmp$y<-c(20,5000,20,5000)

ggplot() + rasterise(geom_point(data=df_sub_sp_cover,aes(x=n_unique_spacer,y=covered),alpha=0.1),dpi=600) + geom_text(data=tmp,aes(x=x,y=y,label=pcent),col="darkred") + scale_y_log10(labels=scales::comma) + ylab("number of bases covered by spacer hits in the target (bp)") + xlab("number of unique spacer hits") + scale_x_log10(labels=scales::comma) + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3), legend.position="bottom")  + geom_vline(xintercept=10,col="darkred") + geom_hline(yintercept=200,col="darkred")
ggsave("Fig_S14_cover-vs-pcent.pdf",width=5,height=5)