library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
source("color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))
#
expanded_reds<-colorRampPalette(brewer.pal(9,"Reds")) # Prep scale
## Load data
df_setsize<-read.delim("Figures/Main/Fig_2/Spacer_set_size_info.tsv",stringsAsFactors = T)
df_setsize$ecosystem <- factor(df_setsize$ecosystem,ordered=T,levels=levels_ecosystem)
df_setsize$type <- factor(df_setsize$type,ordered=T,levels=rev(levels_type))
# Split set size and max coverage into bins for plots
df_setsize$bins_maxcover <- cut(x=df_setsize$max_cover,c(seq(1,10,1),seq(20,100,10),seq(200,1000,100),1E+100),include.lowest = TRUE, right = FALSE)
df_setsize$bins_setsize <- cut(x=df_setsize$spacer_set_size,c(seq(5,50,5),seq(60,200,10),seq(300,3000,100),1E+100),include.lowest = TRUE, right = FALSE)
# Prepare percentage of singleton in set
df_setsize$p_singleton <- df_setsize$n_singleton/df_setsize$spacer_set_size

p1 <- ggplot(df_setsize) + geom_boxplot(aes(x=bins_setsize,y=p_singleton,fill=bins_setsize),outliers=FALSE) + xlab("spacet set size (binned)") + ylab("percentage of singleton spacers in set") + blank_theme + scale_y_continuous(labels=scales::percent,expand=expansion(mult=0.01)) + scale_fill_manual(values=expanded_reds(length(unique(df_setsize$bins_setsize)))) + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),panel.grid.major.y=element_line(colour="lightgrey",linewidth=0.3),legend.position="none")
p2 <- ggplot(df_setsize) + geom_boxplot(aes(x=bins_maxcover,y=p_singleton,fill=bins_maxcover),outliers=FALSE) + xlab("maximum spacer coverage in set (binned)") + ylab("percentage of singleton spacers in set") + blank_theme + scale_y_continuous(labels=scales::percent,expand=expansion(mult=0.01)) + scale_fill_manual(values=expanded_reds(length(unique(df_setsize$bins_maxcover)))) + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),panel.grid.major.y=element_line(colour="lightgrey",linewidth=0.3),legend.position="none")

gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp2$widths<-gp1$widths
pdf("Fig_S7A.pdf",width=6,height=5)
grid.arrange(gp1,gp2,ncol=1,nrow=2)
dev.off()

## Proportion of singletons
# For all (note that this is only for maximum coverage set ≥ 20x, otherwise the differences between ecosystems or sets are linked to coverage, not actual differences)
df_tmp <- df_setsize %>%
  filter(max_cover>=20)
p3 <- ggplot(df_tmp) + geom_boxplot(aes(x="",y=p_singleton),fill="darkgrey") + blank_theme + scale_y_continuous(labels=scales::percent) + ylab("Percentage of singletons in set") + xlab("All sets") 
p4 <- ggplot(df_tmp) + geom_boxplot(aes(x=ecosystem,y=p_singleton,fill=ecosystem)) + blank_theme + scale_y_continuous(labels=scales::percent) + ylab("Percentage of singletons in set") + xlab("Ecosystem")+ scale_fill_manual(values=scale_ecosystem) + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),panel.grid.major.y=element_line(colour="lightgrey",linewidth=0.3),legend.position="none")
p5 <- ggplot(df_tmp) + geom_boxplot(aes(x=type,y=p_singleton,fill=type)) + blank_theme + scale_y_continuous(labels=scales::percent) + ylab("Percentage of singletons in set") + xlab("Array type")+ scale_fill_manual(values=scale_type) + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),panel.grid.major.y=element_line(colour="lightgrey",linewidth=0.3),legend.position="none")


gp3 <- ggplot_gtable(ggplot_build(p3))
gp4 <- ggplot_gtable(ggplot_build(p4))
gp5 <- ggplot_gtable(ggplot_build(p5))
gp3$heights<-gp4$heights
gp5$heights<-gp4$heights


pdf("Fig_S7B.pdf",width=7,height=4)
grid.arrange(gp3,gp4,gp5,ncol=3,nrow=1,widths=c(1.3,4,2.6))
dev.off()


