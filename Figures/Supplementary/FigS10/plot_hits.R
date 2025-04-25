## Load fonts and colors for plots
library(sysfonts)
library(extrafont)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
source("color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"), legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3 Medium"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3 Medium"),legend.key.size = unit(1,"line"))

## Distribution of VR hits across ecosystems
df_vrhits<-read.delim("../../Main/Fig_3/fig_hits_input_vr.tsv",stringsAsFactors = T)
df_vrhits$hit <- as.character(df_vrhits$hit)
df_vrhits[is.na(df_vrhits$hit),]$hit <- "no"
df_vrhits[df_vrhits$hit==1,]$hit <- "yes"
df_vrhits$hit <- factor(df_vrhits$hit,ordered=T,levels=c("no","yes"))
df_vrhits$ecosystem<-factor(df_vrhits$ecosystem,ordered=T,levels=levels_ecosystem)
summary(df_vrhits$hit)
p1 <- ggplot(df_vrhits) + geom_bar(aes(x=ecosystem,fill=hit),col="black",position="fill") + xlab("Virus targets - ecosystem") + ylab("Percentage of sequences with ≥1 hit(s)") + scale_y_continuous(labels=scales::percent) + blank_theme + scale_fill_manual(values=c("white","darkgrey")) + coord_flip() + theme(legend.position="none")
tmp_vrhits <- df_vrhits %>%
  filter(quality=='High-quality' | quality=='Reference')
p2 <- ggplot(tmp_vrhits) + geom_bar(aes(x=ecosystem,fill=hit),col="black",position="fill") + xlab("Virus targets - ecosystem (HQ only)") + ylab("Percentage of sequences with ≥1 hit(s)") + scale_y_continuous(labels=scales::percent) + blank_theme + scale_fill_manual(values=c("white","darkgrey")) + coord_flip() + theme(legend.position="none")
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp2$heights<-gp1$heights
pdf("Fig_SXX_hits_vr-ecosystem-HQ.pdf",width=9,height=4)
grid.arrange(gp1,gp2,ncol=2)
dev.off()

## Distribution of PR hits across ecosystems
df_vrhits<-read.delim("../../Main/Fig_3/fig_hits_input_pr.tsv",stringsAsFactors = T)
df_vrhits$hit <- as.character(df_vrhits$hit)
df_vrhits[is.na(df_vrhits$hit),]$hit <- "no"
df_vrhits[df_vrhits$hit==1,]$hit <- "yes"
df_vrhits$hit <- factor(df_vrhits$hit,ordered=T,levels=c("no","yes"))
df_vrhits$ecosystem<-factor(df_vrhits$ecosystem,ordered=T,levels=levels_ecosystem)
summary(df_vrhits$hit)
p1 <- ggplot(df_vrhits) + geom_bar(aes(x=ecosystem,fill=hit),col="black",position="fill") + xlab("Virus targets - ecosystem") + ylab("Percentage of sequences with ≥1 hit(s)") + scale_y_continuous(labels=scales::percent) + blank_theme + scale_fill_manual(values=c("white","darkgrey")) + coord_flip() + theme(legend.position="none")
tmp_vrhits <- df_vrhits %>%
  filter(quality=='High-quality' | quality=='Reference')
p2 <- ggplot(tmp_vrhits) + geom_bar(aes(x=ecosystem,fill=hit),col="black",position="fill") + xlab("Virus targets - ecosystem (HQ only)") + ylab("Percentage of sequences with ≥1 hit(s)") + scale_y_continuous(labels=scales::percent) + blank_theme + scale_fill_manual(values=c("white","darkgrey")) + coord_flip() + theme(legend.position="none")
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp2$heights<-gp1$heights
pdf("Fig_SXX_hits_pr-ecosystem-HQ.pdf",width=9,height=4)
grid.arrange(gp1,gp2,ncol=2)
dev.off()

### Ecosystem vs Ecosystem hit frequency heatmap for VR
df_vr_hits_eco <- read.delim("vr_ecosystem_heatmap.tsv", stringsAsFactors = T)
summary(df_vr_hits_eco)
df_vr_hits_eco$virus_ecosystem <- factor(df_vr_hits_eco$virus_ecosystem,ordered=T,levels=levels_ecosystem)
df_vr_hits_eco$spacer_ecosystem <- factor(df_vr_hits_eco$spacer_ecosystem,ordered=T,levels=levels_ecosystem)
ggplot(df_vr_hits_eco) + geom_tile(aes(x=virus_ecosystem,y=spacer_ecosystem,fill=p_observation)) + scale_fill_distiller(palette="YlGn",direction=1,na.value="white",lim=c(1,100)) + blank_theme + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
ggsave("heatmap_virus.pdf",height=5,width=5)

## Same heatmap but for hits to PR
df_pr_hits_eco <- read.delim("pr_ecosystem_heatmap.tsv", stringsAsFactors = T)
summary(df_pr_hits_eco)
df_pr_hits_eco <- df_pr_hits_eco %>%
  filter(plasmid_ecosystem!="Other") %>%
  filter(spacer_ecosystem!="Other")
df_pr_hits_eco$plasmid_ecosystem <- factor(df_pr_hits_eco$plasmid_ecosystem,ordered=T,levels=levels_ecosystem)
df_pr_hits_eco$spacer_ecosystem <- factor(df_pr_hits_eco$spacer_ecosystem,ordered=T,levels=levels_ecosystem)
ggplot(df_pr_hits_eco) + geom_tile(aes(x=plasmid_ecosystem,y=spacer_ecosystem,fill=p_observation)) + scale_fill_distiller(palette="YlGn",direction=1,na.value="white",lim=c(1,100)) + blank_theme + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
ggsave("heatmap_plasmid.pdf",height=5,width=5)
