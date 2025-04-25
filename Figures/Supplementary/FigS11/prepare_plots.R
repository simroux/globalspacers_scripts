library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
source("color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))

## We design a new custom scale to better distinguish archaea from bacteria
scale_phylum_bis<-c('d__Archaea;p__Halobacteriota'='#fed976','d__Archaea;p__Methanobacteriota_B'='#fd8d3c','d__Archaea;p__Thermoproteota'='#e31a1c','d__Bacteria;p__Acidobacteriota'='#7fcdbb','d__Bacteria;p__Actinomycetota'='#1d91c0','d__Bacteria;p__Bacillota'='#74a9cf','d__Bacteria;p__Bacillota_A'='#0570b0','d__Bacteria;p__Bacillota_C'='#045a8d','d__Bacteria;p__Bacteroidota'='#006837','d__Bacteria;p__Campylobacterota'='#78c679','d__Bacteria;p__Deinococcota'='#8c6bb1','d__Bacteria;p__Desulfobacterota'='#fa9fb5','d__Bacteria;p__Pseudomonadota'='#dd3497','d__Bacteria;p__Synergistota'='#df65b0','d__Bacteria;p__Verrucomicrobiota'='#66c2a4','d__Viruses;p__Uroviricota'='#525252','other'='#d9d9d9',"none"="#ffffff")

## Load VR data
df_vrhits_a<-read.delim("fig_hits_input_vr-atypical_phages.tsv",stringsAsFactors = T)
df_vrhits_a$hit <- as.character(df_vrhits_a$hit)
df_vrhits_a[is.na(df_vrhits_a$hit),]$hit <- "no"
df_vrhits_a[df_vrhits_a$hit==1,]$hit <- "yes"
df_vrhits_a$hit <- factor(df_vrhits_a$hit,ordered=T,levels=c("no","yes"))
df_vrhits_a$quality <- factor(df_vrhits_a$quality,ordered=T,levels=c("Unsure","Genome fragment","High-quality","Reference"))
summary(df_vrhits_a)
summary(df_vrhits_a$hit)
p1 <- ggplot(df_vrhits_a) + geom_bar(aes(x=custom_tax),col="black") + xlab("Virus targets - class") + ylab("Number of sequences in IMG/VR (log10") + scale_y_log10(labels=scales::comma) + blank_theme + coord_flip() + theme(legend.position="right")
p2 <- ggplot(df_vrhits_a) + geom_bar(aes(x=custom_tax,fill=hit),col="black",position="fill") + xlab("Virus targets - class") + ylab("Percentage of sequences with ≥1 hit(s)") + scale_y_continuous(labels=scales::percent) + blank_theme + scale_fill_manual(values=c("white","darkgrey")) + coord_flip() + theme(legend.position="none")
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp2$heights<-gp1$heights

## Load taxa and eco data
df_atyp <- read.delim("Input_for_atypical_taxa_ecosystem.tsv",stringsAsFactors = T)
summary(df_atyp)
p3 <- ggplot(df_atyp[df_atyp$type=="tax" & df_atyp$cat!="none",]) + geom_bar(aes(x=tax,y=count,fill=cat),stat="identity",col="black",position="fill") + scale_y_continuous(labels=scales::percent) + coord_flip() + scale_fill_manual(values=scale_phylum_bis) + blank_theme + theme(legend.position="bottom")
p4 <- ggplot(df_atyp[df_atyp$type=="eco" & df_atyp$cat!="none",]) + geom_bar(aes(x=tax,y=count,fill=cat),stat="identity",col="black",position="fill") + scale_y_continuous(labels=scales::percent) + coord_flip() + scale_fill_manual(values=scale_ecosystem) + blank_theme + theme(legend.position="bottom")
gp3 <- ggplot_gtable(ggplot_build(p3))
gp4 <- ggplot_gtable(ggplot_build(p4))
## Put all 4 plots on the same pdf
pdf("atypicalphagesandarchaeoviruses_1.pdf",width=12,height=9)
grid.arrange(gp1,gp2,gp3,gp4,nrow=2,ncol=2)
dev.off()

## Counts for supplementary text:
# Number of families with a hit for each class:
# for ≥1 hit
df_atyp %>%
  group_by(custom_tax,host_prediction)

