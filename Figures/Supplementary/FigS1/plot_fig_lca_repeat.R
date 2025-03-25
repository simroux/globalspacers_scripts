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
### show LCA data of arrays
df_lentax <- read.csv("lca_per_type_and_len.tsv", stringsAsFactors = T)
df_lentax$length <- factor(df_lentax$length,ordered=T,levels=c("0-10kb","10-100kb","100-500kb","500kb+"))
df_lentax$lca_rank <- factor(df_lentax$lca_rank,ordered=T,levels=levels_taxolvl)
summary(df_lentax)
## Overview isolates
p_iso<-ggplot(df_lentax[df_lentax$type=="Isolate",]) + geom_bar(aes(x=length,y=count,fill=lca_rank),col="black",position="fill",stat="identity") + scale_y_continuous(labels=scales::percent,expand=expansion(mult=0.01)) + scale_x_discrete(expand=expansion(add=0.5)) + ylab("LCA rank") + xlab("Length of shorter sequence") + scale_fill_manual(values=scale_taxorank) + blank_theme
p_mag<-ggplot(df_lentax[df_lentax$type=="Bin",]) + geom_bar(aes(x=length,y=count,fill=lca_rank),col="black",position="fill",stat="identity") + scale_y_continuous(labels=scales::percent,expand=expansion(mult=0.01)) + scale_x_discrete(expand=expansion(add=0.5)) + ylab("LCA rank") + xlab("Length of shorter sequence")+ scale_fill_manual(values=scale_taxorank) + blank_theme
gpi <- ggplot_gtable(ggplot_build(p_iso))
gpm <- ggplot_gtable(ggplot_build(p_mag))
gpi$widths<-gpi$widths
gpm$widths<-gpm$widths
#### 
df_arraytax <- read.delim("aggregated_stats_subsampling-all.tsv", stringsAsFactors = T)
df_arraytax$rank <- factor(df_arraytax$rank,ordered=T,levels=levels_taxolvl)
plot_r1<-ggplot(df_arraytax[df_arraytax$rank=="phylum" | df_arraytax$rank=="class",]) + geom_line(aes(x=n_sample,y=avg_error,col=rank)) + geom_ribbon(aes(x=n_sample,ymin=avg_error-sd_error,ymax=avg_error+sd_error,group=rank),fill="grey",alpha=0.2) + scale_y_continuous(labels=scales::percent) + ylab("average disagreement rate") + xlab("number of genomes considered") + scale_color_manual(values=scale_taxorank) + blank_theme + scale_x_continuous(breaks=seq(1,15,1),expand=expansion(mult = 0.01))  + geom_vline(xintercept=3.5,col="red",linetype=2,linewidth=0.5,alpha=0.5)
plot_r2<-ggplot(df_arraytax[df_arraytax$rank=="order" | df_arraytax$rank=="family" | df_arraytax$rank=="genus",]) + geom_line(aes(x=n_sample,y=avg_error,col=rank)) + geom_ribbon(aes(x=n_sample,ymin=avg_error-sd_error,ymax=avg_error+sd_error,group=rank),fill="grey",alpha=0.2) + scale_y_continuous(labels=scales::percent) + ylab("average disagreement rate") + xlab("number of genomes considered") + scale_color_manual(values=scale_taxorank) + blank_theme + scale_x_continuous(breaks=seq(1,15,1),expand=expansion(mult = 0.01))  + geom_vline(xintercept=3.5,col="red",linetype=2,linewidth=0.5,alpha=0.5)
gpr1 <- ggplot_gtable(ggplot_build(plot_r1))
gpr2 <- ggplot_gtable(ggplot_build(plot_r2))
## Alternative way:
gpr1$widths<-gpr2$widths

pdf("Fig_S1_draft.pdf")
grid.arrange(arrangeGrob(gpi, gpm, nrow=1, ncol=2), gpr1, gpr2, nrow=3, ncol=1)
dev.off()
