library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
source("../../color_scales.R")
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

### Show the content of the database in terms of repeat, spacer, and repeat types
# left
df_lcao_repeat <- read.delim("../../Main/Fig_1/lca_origin_to_repeat_count-for-ref.tsv", stringsAsFactors = T)
df_lcao_repeat$lca_origin <- factor(df_lcao_repeat$lca_origin,ordered=T,levels=levels_taxo)
p1 <- ggplot(df_lcao_repeat) + geom_bar(aes(x=lca_origin,y=total),col="black",fill="black",stat="identity") + scale_y_continuous(labels=scales::comma) + xlab("Predicted CRISPR repeat origin and taxonomic assignment confidence") + ylab("Number of predicted CRISPR repeats") + coord_flip() + blank_theme
# center
df_lcao_spacer <- read.delim("../../Main/Fig_1/lca_origin_to_spacer_count.tsv", stringsAsFactors = T)
df_lcao_spacer$lca_origin <- factor(df_lcao_spacer$lca_origin,ordered=T,levels=levels_taxo)
p2 <- ggplot(df_lcao_spacer) + geom_bar(aes(x=lca_origin,y=total_spacers),col="black",fill="black",stat="identity")  + scale_y_continuous(labels=scales::comma) + xlab("Predicted CRISPR repeat origin and taxonomic assignment confidence") + ylab("Number of CRISPR spacers identified across SRA") + coord_flip() + blank_theme
# right
df_lcao_type <- read.delim("../../Main/Fig_1/lca_origin_to_type.tsv")
df_lcao_type$type <- factor(df_lcao_type$type,ordered=T,levels=unique(df_lcao_type$type))
df_lcao_type <- df_lcao_type %>%
  separate_wider_delim(type, "-", names = c("type_simple"), too_many = "drop") 
df_lcao_type$type_simple <- factor(df_lcao_type$type_simple,ordered=T,levels=c("Unknown","VI","V","IV","III","II","I"))
df_lcao_type$lca_origin <- factor(df_lcao_type$lca_origin,ordered=T,levels=levels_taxo)
df_lcao_type_forplot <- df_lcao_type %>%
  group_by(type_simple,lca_origin) %>%
  summarise(total=sum(total))
p3 <- ggplot(df_lcao_type_forplot) + geom_bar(aes(x=lca_origin,y=total,fill=type_simple),stat="identity",col="black",position="fill") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_type) + xlab("Predicted CRISPR type") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom")

## Plot all in one panel
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
## Make sure everything is lined up
gp1$heights<-gp3$heights
gp2$heights<-gp3$heights

pdf("Fig_S1B_draft.pdf",width=9,height=2.5)
grid.arrange(gp1,gp2,gp3,nrow=1)
dev.off()