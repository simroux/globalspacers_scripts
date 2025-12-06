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
## Load data
df_sample_detect <-read.delim("../../Main/Fig_1/sample_to_array_counts.tsv",stringsAsFactors = T)
df_sample_detect$sra_run <- df_sample_detect$sra_run_1
df_sample_detect <- df_sample_detect %>%  
  mutate(total_cluster = replace_na(total_cluster, 0))
df_sample_detect$ecosystem_sum <- factor(df_sample_detect$ecosystem_sum,ordered=T,levels=levels_ecosystem)

## Check distribution of n_base to decide of a cutoff when calculating median number of array per 1Gb
ggplot(df_sample_detect) + geom_histogram(aes(x=n_base),bins=200) + scale_x_log10() + geom_vline(xintercept=5e+08)
## We'll use 2.5e+08 for minimum dataset size to be considered in the median distribution

# Left panel is count
df_lp <- df_sample_detect %>%
  group_by(ecosystem_base,ecosystem_sum) %>%
  summarise(total=n()) %>%
  mutate(ecosystem_base=as.character(ecosystem_base)) %>%
  arrange(ecosystem_sum,ecosystem_base) # %>%
#
df_lp$ecosystem_base <- factor(df_lp$ecosystem_base,ordered=T,levels=unique(df_lp$ecosystem_base))
p1 <- ggplot(df_lp,aes(x=ecosystem_base,y=total)) + geom_bar(fill="grey",stat="identity",alpha=0.8,width=0.75,col="black") + scale_y_continuous(labels=scales::comma) + xlab("Ecosystem type") + ylab("Total number of samples") + coord_flip() + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3), legend.position="none")

# Right panel is median detection per 1Gb, with subsampling to get estimation of sample-to-sample variability
df_sta_med_ecod <- read.delim("../../Main/Fig_1/sample_to_array_counts-summarized_by_ecosystem-detailed.tsv",stringsAsFactors = T)
df_sta_med_ecod$ecosystem_sum <- factor(df_sta_med_ecod$ecosystem_sum,ordered=T,levels=levels_ecosystem)
summary(df_sta_med_ecod)
df_rp <- df_sta_med_ecod %>%
  mutate(ecosystem_detailed=as.character(ecosystem_detailed)) %>%
  arrange(ecosystem_sum,ecosystem_detailed) %>%
  mutate(ecosystem_detailed = factor(ecosystem_detailed,ordered=T,levels=df_lp$ecosystem_base))
# 
p2 <- ggplot(df_rp,aes(x=ecosystem_detailed,y=average_median)) + geom_bar(aes(fill=ecosystem_sum),stat="identity",alpha=1,width=0.75,col="black") + geom_linerange(aes(ymin=min_median, ymax=max_median), col="gray30", alpha=0.8) + scale_fill_manual(values=scale_ecosystem) + scale_x_discrete(limits=levels(df_rp$ecosystem_detailed)) + xlab("Ecosystem type") + ylab("Median number of CRISPR array per 1Gb") + coord_flip() + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3), legend.position="none")
## Make this top panel
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
## One way
gp1$widths <- gp2$widths
gp1$heights <- gp2$heights
##
grid.arrange(gp1,gp2,ncol=2)
pdf("Fig_S4_bottom_panel.pdf",width=12,height=6)
grid.arrange(gp1,gp2,ncol=2)
dev.off()

