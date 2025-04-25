library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
source("color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))
## Panel A: distribution of number of spacer per virus-array combination, only for HQ virus, and only for combinations with at least 1 spacer match with 0 or 1 mismatch
df_n_hit <- read.delim("spacer_hit_distribution.tsv",stringsAsFactors = T)
df_n_hit$n_spacers<-factor(df_n_hit$n_spacers,ordered=T,levels=c("1","2","3","4","5","6","7","8","9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90-99",">= 100"))
summary(df_n_hit)
ggplot(df_n_hit) + geom_bar(aes(x=n_spacers,y=n_obs,fill=n_spacers),stat="identity",col="black") + blank_theme + scale_y_continuous(labels=scales::comma,expand=expansion(mult=0.01)) + scale_fill_manual(values=c(rep("gray",9),rep("#de2d26",10))) + theme(legend.position="none") + ylab("number of observations") + xlab("number of spacers connecting array to virus")
ggsave("panel_A_spacernumber.pdf",width=4,height=1.5)


## Panel B: overall coverage vs per-sample coverage for viruses targeted on at least 10% of their length or 1kb
library(ggrastr)
library(grid)
df_hit_samples <- read.delim("spacer_hit_distribution.tsv",stringsAsFactors = T)
## Make categories - in the end we won't use high / very-high, but we will color based on whether the max sample represents more or less of 20% of the total coverage
df_hit_samples$category <- case_when( (df_hit_samples$local_cover/df_hit_samples$global_cover) >= 0.2 ~ "local", (df_hit_samples$local_cover/df_hit_samples$global_cover) < 0.2 ~ "global_only")
df_hit_samples$category <- factor(df_hit_samples$category)
summary(df_hit_samples)
p2 <- ggplot(df_hit_samples) + geom_histogram(aes(x=global_cover,fill=category),col="black",breaks=seq(0,100,2)) + xlim(-0.05,100) + scale_y_continuous(labels=scales::comma) + xlab("") + ylab("Number of observations") + blank_theme + scale_fill_manual(values=c("#66c2a5","#fc8d62")) + theme(legend.position="none")
p3 <- ggplot(df_hit_samples) + geom_histogram(aes(x=local_cover/global_cover,fill=category),col="black",breaks=seq(0,1,0.02)) + scale_y_continuous(labels=scales::comma) + xlab("") + ylab("Number of observations") + blank_theme + coord_flip(xlim=c(0,1))+ scale_fill_manual(values=c("#66c2a5","#fc8d62"))+ theme(legend.position="none")
p1 <- ggplot(df_hit_samples) + rasterise(geom_point(aes(x=global_cover,y=local_cover/global_cover,col=category),alpha=0.1),dpi=600) + scale_y_continuous(labels=scales::percent,lim=c(0,1)) + scale_color_manual(values=c("#66c2a5","#fc8d62")) + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3)) + xlim(0,100) + xlab("virus coverage by spacers (%) across all samples") + ylab("maximum virus coverage in individual sample, relative to total coverage across samples (x-axis)") + theme(legend.position="none")
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp2$widths<-gp1$widths
gp3$heights<-gp1$heights
blank <- grid.rect(gp=gpar(col="white"))
pdf("panel_B_coverage_by_sample.pdf",width=5,height=5)
grid.arrange(gp2,blank,gp1,gp3,nrow=2,ncol=2,widths=c(4,1),heights=c(1,4))
dev.off()

## Some numbers for papers:
# % of global only vs local
df_hit_samples %>%
  group_by(category) %>%
  summarise(count=n()) %>%
  mutate(pcent=count/sum(count)*100)
# Number of high targeting
nrow(df_hit_samples)

## Panels C and D - mismatch profiles
df_profiles<-read.delim("../../../Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/Virus_to_repeat_hits_profile-medium_only_with_cc_with_virusinfo.tsv",stringsAsFactors = T)
df_profiles$bins <- cut(df_profiles$cc,breaks = seq(-1.025,1.025,0.05))
df_profiles$type <- factor(df_profiles$type,ordered=T,levels=c("medium","high"))
colourCount <- length(unique(df_profiles$bins))
summary(df_profiles)

## Get numbers for paper
# Overall proportion
df_profiles %>%
  group_by(category) %>%
  summarise(total=n()) %>%
  mutate(pcent=total/sum(total)*100)
# Proportion by type
df_profiles %>%
  group_by(type,category) %>%
  summarise(total=n()) %>%
  mutate(pcent=total/sum(total)*100)
## First part of figure - the x-y plot to illustrate how the correlation looks
df_profiles_sampled <- df_profiles %>% 
  sample_n(100) %>%
  pivot_longer(cols=n_hit_0:n_hit_3, names_to="n_mismatch", values_to="number") %>%
  unite(code,c("uvig","array")) %>%
  select(c("code","cc","type","number","n_mismatch")) %>%
  group_by(code) %>%
  mutate(z_score = (number-mean(number))/sd(number)) %>%
  ungroup()

p1<-ggplot(df_profiles_sampled) + geom_line(aes(x=n_mismatch,y=z_score,group=code,col=cc),alpha=0.5) + geom_point(aes(x=n_mismatch,y=z_score,group=code,col=cc),alpha=0.5,size=2) + scale_y_continuous(labels=scales::comma,expand=expansion(mult=0.02)) + ylab("Normalized number of spacers (Z-score transform)") + scale_colour_distiller(palette="RdGy",direction=-1,name="Correlation coefficient",lim=c(-1,1)) + xlab("Number of mismatches between spacer and virus") + scale_x_discrete(labels=c("0","1","2","3"),expand=expansion(add=0.1)) + blank_theme + theme(legend.position="bottom")
## next is the overall distribution across all combinations included here
## Note: we use brewer.pal(7) because the help from scale_colour_distiller says " smoothly interpolating 7 colours from any palette" and we try to keep the same scale
p2<-ggplot(df_profiles) + geom_histogram(aes(x=cc,fill=bins),binwidth = 0.05,col="black",alpha=0.8) + scale_fill_manual(values=colorRampPalette(rev(brewer.pal(7, "RdGy")))(colourCount)) + blank_theme  + theme(legend.position="null") + xlab("Correlation coefficient") + ylab("Number of virus-array pairs (≥ 10 spacer hits) XX AND AT LEAST ONE HIGH TARGETING FOR THE VIRUS") + scale_x_continuous(expand = expansion(add=0.03)) + scale_y_continuous(labels=scales::comma,expand=expansion(mult=0.01)) + geom_vline(xintercept=-0.7,linetype="dashed")  + geom_vline(xintercept=0.7,linetype="dashed")
p3 <- ggplot(df_profiles[df_profiles$type=="medium",]) + geom_histogram(aes(x=cc,fill=bins),binwidth = 0.05,col="black",alpha=0.8) + scale_fill_manual(values=colorRampPalette(rev(brewer.pal(7, "RdGy")))(colourCount)) + blank_theme  + theme(legend.position="null") + xlab("Correlation coefficient") + ylab("Number of virus-array pairs (≥ 10 spacer hits) XX AND AT LEAST ONE HIGH TARGETING FOR THE VIRUS") + scale_x_continuous(expand = expansion(add=0.03)) + scale_y_continuous(labels=scales::comma,expand=expansion(mult=0.01)) + geom_vline(xintercept=-0.7,linetype="dashed")  + geom_vline(xintercept=0.7,linetype="dashed") 
p4 <- ggplot(df_profiles[df_profiles$type=="high",]) + geom_histogram(aes(x=cc,fill=bins),binwidth = 0.05,col="black",alpha=0.8) + scale_fill_manual(values=colorRampPalette(rev(brewer.pal(7, "RdGy")))(colourCount)) + blank_theme  + theme(legend.position="null") + xlab("Correlation coefficient") + ylab("Number of virus-array pairs (≥ 10 spacer hits) XX AND AT LEAST ONE HIGH TARGETING FOR THE VIRUS") + scale_x_continuous(expand = expansion(add=0.03)) + scale_y_continuous(labels=scales::comma,expand=expansion(mult=0.01)) + geom_vline(xintercept=-0.7,linetype="dashed")  + geom_vline(xintercept=0.7,linetype="dashed") 
## Finally, add the plasmid for comparison
df_profiles_plasmid<-read.delim("../../../Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/Plasmid_to_repeat_hits_profile-medium_with_cc.tsv",stringsAsFactors = T)
### Distribution of correlation coefficient
df_profiles_plasmid$bins <- cut(df_profiles_plasmid$cc,breaks = seq(-1.025,1.025,0.05))
colourCountp <- length(unique(df_profiles_plasmid$bins))
p5<-ggplot(df_profiles_plasmid) + geom_histogram(aes(x=cc,fill=bins),binwidth = 0.05,col="black",alpha=0.8) + scale_fill_manual(values=colorRampPalette(rev(brewer.pal(7, "RdGy")))(colourCountp)) + blank_theme  + theme(legend.position="null") + xlab("Correlation coefficient") + ylab("Number of virus-array pairs (≥ 10 spacer hits)") + scale_x_continuous(expand = expansion(add=0.03)) + scale_y_continuous(expand=expansion(mult=0.01),labels=scales::comma) + geom_vline(xintercept=-0.7,linetype="dashed")  + geom_vline(xintercept=0.7,linetype="dashed")
## Get numbers for paper
# Overall proportion
df_profiles_plasmid %>%
  group_by(category) %>%
  summarise(total=n()) %>%
  mutate(pcent=total/sum(total)*100)


gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp4 <- ggplot_gtable(ggplot_build(p4))
gp5 <- ggplot_gtable(ggplot_build(p5))
gp2$widths<-gp1$widths
gp3$widths<-gp1$widths
gp4$widths<-gp1$widths
gp5$widths<-gp1$widths
pdf("panels_hits_vs_mismatches.pdf",width=3,height=6)
grid.arrange(gp1,gp2,gp3,gp4,gp5,nrow=5,heights=c(3,1,1,1,1))
dev.off()


### Panel E: get distribution of targeting level and type
df_fc_type <- read.delim("Final_counts_for_type_vs_level-and-dgr.tsv",stringsAsFactors = T)
df_fc_type$type<-factor(df_fc_type$type,ordered=T,levels=rev(c("low","medium","high","negative","intermediary","positive")))
color_scale_type <- c("low"="#ffffd9","medium"="#7fcdbb","high"="#225ea8","negative"="#515151","positive"="#d45c48","intermediary"="white")
p1 <- ggplot(df_fc_type[df_fc_type$metric=="type" & df_fc_type$category=="all",]) + geom_bar(aes(x=value,y=count,fill=type),stat="identity",col="black",position="fill") + blank_theme + coord_flip() + xlab("") + scale_y_continuous(labels=scales::percent) + ylab("Number of observations") + scale_fill_manual(values=color_scale_type)
p2 <- ggplot(df_fc_type[df_fc_type$metric=="race" & df_fc_type$category=="all" & !is.na(df_fc_type$type),]) + geom_bar(aes(x=value,y=count,fill=type),stat="identity",col="black",position="fill") + blank_theme + coord_flip() + xlab("") + scale_y_continuous(labels=scales::percent) + ylab("Number of observations") + scale_fill_manual(values=color_scale_type)
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp2$widths<-gp1$widths
pdf("panel_level_and_type.pdf",width=5,height=1.5)
grid.arrange(gp1,gp2)
dev.off()

### Panel F - Prediction of targeting type
df_knownhost<-read.delim("Summary_known_host_vs_hits.tsv",stringsAsFactors = T)
df_knownhost$n_hits<-factor(df_knownhost$n_hits,ordered=T,levels=unique(df_knownhost$n_hits))
df_knownhost$profile<-factor(df_knownhost$profile,ordered=T,levels=c("negative","unknown","positive"))
df_knownhost$match<-factor(df_knownhost$match,ordered=T,levels=c("other_taxon","unclassified_genus_but_consistent_taxonomy","known_host_genus"))
df_knownhost_f <- df_knownhost %>%
  filter(n_hits!="1-4" & n_hits!="5-9") %>%
  filter(match!="unclassified_genus_but_consistent_taxonomy") %>%
  group_by(profile,match) %>%
  summarise(total=sum(counts))
df_knownhost_f$profile<-factor(df_knownhost_f$profile,ordered=T,levels=c("positive","unknown","negative"))
ggplot(df_knownhost_f) + geom_bar(aes(x=profile,y=total,fill=match),stat="identity",position="fill",color="black") + xlab("Number of spacer hits for phage-repeat pair") + ylab("Percentage of observation (ignoring unclassified taxa, unless taxonomies\nare already inconsistent for ranks where both are available)") + scale_y_continuous(labels=scales::percent) + coord_flip() + blank_theme + scale_fill_manual(values=c("#feb24c","#31a354")) + theme(legend.position="top")
ggsave("known_host_panel.pdf",width=5,height=2)

#### Panel G - Link between mismatch profile and DGR
ggplot(df_fc_type[df_fc_type$metric=="race" & df_fc_type$category=="dgr" & !is.na(df_fc_type$type),]) + geom_bar(aes(x=value,y=count,fill=type),stat="identity",col="black",position="fill") + blank_theme + coord_flip() + xlab("") + scale_y_continuous(labels=scales::percent) + ylab("Number of observations") + scale_fill_manual(values=color_scale_type) + theme(legend.position="none")
ggsave("dgr_bars.pdf",width=5,height=1.3)