library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
source("../../color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))
## Panel A: distribution of number of spacer per virus-array combination, only for HQ virus, and only for combinations with at least 1 spacer match with 0 or 1 mismatch
df_n_hit <- read.delim("spacer_hit_distribution.tsv",stringsAsFactors = T)
df_n_hit$n_spacers<-factor(df_n_hit$n_spacers,ordered=T,levels=c("1","2","3","4","5","6","7","8","9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90-99",">= 100"))
summary(df_n_hit)
ggplot(df_n_hit) + geom_bar(aes(x=n_spacers,y=n_obs,fill=n_spacers),stat="identity",col="black") + blank_theme + scale_y_continuous(labels=scales::comma,expand=expansion(mult=0.01)) + scale_fill_manual(values=c(rep("gray",9),rep("#de2d26",10))) + theme(legend.position="none") + ylab("number of observations") + xlab("number of spacers connecting array to virus")
ggsave("panel_A_spacernumber.pdf",width=4,height=1.5)

## Panel B: Distribution of consistent/inconsistent for viruses with known hosts based on number of spacer hits
df_knownhost<-read.delim("Summary_known_host_vs_hits.tsv",stringsAsFactors = T)
df_knownhost$n_hits<-factor(df_knownhost$n_hits,ordered=T,levels=unique(df_knownhost$n_hits))
df_knownhost$match<-factor(df_knownhost$match,ordered=T,levels=c("other_taxon","unclassified_genus_but_consistent_taxonomy","known_host_genus"))
summary(df_knownhost)
tmp <- df_knownhost %>%
  group_by(n_hits,match) %>%
  summarise(total=sum(counts))
ggplot(tmp[tmp$match!="unclassified_genus_but_consistent_taxonomy",]) + geom_bar(aes(x=n_hits,y=total,fill=match),stat="identity",position="fill",color="black") + xlab("Number of spacer hits for phage-repeat pair") + ylab("Percentage of observation (ignoring unclassified taxa, unless taxonomies\nare already inconsistent for ranks where both are available)") + scale_y_continuous(labels=scales::percent) + coord_flip() + blank_theme + scale_fill_manual(values=c("#feb24c","#31a354")) + theme(legend.position="top")
ggsave("panel_B_known_hosts.pdf",width=4,height=2)
# Get numbers for paper
tmp %>%
  filter(match!="unclassified_genus_but_consistent_taxonomy") %>%
  group_by(n_hits, match) %>%
  summarise(total=sum(total)) %>%
  mutate(pcent=total/sum(total)*100)


## Panel C: Frequency of DGR across inconsistent virus-repeat pairs
df_known_raw <- read.delim("../../../Analyses/Target_IMGVR_IMGPR/Known_hosts/Consistency_known_host_with_hits.tsv",stringsAsFactors = T)
summary(df_known_raw)
length(levels(df_known_raw$uvig))
tmp <- df_known_raw %>%
  filter(n_hits>=10) %>%
  filter(quality=="Reference" | quality=="High-quality")
ggplot(tmp[tmp$result=="inconsistent",]) + geom_bar(aes(x=result,fill=dgr),position="fill",col="black") + scale_y_continuous(labels=scales::percent) + coord_flip() + scale_fill_manual(values=c("#515151","#d45c48")) + theme(legend.position="bottom") + blank_theme
ggsave("panel_C_dgr.pdf",width=4,height=0.82)
## For counts
tmp %>%
  group_by(result,dgr) %>%
  summarise(total=n()) %>%
  mutate(pcent=total/sum(total)*100)
## 57.6% DGR in the inconsistent category
tmp %>%
  group_by(dgr) %>%
  summarise(total=n()) %>%
  mutate(pcent=total/sum(total)*100)
## 16.8% DGR overall

## Panel D: Distribution of near-exact vs distant matches
theme_d<-theme_classic() + theme(text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"))
df_profile_all <- read.delim("../../../Analyses/Target_IMGVR_IMGPR/Beyond_near_exact/Virus_to_repeat_hits_profile-medium_only_with_cat_with_virusinfo.tsv",stringsAsFactors = T)
df_profile_all$category<-factor(df_profile_all$category,ordered=T,levels=c("mostly_inexact","mixed","mostly_exact"))
scale_category<-c("mostly_inexact"="#f1a340","mixed"="#f7f7f7","mostly_exact"="#998ec3")
tmp_best <- df_profile_all %>%
  arrange(desc(total)) %>%
  group_by(uvig) %>%
  summarise(category=first(category),dgr=first(dgr),known_host=first(known_host),ratio=first(ratio))
p1<-ggplot(tmp_best) + geom_histogram(aes(x=ratio,fill=category),binwidth=0.025,center=0.0125,col="black") + scale_fill_manual(values=scale_category) + theme_d + theme(legend.position="top") + ylab("number of observations")
p2<-ggplot(df_profile_all[!is.na(df_profile_all$known_host) & (df_profile_all$known_host=="consistent" | df_profile_all$known_host=="inconsistent"),]) + geom_histogram(aes(x=ratio,fill=category),binwidth=0.025,center=0.0125,col="black") + facet_grid(rows=vars(known_host), scales="free_y") + scale_fill_manual(values=scale_category) + theme_d + theme(legend.position="bottom") + ylab("number of observations")
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp1$widths<-gp2$widths
gp1$heights<-gp2$heights
pdf("panel_D_distribution.pdf",width=6,height=4)
grid.arrange(gp1,gp2,nrow=2,ncol=1)
dev.off()

## Panel E: how many distinct repeats with 10 spacer hits or more, and category distribution, for DGR+ vs DGR-
# first, the number of repeats
tmp1<-df_profile_all %>%
  group_by(dgr,uvig) %>%
  summarise(n_repeat=n()) 
p1<-ggplot(tmp1) + geom_boxplot(aes(x=dgr,y=n_repeat,fill=dgr),outliers=F) + ylab("Number of distinct repeats") + xlab("dgr")  + coord_flip() + scale_fill_manual(values=c("#515151","#d45c48")) + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3), legend.position="bottom")
# Compare median
tmp1 %>%
  group_by(dgr) %>%
  summarise(med_n=median(n_repeat))
# second, the category of these hits
p2 <- ggplot(df_profile_all) + geom_bar(aes(x=dgr,fill=forcats::fct_rev(category)),position="fill",col="black",width=0.7) + scale_fill_manual(values=scale_category)+ blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3), legend.position="bottom") + ylab("percentage of virus-repeat pairs") + scale_y_continuous(labels=scales::percent) + coord_flip()
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp1$heights<-gp2$heights
pdf("panel_E_DGR.pdf",width=5,height=1.8)
grid.arrange(gp1,gp2,ncol=2,widths=c(1.5,2))
dev.off()
