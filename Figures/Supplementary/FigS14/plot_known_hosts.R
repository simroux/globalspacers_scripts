library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
source("../../color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))

df_knownhost<-read.delim("../../Main/Fig_4/Summary_known_host_vs_hits.tsv",stringsAsFactors = T)
df_knownhost$ref_host_source<-factor(df_knownhost$ref_host_source,ordered=T,levels=unique(df_knownhost$ref_host_source))
df_knownhost$n_hits<-factor(df_knownhost$n_hits,ordered=T,levels=unique(df_knownhost$n_hits))
df_knownhost$profile<-factor(df_knownhost$profile,ordered=T,levels=c("negative","unknown","positive"))
df_knownhost$match<-factor(df_knownhost$match,ordered=T,levels=c("other_taxon","unclassified_genus_but_consistent_taxonomy","known_host_genus"))
summary(df_knownhost)

## Consistent vs different taxon based on number of hits
## Make counts by number of hits
tmp <- df_knownhost %>%
  group_by(n_hits,match) %>%
  summarise(total=sum(counts))

p1 <- ggplot(tmp) + geom_bar(aes(x=n_hits,y=total,fill=match),stat="identity",position="fill",color="black") + xlab("Number of spacer hits for phage-repeat pair") + ylab("Percentage of observation") + scale_y_continuous(labels=scales::percent) + coord_flip() + blank_theme + scale_fill_manual(values=c("#feb24c","#7fcdbb","#31a354")) + theme(legend.position="top")
p2 <- ggplot(tmp[tmp$match!="unclassified_genus_but_consistent_taxonomy",]) + geom_bar(aes(x=n_hits,y=total,fill=match),stat="identity",position="fill",color="black") + xlab("Number of spacer hits for phage-repeat pair") + ylab("Percentage of observation (ignoring unclassified taxa, unless taxonomies\nare already inconsistent for ranks where both are available)") + scale_y_continuous(labels=scales::percent) + coord_flip() + blank_theme + scale_fill_manual(values=c("#feb24c","#31a354")) + theme(legend.position="top")

## Consistent vs inconsistent when more than 10 hits, and based on profile
## Make counts by profile
tmp2 <- df_knownhost %>%
  filter(n_hits!="1-4" & n_hits!="5-9") %>%
  group_by(profile,match) %>%
  summarise(total=sum(counts))

p3<-ggplot(tmp2) + geom_bar(aes(x=profile,y=total,fill=match),stat="identity",position="fill",color="black") + xlab("Number of spacer hits for phage-repeat pair") + ylab("Percentage of observation") + scale_y_continuous(labels=scales::percent) + coord_flip() + blank_theme + scale_fill_manual(values=c("#feb24c","#7fcdbb","#31a354")) + theme(legend.position="top")
p4<-ggplot(tmp2[tmp2$match!="unclassified_genus_but_consistent_taxonomy",]) + geom_bar(aes(x=profile,y=total,fill=match),stat="identity",position="fill",color="black") + xlab("Number of spacer hits for phage-repeat pair") + ylab("Percentage of observation (ignoring unclassified taxa, unless taxonomies\nare already inconsistent for ranks where both are available)") + scale_y_continuous(labels=scales::percent) + coord_flip() + blank_theme + scale_fill_manual(values=c("#feb24c","#31a354")) + theme(legend.position="top")

## Plot both atop of one another in a single pdf
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp4 <- ggplot_gtable(ggplot_build(p4))
gp2$widths<-gp1$widths
gp2$heights<-gp1$heights
gp3$widths<-gp1$widths
gp3$heights<-gp1$heights
gp4$widths<-gp1$widths
gp4$heights<-gp1$heights
pdf("Known_hosts_vs_hits.pdf",width=8,height=5)
grid.arrange(gp1,gp2,gp3,gp4,ncol=2,nrow=2,heights=c(5,3))
dev.off()


## Stats: 10 and more vs less than 10 hits
for_test <- tmp %>%
  mutate(n_hits=case_when(n_hits == "1-4" ~ "less_than_10", n_hits == "5-9" ~ "less_than_10", n_hits == "10-49" ~ "more_than_10", n_hits == "50-99" ~ "more_than_10", n_hits == "100-and-more" ~ "more_than_10")) %>%
  rename(n_obs=total) %>%
  group_by(n_hits,match) %>%
  summarise(n_obs=sum(n_obs)) %>%
  mutate(total=sum(n_obs)) %>%
  filter(match=="known_host_genus")
prop.test(x = c(for_test[for_test$n_hits=="less_than_10",]$n_obs, for_test[for_test$n_hits=="more_than_10",]$n_obs),
          n = c(for_test[for_test$n_hits=="less_than_10",]$total, for_test[for_test$n_hits=="more_than_10",]$total),
          correct = FALSE)

## Stats: negative vs positive or unknown
for_test <- tmp2 %>%
  mutate(profile=case_when(profile == "negative" ~ "negative", profile == "unknown" ~ "other", profile == "positive" ~ "other")) %>%
  rename(n_obs=total) %>%
  group_by(profile,match) %>%
  summarise(n_obs=sum(n_obs)) %>%
  mutate(total=sum(n_obs)) %>%
  filter(match=="known_host_genus")
prop.test(x = c(for_test[for_test$profile=="negative",]$n_obs, for_test[for_test$profile=="other",]$n_obs),
          n = c(for_test[for_test$profile=="negative",]$total, for_test[for_test$profile=="other",]$total),
          correct = FALSE)
