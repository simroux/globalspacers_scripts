library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(stringr)
# loadfonts() ## Not required every time
source("color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))

## Panel A - requires loading of 1Gb+ table
## Load VR data
df_vrhits<-read.delim("fig_hits_input_vr.tsv",stringsAsFactors = T)
df_vrhits$hit <- as.character(df_vrhits$hit)
df_vrhits[is.na(df_vrhits$hit),]$hit <- "no"
df_vrhits[df_vrhits$hit==1,]$hit <- "yes"
df_vrhits$hit <- factor(df_vrhits$hit,ordered=T,levels=c("no","yes"))
summary(df_vrhits$hit)
## Load PR data
df_prhits<-read.delim("fig_hits_input_pr.tsv",stringsAsFactors = T)
df_prhits$hit <- as.character(df_prhits$hit)
df_prhits[is.na(df_prhits$hit),]$hit <- "no"
df_prhits[df_prhits$hit==1,]$hit <- "yes"
df_prhits$hit <- factor(df_prhits$hit,ordered=T,levels=c("no","yes"))
summary(df_prhits$hit)
## Count and plot
df_stats <- data.frame()
# Count all hits
tmp <- df_vrhits %>%
  group_by(hit) %>%
  summarise(count=n()) %>%
  mutate(category=c("vr_all"))
df_stats <- rbind(df_stats,tmp)
# Count riboviria hits
tmp <- df_vrhits %>%
  group_by(hit) %>%
  filter(str_detect(class, "^r__Riboviria")) %>%
  summarise(count=n()) %>%
  mutate(category=c("vr_riboviria"))
df_stats <- rbind(df_stats,tmp)
# Count all non-caudo phage hits
tmp <- df_vrhits %>%
  filter(host_from_taxo=="phage" | host_from_taxo=="archaea") %>%
  filter(class!="r__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;c__Caudoviricetes") %>%
  group_by(hit) %>%
  summarise(count=n()) %>%
  mutate(category=c("vr_phage_non-caudovirus"))
df_stats <- rbind(df_stats,tmp)
# Count all non-caudo phage hits
tmp <- df_vrhits %>%
  filter(class=="r__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;c__Caudoviricetes") %>%
  group_by(hit) %>%
  summarise(count=n()) %>%
  mutate(category=c("vr_caudovirus"))
df_stats <- rbind(df_stats,tmp)
# Count all caudo HQ hits
tmp <- df_vrhits %>%
  filter(class=="r__Duplodnaviria;k__Heunggongvirae;p__Uroviricota;c__Caudoviricetes") %>%
  filter(quality=="High-quality" | quality=="Reference") %>%
  group_by(hit) %>%
  summarise(count=n()) %>%
  mutate(category=c("vr_hq_caudovirus"))
df_stats <- rbind(df_stats,tmp)
# Count all hits for plasmids
tmp <- df_prhits %>%
  group_by(hit) %>%
  summarise(count=n()) %>%
  mutate(category=c("pr_all"))
df_stats <- rbind(df_stats,tmp)
# Count all near-complete hits
tmp <- df_prhits %>%
  filter(putatively_complete=="Yes") %>%
  group_by(hit) %>%
  summarise(count=n()) %>%
  mutate(category=c("pr_near_complete"))
df_stats <- rbind(df_stats,tmp)
# Count all mpf hits
tmp <- df_prhits %>%
  filter(has_mpf=="yes") %>%
  group_by(hit) %>%
  summarise(count=n()) %>%
  mutate(category=c("pr_has_mpf"))
df_stats <- rbind(df_stats,tmp)
df_stats$category<-factor(df_stats$category,ordered=T,levels=rev(c("vr_all","vr_riboviria","vr_phage_non-caudovirus","vr_caudovirus","vr_hq_caudovirus","pr_all","pr_near_complete","pr_has_mpf")))
## Export plots
pdf("Panel_A.pdf",width=3,height=2.5)
ggplot(df_stats) + geom_bar(aes(x=category,y=count,fill=hit),stat="identity",col="black",position="fill") + xlab("plasmid category") + ylab("percentage of potential targets") + scale_y_continuous(labels=scales::percent) + blank_theme + scale_fill_manual(values=c("white","darkgrey")) + coord_flip() + theme(legend.position="none")
dev.off()

## Panel B - PAM // Summary_PAM_detection_by_type.tsv
df_pam_type <- read.delim("Summary_PAM_detection_by_type.tsv",stringsAsFactors = T)
summary(df_pam_type)
tmp <- df_pam_type %>%
  group_by(Lvl_1_type,Expected) %>%
  summarise(total=sum(Count)) %>%
  mutate(Lvl_1_type=factor(Lvl_1_type,ordered=T,levels=levels_type))

ggplot(tmp) + geom_bar(aes(x=Lvl_1_type,y=total,fill=Lvl_1_type,alpha=Expected),col="black",stat="identity",position="fill") + scale_x_discrete(drop=F) + scale_y_continuous(labels=scales::percent, expand=c(0.01,0.01)) + xlab("CRISPR array (predicted) type") + ylab("Percentage of arrays") + coord_flip() + scale_fill_manual(values=scale_type) + scale_alpha_manual(values=c(0,0.3,1)) + blank_theme + theme(legend.position="bottom")
ggsave("Panel_B_PAM.pdf",width=3,height=4)

# Overall stats
df_pam_type %>%
  group_by(Expected) %>%
  summarise(count=sum(Count)) %>%
  mutate(percent=count/sum(count)*100) %>%
  mutate(inverse=100-percent)


## Panel C - Hit by spacer type and ecosystem
df_sp_typ_hit <- read.delim("counts_stats_clean.tsv",stringsAsFactors = T)
df_sp_typ_hit$type<-factor(df_sp_typ_hit$type,ordered=T,levels=c("none","pr","vr","vr_and_pr"))
hit_type_scale<-c("none"="white","pr"="#d6604d","vr"="#4393c3","vr_and_pr"="#ffc13b")
## Numbers: total hit rare vs common, single vs multiple samples
df_sp_typ_hit %>%
  filter(category=="sp_alpha" | category=="sp_beta") %>%
  mutate(type=case_when(type!="none" ~ "hit", type=="none" ~ "none")) %>% ## Transform all the not "none" into hits
  group_by(value,type) %>%
  summarise(total = sum(count)) %>%
  ungroup(type) %>%
  mutate(pcent=total/sum(total))
## Numbers: overall human vs other (mentioned in text)
df_sp_typ_hit_plot_eco %>%
  mutate(type=case_when(type=="none" ~ "none", .default="hit"))  %>%
  group_by(value,type) %>%
  summarise(total = sum(count)) %>%
  mutate(pcent=total/sum(total)*100)
## Prefilter a bit and simplify the input
## .. for spacer type
df_sp_typ_hit_plot_sptype <- df_sp_typ_hit %>%
  filter((category=="sp_alpha" & value=="rare") | (category=="sp_alpha_and_beta" & (value=="common_multi_sample" | value=="common_single_sample"))) %>%
  mutate(value=factor(value,ordered=TRUE,levels=rev(c("rare","common_single_sample","common_multi_sample"))))
## .. for ecosystem type
df_sp_typ_hit_plot_eco <- df_sp_typ_hit %>%
  filter(category=="eco") %>%
  mutate(value=case_when(value=="Human-associated_Digestive-system" ~ "Human", value=="Human-associated_Other" ~ "Human", value=="mixed" ~ "none", .default="Other"))  %>%
  filter(value!="none") %>%
  mutate(value=factor(value,ordered=TRUE,levels=c("Other","Human"))) %>%
  group_by(value,type) %>%
  summarise(count=sum(count))
## Put all in the same data frame to get a single bar chart in the end
df_sp_typ_hit_plot <- df_sp_typ_hit_plot_sptype %>%
  select("value","type","count") %>%
  mutate(value=as.character(value)) 
df_sp_typ_hit_plot <- df_sp_typ_hit_plot_eco %>%
  mutate(value=as.character(value)) %>%
  bind_rows(df_sp_typ_hit_plot)
df_sp_typ_hit_plot$value <- factor(df_sp_typ_hit_plot$value,ordered=T,levels=c("Other","Human","common_multi_sample","common_single_sample","rare"))
ggplot(df_sp_typ_hit_plot) + geom_bar(aes(x=value,y=count,fill=type),col="black",stat="identity",position="fill") + blank_theme + coord_flip() + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=hit_type_scale) + theme(legend.position="none") + ylab("Percentage of spacers") + xlab("Spacer/Array/Ecosystem type")  
ggsave("Panel_C_spacer_type_ecosystem.pdf",width=4,height=2.25)

## Panel D - Number of distinct high-quality vOTU per spacer
df_dist_nvotu <- read.delim("multihit_votu_totaldistrib_hqonly.tsv",stringsAsFactors = T)
df_dist_nvotu$category <- factor(df_dist_nvotu$category,ordered=T,levels=unique(df_dist_nvotu$category))
ggplot(df_dist_nvotu) + geom_bar(aes(x=category,y=count),stat="identity",col="black",fill="grey") + scale_y_continuous(labels=scales::comma,breaks=c(0,5000000,10000000)) + xlab("Number of distinct vOTUs targeted") + ylab("Number of unique spacers") + blank_theme + theme(axis.text.x=element_text(angle=60,vjust=1,hjust=1))
ggsave("panel_D_spacer_nvotu.pdf",width=2.75,height=1.5)

### Panel E: x-y plot showing the maximum distance between pairs of high-quality vOTU matched by the same spacer
library(ggrastr)
library(grid)
df_hit_xy <- read.delim("multihit_votu_vs_distance_sample_hqonly.tsv")
summary(df_hit_xy)
df_hit_xy <- df_hit_xy %>%
  filter(n_hq_uvig>1)
df_hit_xy_dots <- df_hit_xy %>%
  sample_n(100000)
max_x<-max(df_hit_xy$n_hq_uvig)
p2 <- ggplot(df_hit_xy) + geom_histogram(aes(x=n_hq_uvig),col="black",bins=20) + scale_x_log10(breaks=c(2,10,100,1000),lim=c(1.5,max_x)) + scale_y_continuous(labels=scales::comma) + xlab("") + ylab("Number of observations") + blank_theme
p3 <- ggplot(df_hit_xy) + geom_histogram(aes(x=max_dist),col="black") + scale_x_continuous(labels=scales::percent) + scale_y_continuous(labels=scales::comma) + xlab("") + ylab("Number of observations") + blank_theme + coord_flip(xlim=c(-0.018,1.018)) 
p1 <- ggplot(df_hit_xy_dots) + rasterise(geom_point(aes(x=n_hq_uvig,y=max_dist),col="black",alpha=0.01),dpi=600) + scale_x_log10(breaks=c(2,10,100,1000),lim=c(1.5,max_x)) + scale_y_continuous(labels=scales::percent,lim=c(-0.018,1.018)) + xlab("Number of targeted viruses (high-quality genomes only)")  + ylab("Maximum distance between targets (ANI %)") + blank_theme  + theme(panel.grid.major.x=element_line(colour="lightgrey",linewidth=0.3))
##
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp2$widths<-gp1$widths
gp3$heights<-gp1$heights
blank <- grid.rect(gp=gpar(col="white"))
pdf("panel_E_xyplot_dist-vs-nhits.pdf",width=5,height=5)
grid.arrange(gp2,blank,gp1,gp3,nrow=2,ncol=2,widths=c(4,1),heights=c(1,4))
dev.off()