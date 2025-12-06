library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
source("../../color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))

################## Panel B
df_repeat<-read.delim("../../../Data/Spacer_db/Array_info_filtered_for_db-Oct24-25.tsv",stringsAsFactors = T)
## Get order of lca
df_repeat$lca_origin <- factor(df_repeat$lca_origin,ordered=T,levels=rev(c("Genome_high-confidence","Genome_medium-confidence","Genome_low-confidence","Contig","No-information_long","No-information_short")))
## Simplify taxo to phyla-level + other
df_repeat <- df_repeat %>%
  separate_wider_delim(lca_class, ";", names = c("lca_domain","phylum"), too_many = "drop", cols_remove=FALSE) %>%
  unite("lca_phylum", lca_domain:phylum, remove=FALSE)
## Simplify lca_origin
df_repeat <- df_repeat %>% 
  mutate(lca_origin_simp = case_when(str_detect(lca_origin,"^Genome") ~ "Genome", .default="Metagenome_only")) %>%
  mutate(lca_origin_simp = factor(lca_origin_simp,ordered=T,levels=rev(c("Genome","Metagenome_only")))) 
# Create an "other" category for < 5%
n_affi <- df_repeat %>%
  filter(lca_phylum!="NA_NA") %>%
  summarise(n=n()) %>%
  pull(n)
min_count <- n_affi * 0.01
phylum_to_other <- df_repeat %>%
  filter(lca_phylum!="NA_NA") %>%
  group_by(lca_phylum) %>%
  summarise(n=n()) %>%
  filter(n<min_count) %>%
  pull(lca_phylum)
df_repeat[df_repeat$lca_phylum %in% phylum_to_other,]$lca_phylum<-paste(df_repeat[df_repeat$lca_phylum %in% phylum_to_other,]$lca_domain,";Other",sep="")
df_repeat$lca_phylum <- factor(df_repeat$lca_phylum,ordered=T,levels=c("d__Viruses_p__Uroviricota","d__Archaea;Other","d__Archaea_p__Thermoproteota","d__Archaea_p__Halobacteriota","d__Bacteria_p__unclassified","d__Bacteria;Other","d__Bacteria_p__Verrucomicrobiota","d__Bacteria_p__Pseudomonadota","d__Bacteria_p__Planctomycetota","d__Bacteria_p__Desulfobacterota","d__Bacteria_p__Cyanobacteriota","d__Bacteria_p__Chloroflexota","d__Bacteria_p__Campylobacterota","d__Bacteria_p__Bacteroidota","d__Bacteria_p__Bacillota_C","d__Bacteria_p__Bacillota_A","d__Bacteria_p__Bacillota","d__Bacteria_p__Actinomycetota","d__Bacteria_p__Acidobacteriota"))
summary(df_repeat$lca_phylum)

## How many distinct repeats total in the database, and by origin - to draw circles in panel B
tmp <- df_repeat %>%
  group_by(lca_origin_simp) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  mutate(lca_origin_simp=factor(lca_origin_simp,ordere=T,levels=c("Genome","Metagenome_only")))
ggplot(tmp) + geom_point(aes(x=lca_origin_simp,y=5,fill=lca_origin_simp,size=pcent),pch=21) + scale_fill_manual(values=c("black","grey")) + scale_radius(limits=c(0,100),range=c(1,60)) + blank_theme + theme(legend.position="none")

## Taxonomy breakdown: separately for genomes and metagenomes
ggplot(df_repeat[df_repeat$lca_origin_simp=="Genome",]) + geom_bar(aes(x=1,fill=lca_phylum),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_phylum) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom") + guides(fill=guide_legend(ncol=3))
ggplot(df_repeat[df_repeat$lca_origin_simp=="Metagenome_only",]) + geom_bar(aes(x=1,fill=lca_phylum),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_phylum) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom") + guides(fill=guide_legend(ncol=3))
p21 <- ggplot(df_repeat[df_repeat$lca_origin_simp=="Genome",]) + geom_bar(aes(x=1,fill=lca_phylum),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_phylum) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom") + guides(fill=guide_legend(ncol=3))
p22 <- ggplot(df_repeat[df_repeat$lca_origin_simp=="Metagenome_only",]) + geom_bar(aes(x=1,fill=lca_phylum),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_phylum) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom") + guides(fill=guide_legend(ncol=3))
## Plot the two on one panel
gp21 <- ggplot_gtable(ggplot_build(p21))
gp22 <- ggplot_gtable(ggplot_build(p22))
## Make sure gp21 and gp22 are lined up
gp22$widths<-gp21$widths
gp22$heights<-gp21$heights

pdf("Panel_B-fig.pdf",width=7,height=6)
grid.arrange(gp21,gp22,nrow=2)
dev.off()

################## Panel D
df_repeat_detect <-read.delim("array_to_sample_counts.tsv",stringsAsFactors = T)
## Adjust for the weirdness of SQL (first column repeat_cluster is from the left part of the join so incomplete, and zeroes are NA so we restore them as zeroes)
df_repeat_detect$repeat_cluster <- df_repeat_detect$repeat_cluster_1
df_repeat_detect <- df_repeat_detect %>%  
  mutate(total_sample = replace_na(total_sample, 0))
## Overall counts: we'll divide the the repeats in groups based on prevalence (how many sample) when found in at least 5 reads
df_repeat_detect <- df_repeat_detect %>%
  mutate(sample_groups = cut(total_sample,breaks=c(0,1,10,50,100,1000,999999),right=F)) %>%
  mutate(sample_groups = recode(sample_groups, "[0,1)"="0","[1,10)"="1-9","[10,50)"="10-49","[50,100)"="50-99","[100,1e+03)"="100-999","[1e+03,1e+06)"=">=1000", .default="NA")) # %>%
## Keep the order of the categories as we want them
df_repeat_detect$lca_origin <- factor(df_repeat_detect$lca_origin,ordered=T,levels=levels_taxo)
summary(df_repeat_detect)
## Prep circles
df_spacer_count <- read.delim("lca_origin_to_spacer_count",stringsAsFactors = T)
df_spacer_count <- df_spacer_count %>% 
  mutate(lca_origin_simp = case_when(str_detect(lca_origin,"^Genome") ~ "Genome", .default="Metagenome_only")) %>%
  mutate(lca_origin_simp = factor(lca_origin_simp,ordered=T,levels=rev(c("Genome","Metagenome_only")))) 
tmp <- df_spacer_count %>%
  group_by(lca_origin_simp) %>%
  summarise(count=sum(total_spacers)) %>%
  mutate(pcent = count / sum(count) *100) %>%
  mutate(lca_origin_simp=factor(lca_origin_simp,ordere=T,levels=c("Genome","Metagenome_only")))
ggplot(tmp) + geom_point(aes(x=lca_origin_simp,y=5,fill=lca_origin_simp,size=pcent),pch=21) + scale_fill_manual(values=c("black","grey")) + scale_radius(limits=c(0,100),range=c(1,60)) + blank_theme + theme(legend.position="none")
## Prep bar charts
p31 <- ggplot(df_repeat_detect[df_repeat_detect$lca_origin_simp=="Genome",]) + geom_bar(aes(x=1,fill=sample_groups),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent)+ scale_fill_manual(values=c("#bdbdbd",brewer.pal("YlGnBu",n=5)))  + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom") + guides(fill=guide_legend(ncol=3))
ggplot(df_repeat_detect[df_repeat_detect$lca_origin_simp=="Metagenome_only",]) + geom_bar(aes(x=1,fill=sample_groups),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent)+ scale_fill_manual(values=c("#bdbdbd",brewer.pal("YlGnBu",n=5)))  + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom") + guides(fill=guide_legend(ncol=3))
p32 <- ggplot(df_repeat_detect[df_repeat_detect$lca_origin_simp=="Metagenome_only",]) + geom_bar(aes(x=1,fill=sample_groups),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent)+ scale_fill_manual(values=c("#bdbdbd",brewer.pal("YlGnBu",n=5)))  + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom") + guides(fill=guide_legend(ncol=3))

## Plot the two on one panel
gp31 <- ggplot_gtable(ggplot_build(p31))
gp32 <- ggplot_gtable(ggplot_build(p32))
## Make sure gp21 and gp22 are lined up
gp32$widths<-gp31$widths
gp32$heights<-gp31$heights

pdf("Panel_D-fig.pdf",width=7,height=3.2)
grid.arrange(gp31,gp32,nrow=2)
dev.off()


# Add some numbers
## How many arrays with at least 1 detection
df_repeat_detect %>%
  filter(total_sample>0) %>%
  group_by(repeat_cluster) %>%
  nrow()
## Which type of arrays have 0 detection
df_repeat_detect %>%
  filter(total_sample==0) %>%
  group_by(lca_origin) %>%
  summarise(n=n(), pcent = n() / nrow(.) *100)
## How many repeats with more than 10 detectopms for repeats identified in genomes
df_repeat_detect %>%
  filter(lca_origin=="Genome_high-confidence" | lca_origin=="Genome_medium-confidence" | lca_origin=="Genome_low-confidence") %>%
  mutate(total_sample = case_when(total_sample>=10 ~ "more_than_10", total_sample<10 ~ "less_than_10")) %>%
  group_by(total_sample) %>%
  summarise(n=n(), pcent = n() / nrow(.) *100)

## Somme additional stats around whether spacers are found in multiple arrays and/or genera
df_spacer_array <- read.delim("spacer_to_array_counts.tsv",stringsAsFactors = T)
df_spacer_array %>%
  reframe(n_array=n_array, count=count, pcent = count / sum(count) *100) %>%
  filter(n_array==1)
df_spacer_genus <- read.delim("spacer_to_genus_counts.tsv",stringsAsFactors = T)
df_spacer_genus$n_genus<-as.factor(df_spacer_genus$n_genus)
df_spacer_genus %>%
  filter(n_array>1) %>%
  group_by(n_genus) %>%
  summarise(total=sum(count)) %>%
  reframe(n_genus = n_genus, total = total, pcent = total/sum(total) * 100)

################## Panel C
df_sta_med_eco <- read.delim("sample_to_array_counts-summarized_by_ecosystem.tsv",stringsAsFactors = T)
summary(df_sta_med_eco)
df_sta_med_eco <- df_sta_med_eco %>%
  filter(ecosystem!="Unknown") %>%
  filter(ecosystem!="Other") ## Also removing others, not really meaningful, and it's detailed in Fig. S4
df_sta_med_eco$ecosystem <- factor(df_sta_med_eco$ecosystem,ordered=T,levels=levels_ecosystem)
ggplot(df_sta_med_eco,aes(x=ecosystem,y=average_median)) + geom_bar(aes(fill=ecosystem),stat="identity",alpha=1,width=0.75,col="black") + geom_linerange(aes(ymin=min_median, ymax=max_median), col="gray30", alpha=0.8) + scale_fill_manual(values=scale_ecosystem) + xlab("Ecosystem type") + ylab("Median number of CRISPR array per 1Gb") + coord_flip() + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3), legend.position="none")
ggsave("Panel_C_ecosystem_perGb.pdf",width=5.5,height=3.52)
