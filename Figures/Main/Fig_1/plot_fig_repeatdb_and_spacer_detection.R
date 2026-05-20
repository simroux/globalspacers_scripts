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
df_repeat<-read.delim("../../../Data/Spacer_db/Array_info_filtered_for_db-Apr23-26.tsv",stringsAsFactors = T)
## Transform lca_origin in repeat source
# This just requires merging "Contig" with "No-information_long" since these are "Long_contig", and renaming "No-information_short" into "Short_contig"
df_repeat$source_repeat <- as.character(df_repeat$lca_origin)
df_repeat[df_repeat$lca_origin=="Contig",]$source_repeat <- "Long_contig"
df_repeat[df_repeat$lca_origin=="No-information_long",]$source_repeat <- "Long_contig"
df_repeat[df_repeat$lca_origin=="No-information_short",]$source_repeat <- "Short_contig"
df_repeat$source_repeat <- factor(df_repeat$source_repeat,ordered=T,levels=rev(c("Genome_high-confidence","Genome_medium-confidence","Genome_low-confidence","Long_contig","Short_contig")))
## Simplify lca_origin
df_repeat <- df_repeat %>% 
  mutate(lca_origin_simp = case_when(str_detect(lca_origin,"^Genome") ~ "Genome", .default="Metagenome_only")) %>%
  mutate(lca_origin_simp = factor(lca_origin_simp,ordered=T,levels=rev(c("Genome","Metagenome_only")))) 
## Simplify CRISPR type
df_repeat <- df_repeat %>%
  separate_wider_delim(type, "-", names = c("type_simple"), too_many = "drop", cols_remove = F) 
df_repeat$type_simple <- factor(df_repeat$type_simple,ordered=T,levels=c("Unknown","VI","V","IV","III","II","I"))
## Simplify taxo to phyla-level + other
df_repeat <- df_repeat %>%
  separate_wider_delim(lca_class, ";", names = c("lca_domain","phylum"), too_many = "drop", cols_remove=FALSE) %>%
  unite("lca_phylum", lca_domain:phylum, remove=FALSE)
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


## How many distinct repeats total in the database, and by origin
df_repeat %>%
  summarise(count=n())
df_repeat %>%
  group_by(lca_origin) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100)
df_repeat %>%
  group_by(lca_origin_simp) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100)

## How many are just short contigs
df_repeat %>%
  filter(lca_origin_simp=="Metagenome_only") %>%
  group_by(lca_origin) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100)

## How many phyla / classes in genome
# Counts for text:
df_repeat %>%
  filter(str_detect(lca_origin,"^Genome")) %>%
  filter(!is.na(lca_phylum)) %>%
  group_by(lca_domain,phylum) %>%
  filter(!grepl("p__unclassified",phylum)) %>%
  filter(!grepl("d__Viruses",lca_domain)) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=500) 
## 159 phyla, since we remove Viruses and unclassified
df_repeat %>%
  filter(str_detect(lca_origin,"^Genome")) %>%
  filter(!is.na(lca_family)) %>%
  group_by(lca_family) %>%
  filter(!grepl("f__unclassified",lca_family)) %>%
  filter(!grepl("d__Viruses",lca_family)) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=5000) 
## 2,438 families
df_repeat %>%
  filter(str_detect(lca_origin,"^Genome")) %>%
  filter(!is.na(lca_genus)) %>%
  group_by(lca_genus) %>%
  filter(!grepl("g__unclassified",lca_genus)) %>%
  filter(!grepl("d__Viruses",lca_genus)) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=50000) 
# 8,135

## Also look at the distribution of CRISPR types in genome vs metagenome
df_repeat %>%
  group_by(lca_origin_simp,type_simple) %>%
  summarise(count=n()) %>%
  mutate(pcent=count/sum(count)*100)
df_repeat %>%
  filter(lca_origin_simp=="Genome") %>%
  group_by(type) %>%
  summarise(count=n()) %>%
  mutate(pcent=count/sum(count)*100) %>%
  arrange(desc(pcent))
# 
df_repeat %>%
  filter(lca_origin_simp=="Genome") %>%
  filter(lca_domain=="d__Archaea") %>%
  group_by(type) %>%
  summarise(count=n()) %>%
  mutate(pcent=count/sum(count)*100) %>%
  arrange(desc(pcent)) %>%
  print(n=50)

## Same without no-info short
df_repeat %>%
  filter(lca_origin!="No-information_short") %>%
  summarise(count=n(), pcent = n() / nrow(.) *100)
df_repeat %>%
  filter(lca_origin!="No-information_short" & lca_origin!="No-information_long") %>%
  group_by(lca_origin) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100)
## df_repeat will be used for the number of distinct repeats, and the breakdown by type and taxonomic assignment


## 
df_lcao_spacer <- read.delim("./lca_origin_to_spacer_count.tsv", stringsAsFactors = T)
## Transform lca_origin in repeat source (same as above)
df_lcao_spacer$source_repeat <- as.character(df_lcao_spacer$lca_origin)
df_lcao_spacer[df_lcao_spacer$lca_origin=="Contig",]$source_repeat <- "Long_contig"
df_lcao_spacer[df_lcao_spacer$lca_origin=="No-information_long",]$source_repeat <- "Long_contig"
df_lcao_spacer[df_lcao_spacer$lca_origin=="No-information_short",]$source_repeat <- "Short_contig"
df_lcao_spacer$source_repeat <- factor(df_lcao_spacer$source_repeat,ordered=T,levels=rev(c("Genome_high-confidence","Genome_medium-confidence","Genome_low-confidence","Long_contig","Short_contig")))
## Make a group by to avoid an ugly line in your bar chart
df_lcao_spacer <- df_lcao_spacer %>%
  group_by(source_repeat) %>%
  summarise(total_spacers=sum(total_spacers))  
## Simplify lca_origin
df_lcao_spacer <- df_lcao_spacer %>% 
  mutate(lca_origin_simp = case_when(str_detect(source_repeat,"^Genome") ~ "Genome", .default="Metagenome_only")) %>%
  mutate(lca_origin_simp = factor(lca_origin_simp,ordered=T,levels=rev(c("Genome","Metagenome_only"))))
summary(df_lcao_spacer)
## Make counts for papers: how many spacers in genomes vs metagenome only
df_lcao_spacer %>%
  group_by(lca_origin_simp) %>%
  summarise(count=sum(total_spacers)) %>%
  mutate(pcent=count/sum(count)*100)
## and how many in short contigs specifically
df_lcao_spacer %>%
  group_by(source_repeat) %>%
  summarise(count=sum(total_spacers)) %>%
  mutate(pcent=count/sum(count)*100)

## And same for the number of sets
df_repeat_detect <-read.delim("./array_to_sample_counts.tsv",stringsAsFactors = T)
## Adjust for the weirdness of SQL (first column repeat_cluster is from the left part of the join so incomplete, and zeroes are NA so we restore them as zeroes)
df_repeat_detect$repeat_cluster <- df_repeat_detect$repeat_cluster_1
df_repeat_detect <- df_repeat_detect %>%  
  mutate(total_sample = replace_na(total_sample, 0))
## Overall counts: we'll divide the the repeats in groups based on prevalence (how many sample) when found in at least 5 reads
df_repeat_detect <- df_repeat_detect %>%
  mutate(sample_groups = cut(total_sample,breaks=c(0,1,10,50,100,1000,999999),right=F)) %>%
  mutate(sample_groups = recode(sample_groups, "[0,1)"="0","[1,10)"="1-9","[10,50)"="10-49","[50,100)"="50-99","[100,1e+03)"="100-999","[1e+03,1e+06)"=">=1000", .default="NA")) 
## Transform lca_origin in repeat source (same as above)
df_repeat_detect$source_repeat <- as.character(df_repeat_detect$lca_origin)
df_repeat_detect[df_repeat_detect$lca_origin=="Contig",]$source_repeat <- "Long_contig"
df_repeat_detect[df_repeat_detect$lca_origin=="No-information_long",]$source_repeat <- "Long_contig"
df_repeat_detect[df_repeat_detect$lca_origin=="No-information_short",]$source_repeat <- "Short_contig"
df_repeat_detect$source_repeat <- factor(df_repeat_detect$source_repeat,ordered=T,levels=rev(c("Genome_high-confidence","Genome_medium-confidence","Genome_low-confidence","Long_contig","Short_contig")))
summary(df_repeat_detect)

## Overall count: how many repeats have 0 sets (in short contigs)
df_repeat_detect %>%
  filter(sample_groups==0) %>%
  group_by(source_repeat) %>%
  summarise(n=n())



########## FINAL PREPARATION PANELS B AND C TOGETHER ###########
# Number of distinct repeats
p1 <- ggplot(df_repeat) + geom_bar(aes(x=source_repeat,fill=lca_origin_simp),col="black",width=0.8) + scale_y_continuous(labels=scales::comma) + xlab("Predicted CRISPR repeat origin and taxonomic assignment confidence") + ylab("Number of predicted CRISPR repeats") + coord_flip() + blank_theme + scale_fill_manual(values=c("gray","black")) + theme(legend.position="bottom")
# Number of spacers associated with these repeats
p2 <- ggplot(df_lcao_spacer) + geom_bar(aes(x=source_repeat,y=total_spacers,fill=lca_origin_simp),col="black",stat="identity",width=0.8) + xlab("Predicted CRISPR repeat origin and taxonomic assignment confidence") + ylab("Number of CRISPR spacers identified across SRA") + coord_flip() + blank_theme + scale_fill_manual(values=c("gray","black")) + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom")
# Distribution of CRISPR types
p3 <- ggplot(df_repeat) + geom_bar(aes(x=source_repeat,fill=type_simple),col="black",position="fill",width=0.8) + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_type) + xlab("Predicted CRISPR type") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom")
# Distribution of taxonomic assignments
p4 <- ggplot(df_repeat) + geom_bar(aes(x=source_repeat,fill=lca_phylum),col="black",position="fill",width=0.8) + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_phylum) + xlab("Taxononomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom")
p5 <- ggplot(df_repeat_detect) + geom_bar(aes(x=source_repeat,fill=sample_groups),position="fill",col="black",width=0.8) + blank_theme + coord_flip() + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=c("#bdbdbd",brewer.pal("YlGnBu",n=5))) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),legend.position="bottom")
# 
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp4 <- ggplot_gtable(ggplot_build(p4))
gp5 <- ggplot_gtable(ggplot_build(p5))
## Make sure everything is lined up
gp1$heights<-gp4$heights
gp2$heights<-gp4$heights
gp3$heights<-gp4$heights
gp5$heights<-gp4$heights
# 
pdf("Panels_B_and_C_distribution.pdf",width=11,height=2.76)
grid.arrange(gp1,gp3,gp4,gp2,gp5,nrow=1,widths=c(1.32,1,1,0.67,1))
dev.off()




################## Panel D
df_sta_med_eco <- read.delim("sample_to_array_counts-summarized_by_ecosystem.tsv",stringsAsFactors = T)
summary(df_sta_med_eco)
df_sta_med_eco <- df_sta_med_eco %>%
  filter(ecosystem!="Unknown") %>%
  filter(ecosystem!="Other") ## Also removing others, not really meaningful, and it's detailed in Fig. S4
df_sta_med_eco$ecosystem <- factor(df_sta_med_eco$ecosystem,ordered=T,levels=levels_ecosystem)
ggplot(df_sta_med_eco,aes(x=ecosystem,y=average_median)) + geom_bar(aes(fill=ecosystem),stat="identity",alpha=1,width=0.75,col="black") + geom_linerange(aes(ymin=min_median, ymax=max_median), col="gray30", alpha=0.8) + scale_fill_manual(values=scale_ecosystem) + xlab("Ecosystem type") + ylab("Median number of CRISPR array per 1Gb") + coord_flip() + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3), legend.position="none")
ggsave("Panel_D_ecosystem_perGb.pdf",width=5.5,height=3.52)
