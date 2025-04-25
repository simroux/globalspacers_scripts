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
## Taxon vs ecosystem
df_genus_eco<-read.delim("genus_to_ecosystem_counts.tsv",stringsAsFactors = T)
df_genus_eco$ecosystem_sum<-factor(df_genus_eco$ecosystem_sum,ordered=T,levels=levels_ecosystem)
df_genus_eco <- df_genus_eco %>%
  separate_wider_delim(lca_genus, ";", names = c("domain","phylum","class","order","family","genus"), too_many = "drop", cols_remove=FALSE) 
summary(df_genus_eco)
# Remove unclassified genera and unknown ecosystems
df_genus_eco <- df_genus_eco %>%
  filter(ecosystem_sum!="Unknown") %>%
  filter(ecosystem_sum!="Other") %>%
  filter(genus!="g__unclassified")
# Calculate frequency of genera per ecosystem
df_genus_eco <- df_genus_eco %>%
  group_by(ecosystem_sum, lca_genus) %>%
  summarise(n=sum(total_clusters)) %>%
  mutate(freq = n / sum(n) * 100) 
# Prepare the rank for each genus in each ecosystem (and select which ones should be highlighted)
df_genus_eco <- df_genus_eco %>% 
  group_by(ecosystem_sum) %>% 
  mutate(rank = rank(desc(freq), ties.method = "first")) %>%
  mutate(width = case_when(rank<=3 ~ 1, .default = 0))
# Pull all genera that are > 2% OR in the top 3
genus_kept <- df_genus_eco %>%
  filter(freq>=2 | rank <=3) %>%
  group_by(lca_genus) %>%
  pull(lca_genus)
df_genus_eco$lca_genus <- as.character(df_genus_eco$lca_genus)
df_genus_eco[!df_genus_eco$lca_genus %in% genus_kept,]$lca_genus<-"Other"
df_genus_eco$lca_genus <- factor(df_genus_eco$lca_genus,ordered=T)
# Make it ordered in the "right" way, i.e. alphabetical from top to bottom, not bottom to top
df_genus_eco$lca_genus <- factor(df_genus_eco$lca_genus,ordered=T,levels=rev(levels(df_genus_eco$lca_genus)))
# Make up extra values for zeroes (otherwise missing and it's ugly)
df_genus_eco <- df_genus_eco %>%
  ungroup() %>%
  complete(ecosystem_sum,lca_genus,fill = list(n = NA, freq = NA, rank=NA, width=NA)) 
# Some summary to check
summary(df_genus_eco$lca_genus)
summary(df_genus_eco)
# Remove the 'other' category so that we don't have to deal with this that is not really a genus or an ecosystem anyway
df_genus_eco <- df_genus_eco %>%
  filter(lca_genus!="Other") %>%
  filter(ecosystem_sum!="Other") %>%
  filter(ecosystem_sum!="Unknown")
# Set anything under 0.05% as 0 for clarity
df_genus_eco <- df_genus_eco %>%
  mutate(freq = case_when(freq>=0.05 ~ freq, freq < 0.05 ~ NA ))
# Get a prettier taxonomy string for selected genera
df_genus_eco <- df_genus_eco %>%
  separate_wider_delim(lca_genus, ";", names = c("domain","phylum","class","order","family","genus"), too_many = "drop", cols_remove=FALSE) 
df_genus_eco$pretty_genus <- as.factor(paste(df_genus_eco$domain,df_genus_eco$phylum,df_genus_eco$genus,sep=";"))
# Finally, we get a nice plot
ggplot(df_genus_eco,aes(x=ecosystem_sum,y=pretty_genus)) + geom_tile(aes(fill=freq)) + geom_tile(aes(linewidth=width,col=width_col),alpha=0,col="black") + blank_theme + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)) + scale_fill_distiller(palette="GnBu",direction=1,name="Frequency in ecosystem (%)",na.value="white") + scale_linewidth(range=c(0,1)) + guides(linewidth = "none") + theme(legend.position="bottom") + ylab("CRISPR repeat - predicted genus") + xlab("Ecosystem type")
## Test palettes: GnBu
ggsave("Fig_S5_draft_Genus_vs_ecosystem_vs_array_frequency.pdf",width=6,height=7)

