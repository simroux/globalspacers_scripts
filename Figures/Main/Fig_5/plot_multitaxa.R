library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
# loadfonts()
source("../../color_scales.R")
nice_theme<-theme(text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))
nice_theme_vertical<-theme(text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),legend.key.size = unit(1,"line"))
# Nice stripped theme
nice_theme_blank<-theme_classic() + theme(text=element_text(color="black",size=6,family="Source Sans 3"),axis.line = element_blank(),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=6,family="Source Sans 3"))

## Panel A - 3 breakdown of multi-taxa
scale_broad_cat<-c("single_taxon"="#3288bd","multi_genus"="#006837","multi_family"="#d9f0a3","multi_order"="#fed976","multi_class"="#f46d43","multi_phylum"="#d7301f","multi_domain"="#000000","mixed"="#e0e0e0","no_high"="#e0e0e0","no_positive"="#ffffff")
df_base<-read.delim("Counts_for_multitaxa_frequency.tsv",stringsAsFactors = T)
df_base$type<-factor(df_base$type,ordered=T,levels=rev(c("hc","hc_high","hc_high_mis")))
df_base$cat<-factor(df_base$cat,ordered=T,levels=c("mixed","no_high","no_positive","single_taxon","multi_genus","multi_family","multi_order","multi_class","multi_phylum","multi_domain"))
summary(df_base)
ggplot(df_base) + geom_bar(aes(x=type,y=count,fill=cat),col="black",stat="identity",position="fill") + blank_theme + scale_y_continuous(labels=scales::percent) + xlab("") + ylab("Percentage of uvigs") + coord_flip() + scale_fill_manual(values=scale_broad_cat)
ggsave("panel_multitaxa_frequency.pdf",width=6,height=2.5)

## Panel B - link to metadata
df_meta<-read.delim("Multitaxa_features.tsv",stringsAsFactors = T)
## First, we select only the categories we are interested in, and we transform everything into ordered factors
df_meta_filt <- df_meta %>%
  filter(category=="single_class" | category=="one_class_high_other_not_high" | category=="one_class_highpos_other_class_high" | category=="multi_class_with_high_pos") %>%
  mutate(category=factor(category,ordered=T,levels=rev(c("single_class","one_class_high_other_not_high","one_class_highpos_other_class_high","multi_class_with_high_pos")))) %>%
  mutate(lifestyle=case_when(lifestyle=="temperate" ~ "temperate", .default="other")) %>%
  mutate(lifestyle=factor(lifestyle,ordered=T,levels=c("other","temperate"))) %>%
  mutate(dgr=factor(dgr,ordered=T,levels=c("no","yes"))) %>%
  mutate(eco_simp=case_when(ecosystem=="Human-associated_Digestive-system" ~ "Human", ecosystem=="Human-associated_Other" ~ "Human", .default="Other")) %>%
  mutate(eco_simp=factor(eco_simp,ordered=T,levels=c("Other","Human")))
## Next, we will work only with the HQ because otherwise absence of DGR / temperate / Acr not as reliable
df_meta_filt_hq <- df_meta_filt %>%
  filter(quality=="Reference" | quality=="High-quality")
## Then we will do a consensus by vOTU with a 2/3rd majority rule
smart_consensus <- function(x) {
  ux <- unique(x)
  if (length(ux)==1){
    toString(ux)
  }
  else{
    test <- table(x)
    max <- max(test)
    total <- sum(test)
    ratio <- max/total
    if (ratio>=0.66){
      names(test)[which.max(test)]
    }
    else{
      "mixed"
    }
  }
}
df_meta_votu_hq <- df_meta_filt_hq %>%
  group_by(votu) %>%
  summarise(category=smart_consensus(category),quality=smart_consensus(quality),lifestyle=smart_consensus(lifestyle),fibers=smart_consensus(fibers),acr=smart_consensus(acr),eco_simp=smart_consensus(eco_simp),dgr=smart_consensus(dgr),total=n())
## Check how many time we had a consensus vs not for the main parameters we are interested in
df_meta_votu_hq %>%
  group_by(category) %>%
  summarise(total=n()) %>%
  mutate(ratio=total/sum(total)*100) %>%
  arrange(ratio)
# For categories, it's only 0.3%
df_meta_votu_hq %>%
  group_by(dgr) %>%
  summarise(total=n()) %>%
  mutate(ratio=total/sum(total)*100) %>%
  arrange(ratio)
# same for DGR
df_meta_votu_hq %>%
  group_by(lifestyle) %>%
  summarise(total=n()) %>%
  mutate(ratio=total/sum(total)*100) %>%
  arrange(ratio)
# 1.18% for lifestyle, still ok
df_meta_votu_hq %>%
  group_by(eco_simp) %>%
  summarise(total=n()) %>%
  mutate(ratio=total/sum(total)*100) %>%
  arrange(ratio)
# and 0.9% for ecosystem, so we can go ahead
## 
df_meta_votu_hq_forplot <- df_meta_votu_hq %>%
  filter(category!="mixed") %>% ## For pretty plots, we avoid cases where we could not find a 2/3rd majority for the targeting category, this is relatively rare (see above)
  mutate(category=factor(category,ordered=T,levels=rev(c("single_class","one_class_high_other_not_high","one_class_highpos_other_class_high","multi_class_with_high_pos")))) %>%
  mutate(dgr=factor(dgr,ordered=T,levels=c("no","yes"))) %>%
  mutate(eco_simp=factor(eco_simp,ordered=T,levels=c("Other","Human"))) %>%
  mutate(lifestyle=factor(lifestyle,ordered=T,levels=c("other","temperate")))
## For each plot, for simplicity, we ignore the mixed features (again, they are relatively rare, see above) - because we switched to factors, these mixed are now NA
p1 <- ggplot(df_meta_votu_hq_forplot[!is.na(df_meta_votu_hq_forplot$lifestyle),]) + geom_bar(aes(x=category,fill=lifestyle),col="black",position="fill") + coord_flip() + blank_theme + scale_y_continuous(labels=scales::percent) + ylab("Percentage of viruses") + scale_fill_manual(values=c("#ffffff","#fc8d62"))
p2 <- ggplot(df_meta_votu_hq_forplot[!is.na(df_meta_votu_hq_forplot$eco_simp),]) + geom_bar(aes(x=category,fill=eco_simp),col="black",position="fill") + coord_flip() + blank_theme + scale_y_continuous(labels=scales::percent) + ylab("Percentage of viruses") + scale_fill_manual(values=c("#ffffff","#7570b3"))
p3 <- ggplot(df_meta_votu_hq_forplot[!is.na(df_meta_votu_hq_forplot$dgr),]) + geom_bar(aes(x=category,fill=dgr),col="black",position="fill") + coord_flip() + blank_theme + scale_y_continuous(labels=scales::percent) + ylab("Percentage of viruses") + scale_fill_manual(values=c("#ffffff","#66c2a5"))
## Prepare the plot
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp2$widths<-gp1$widths
gp3$widths<-gp1$widths
pdf("multitaxa_features.pdf",width=4.5,height=2.8)
grid.arrange(gp1,gp2,gp3,nrow=3)
dev.off()


## Confirming that the percentage of DGR is higher in categories 2 and 3
for_test <- df_meta_votu_hq_forplot %>%
  filter(!is.na(df_meta_votu_hq_forplot$dgr)) %>%
  group_by(category,dgr) %>%
  summarise(n_obs=n()) %>%
  mutate(total=sum(n_obs)) %>%
  filter(dgr=="yes")
## one_class_high_other_not_high
prop.test(x = c(for_test[for_test$category=="one_class_high_other_not_high",]$n_obs, for_test[for_test$category=="single_class",]$n_obs),
          n = c(for_test[for_test$category=="one_class_high_other_not_high",]$total, for_test[for_test$category=="single_class",]$total),
          correct = FALSE)
prop.test(x = c(for_test[for_test$category=="one_class_high_other_not_high",]$n_obs, for_test[for_test$category=="multi_class_with_high_pos",]$n_obs),
          n = c(for_test[for_test$category=="one_class_high_other_not_high",]$total, for_test[for_test$category=="multi_class_with_high_pos",]$total),
          correct = FALSE)
## one_class_highpos_other_class_high
prop.test(x = c(for_test[for_test$category=="one_class_highpos_other_class_high",]$n_obs, for_test[for_test$category=="single_class",]$n_obs),
          n = c(for_test[for_test$category=="one_class_highpos_other_class_high",]$total, for_test[for_test$category=="single_class",]$total),
          correct = FALSE)
prop.test(x = c(for_test[for_test$category=="one_class_highpos_other_class_high",]$n_obs, for_test[for_test$category=="multi_class_with_high_pos",]$n_obs),
          n = c(for_test[for_test$category=="one_class_highpos_other_class_high",]$total, for_test[for_test$category=="multi_class_with_high_pos",]$total),
          correct = FALSE)
