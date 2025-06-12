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
df_repeat<-read.delim("Data/Spacer_db/Array_info_filtered_for_db-Nov1-24.tsv",stringsAsFactors = T)
## Get order of lca
df_repeat$lca_origin <- factor(df_repeat$lca_origin,ordered=T,levels=rev(c("Genome_high-confidence","Genome_medium-confidence","Genome_low-confidence","Contig","No-information_long","No-information_short")))
## Simplify CRISPR type
df_repeat <- df_repeat %>%
  separate_wider_delim(type, "-", names = c("type_simple"), too_many = "drop") 
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
## Same without no-info short
df_repeat %>%
  filter(lca_origin!="No-information_short") %>%
  summarise(count=n(), pcent = n() / nrow(.) *100)
df_repeat %>%
  filter(lca_origin!="No-information_short" & lca_origin!="No-information_long") %>%
  group_by(lca_origin) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100)

## Distribution of taxonomic information for global repeat database
ggplot(df_repeat) + geom_bar(aes(x=1,fill=lca_origin),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_taxo) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
p1<-ggplot(df_repeat) + geom_bar(aes(x=1,fill=lca_origin),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_taxo) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())

## Distribution of CRISPR types for global repeat database
df_repeat %>%
  group_by(type_simple) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=10)
# Count and percentage without unknowns
df_repeat %>%
  filter(type_simple!="Unknown") %>%
  group_by(type_simple) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=10)
# Verify that overall, we mostly don't know the type for very short contigs - IT IS NOT THE CASE
df_repeat %>%
  filter(lca_origin=="No-information_short") %>%
  group_by(type_simple) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=10)
df_repeat %>%
  filter(lca_origin!="No-information_short") %>%
  group_by(type_simple) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=10)
# 
ggplot(df_repeat) + geom_bar(aes(x=1,fill=type_simple),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_type) + xlab("Predicted CRISPR type") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())
p2 <- ggplot(df_repeat) + geom_bar(aes(x=1,fill=type_simple),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_type) + xlab("Predicted CRISPR type") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank()) 

## Distribution across taxa
# Counts for text:
df_repeat %>%
  filter(!is.na(lca_phylum)) %>%
  group_by(lca_domain,phylum) %>%
  filter(!grepl("p__unclassified",phylum)) %>%
  filter(!grepl("d__Viruses",lca_domain)) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=500) 
## We remove Viruses and unclassified
df_repeat %>%
  filter(!is.na(lca_family)) %>%
  group_by(lca_family) %>%
  filter(!grepl("f__unclassified",lca_family)) %>%
  filter(!grepl("d__Viruses",lca_family)) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=5000) 
df_repeat %>%
  filter(!is.na(lca_genus)) %>%
  group_by(lca_genus) %>%
  filter(!grepl("g__unclassified",lca_genus)) %>%
  filter(!grepl("d__Viruses",lca_genus)) %>%
  summarise(count=n(), pcent = n() / nrow(.) *100) %>%
  print(n=50000) 


# Plot
ggplot(df_repeat[!is.na(df_repeat$lca_phylum),]) + geom_bar(aes(x=1,fill=lca_phylum),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_phylum) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank()) + guides(fill=guide_legend(ncol=3))
p3 <- ggplot(df_repeat[!is.na(df_repeat$lca_phylum),]) + geom_bar(aes(x=1,fill=lca_phylum),position="fill",col="black") + blank_theme + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=scale_phylum) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + coord_flip() + theme(axis.ticks.y = element_blank(),axis.text.y = element_blank())  + guides(fill=guide_legend(ncol=3))

## Plot all three on one panel
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
## One way
## Alternative way:
gp1$widths<-gp3$widths
gp2$widths<-gp3$widths

grid.arrange(gp1,gp2,gp3,nrow=3)


pdf("Figures/Fig_1_split/Panel_B-fig.pdf",width=9,height=2)
grid.arrange(gp2,gp1,gp3,nrow=3)
dev.off()
pdf("Figures/Fig_1_split/Panel_B-fig_for_legend.pdf",width=9,height=6)
grid.arrange(gp2,gp1,gp3,nrow=3)
dev.off()

################## Panel C

df_repeat_detect <-read.delim("array_to_sample_counts.tsv",stringsAsFactors = T)
## Adjust for the weirdness of SQL (first column repeat_cluster is from the left part of the join so incomplete, and zeroes are NA so we restore them as zeroes)
df_repeat_detect$repeat_cluster <- df_repeat_detect$repeat_cluster_1
df_repeat_detect <- df_repeat_detect %>%  
  mutate(total_sample = replace_na(total_sample, 0))
## Simplify the CRISPR type
df_repeat_detect <- df_repeat_detect %>%
  separate_wider_delim(type, "-", names = c("type_simple"), too_many = "drop", cols_remove=FALSE) 
df_repeat_detect$type_simple <- factor(df_repeat_detect$type_simple,ordered=T,levels=c("Unknown","VI","V","IV","III","II","I"))
## Overall counts: we'll divide the the repeats in groups based on prevalence (how many sample) when found in at least 5 reads
df_repeat_detect <- df_repeat_detect %>%
  mutate(sample_groups = cut(total_sample,breaks=c(0,1,10,50,100,1000,999999),right=F)) %>%
  mutate(sample_groups = recode(sample_groups, "[0,1)"="0","[1,10)"="1-9","[10,50)"="10-49","[50,100)"="50-99","[100,1e+03)"="100-999","[1e+03,1e+06)"=">=1000", .default="NA")) # %>%
## Keep the order of the categories as we want them
df_repeat_detect$lca_origin <- factor(df_repeat_detect$lca_origin,ordered=T,levels=levels_taxo)
summary(df_repeat_detect)
# ggplot(df_repeat_detect) + geom_bar(aes(x=1,fill=sample_groups),position="fill") + nice_theme + coord_flip()
ggplot(df_repeat_detect) + geom_bar(aes(x=lca_origin,fill=sample_groups),position="fill",col="black") + blank_theme + coord_flip() + scale_y_continuous(labels=scales::percent) + scale_fill_manual(values=c("#bdbdbd",brewer.pal("YlGnBu",n=5))) + xlab("Taxonomic assignment") + ylab("Percentage of CRISPR repeats") + theme(legend.position="bottom")
ggsave("Figures/Fig_1_split/Panel_C_Distribution_sample_n_by_repeat.pdf",width=3.5,height=2.4)

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

################## Panel D
df_sample_detect <-read.delim("sample_to_array_counts.tsv",stringsAsFactors = T)
df_sample_detect$sra_run <- df_sample_detect$sra_run_1
df_sample_detect <- df_sample_detect %>%  
  mutate(total_cluster = replace_na(total_cluster, 0))
df_sample_detect$ecosystem_sum <- factor(df_sample_detect$ecosystem_sum,ordered=T,levels=levels_ecosystem)
## Calculate the number of CRISPR arrays (i.e. unique repeats) normalized for 1E+09 bp, aggregate by ecosystem, and get the 99% confidence interval
new.dat<-data.frame(n_base=1E+09) ## This will be the number for which we will predict this (i.e. 1Gb)
mean_w_err<-data.frame(eco=character(),total_mean=numeric(),total_lower=numeric(),total_upper=numeric(),rsq=numeric(),pval=numeric(),est=numeric(),stderr=numeric())
for (eco in levels(df_sample_detect$ecosystem_sum)){
  print(eco)
  test<-df_sample_detect[df_sample_detect$ecosystem_sum==eco,]
  print(nrow(test))
  lm.model<-lm(total_cluster ~ n_base + 0, data=test) ## Linear model, forcing it through zero
  toto <- summary(lm.model)
  titi <- predict(lm.model, newdata=new.dat, interval='confidence',level=0.99) ## This will give us a predicted value for 1E+09 with a confidence interval for this predicted value. We use "confidence" here, i.e. we have the confidence interval of the mean, not for any individual predicted value (see https://rpubs.com/aaronsc32/regression-confidence-prediction-intervals) ## The alternative is "interval='prediction'"
  mean_w_err<-rbind(mean_w_err,data.frame(eco,titi[1],titi[2],titi[3],toto$r.squared,toto$coefficients[4],toto$coefficients[1],toto$coefficients[2]))
}
colnames(mean_w_err) <- c("eco","mean","lower","upper","rsquared","pvalue","estimate","stderr")
mean_w_err$eco<-factor(mean_w_err$eco,ordered=T,levels=levels_ecosystem)
export_tab <- df_sample_detect %>%
  group_by(ecosystem_sum) %>%
  summarise(n_sample = n()) %>%
  left_join(mean_w_err, by = join_by(ecosystem_sum==eco)) %>%
  arrange(desc(mean))
write.csv(export_tab,file="Supp_mat/Sup_table_number_of_array_per_sample_ecosum.csv",quote=F)
# 
ggplot(mean_w_err,aes(x=eco,y=mean)) + geom_bar(aes(fill=eco),stat="identity",alpha=1,width=0.75,col="black") + geom_linerange(aes(ymin=lower, ymax=upper), col="gray30", alpha=0.8) + scale_fill_manual(values=scale_ecosystem) + xlab("Ecosystem type") + ylab("Average number of CRISPR array (per 1Gb)") + coord_flip() + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3), legend.position="none")
ggsave("Figures/Fig_1_split/Panel_D_ecosystem_perGb.pdf",width=6,height=3)


### Calculate Z-scores to compare thermal springs and engineered to all others regression coefficients
z_scores <- data.frame(eco_1=character(),eco_2=character(),z_stat=numeric(),pvalue=numeric())
for (i in c(4,7)){
  for (j in c(1:3,5:6,8:12)){
    print(paste(i," vs ",j))
    z_res <- abs(mean_w_err[i,"estimate"]-mean_w_err[j,"estimate"])/sqrt(mean_w_err[i,"stderr"]*mean_w_err[i,"stderr"]+mean_w_err[j,"stderr"]*mean_w_err[j,"stderr"])
    z_pv <- 1 - cdf(Z, z_res) + cdf(Z, -z_res)
    z_scores<-rbind(z_scores,data.frame(mean_w_err[i,"eco"],mean_w_err[j,"eco"],z_res,z_pv))
  }
}