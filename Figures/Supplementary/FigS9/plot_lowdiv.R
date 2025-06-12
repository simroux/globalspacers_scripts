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
## Left panels - fraction of low-diversity sets across ecosystems and taxa
df_lowdiv <- read.delim("input_lowdiv_barcharts.Jan3.tsv", stringsAsFactors = T)
summary(df_lowdiv)
df_lowdiv_eco <- df_lowdiv %>%
  filter(category=="ecosystem") %>%
  filter(!is.na(average_sub)) %>%
  mutate(group=factor(group,ordered=T,levels=rev(levels_ecosystem)))
df_lowdiv_tax <- df_lowdiv %>%
  filter(category=="taxon") %>%
  arrange(desc(average_lowdiv)) %>%
  filter(row_number() <=5 | row_number() >= (max(row_number()) - 5)) %>%
  mutate(group=factor(group,ordered=T,levels=rev(unique(group))))
# Prepare scale for class
expanded_purples<-colorRampPalette(brewer.pal(9,"Purples")) # Prep scale
# all
p1<-ggplot(df_lowdiv[df_lowdiv$category=="all",]) + geom_bar(aes(x=group,y=average_sub),col="black",fill="gray",stat="identity") + geom_errorbar(aes(x=group,ymin=average_sub-stdev_sub,ymax=average_sub+stdev_sub),width=0) + scale_y_continuous(labels=scales::percent,lim=c(0,0.16)) + blank_theme + coord_flip() + xlab("") + ylab("Average percentage of low-diversity sets") + theme(panel.grid.major.x=element_line(colour="lightgrey",linewidth=0.3))
# By ecosystem
p2<-ggplot(df_lowdiv_eco) + geom_bar(aes(x=group,y=average_sub,fill=group),col="black",stat="identity") + geom_errorbar(aes(x=group,ymin=average_sub-stdev_sub,ymax=average_sub+stdev_sub),width=0) + scale_y_continuous(labels=scales::percent,lim=c(0,0.16)) + blank_theme + coord_flip() + xlab("") + ylab("Average percentage of low-diversity sets") + theme(legend.position="none",panel.grid.major.x=element_line(colour="lightgrey",linewidth=0.3)) + scale_fill_manual(values=scale_ecosystem)
# By taxon
p3<-ggplot(df_lowdiv_tax) + geom_bar(aes(x=group,y=average_sub,fill=group),col="black",stat="identity") + geom_errorbar(aes(x=group,ymin=average_sub-stdev_sub,ymax=average_sub+stdev_sub),width=0) + scale_y_continuous(labels=scales::percent,lim=c(0,0.16)) + blank_theme + coord_flip() + xlab("") + ylab("Average percentage of low-diversity sets") + theme(legend.position="none",panel.grid.major.x=element_line(colour="lightgrey",linewidth=0.3))  + scale_fill_manual(values=expanded_purples(length(unique(df_lowdiv_tax$group))))
## Plot both atop of one another in a single pdf
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp1$widths<-gp3$widths
gp2$widths<-gp3$widths
pdf("FigsS9_left.pdf",width=5,height=7)
grid.arrange(gp1,gp2,gp3,nrow=3,heights=c(1.25,5,5.5))
dev.off()
### Also stats tests
df_lowdiv_fortest<-read.delim("all_counts.Jan3.tsv")
df_lowdiv_fortest <- df_lowdiv_fortest %>%
  filter(!is.na(n_lowdiv)) %>%
  filter(sampled=="complete")
ref <- df_lowdiv_fortest[1,]
df_result <- data.frame(ref)
df_result$ztestless <- NA
df_result$ztestgreat <- NA
df_result$ztesttwosided <- NA
for (i in 2:nrow(df_lowdiv_fortest)){
  print(i)
  vec <- df_lowdiv_fortest[i,]
  print(vec)
  test <- prop.test(x = c(ref$n_lowdiv, vec$n_lowdiv), n = c(ref$n_obs, vec$n_obs), correct = FALSE, alternative="less")
  vec$ztestless <- test$p.value
  test <- prop.test(x = c(ref$n_lowdiv, vec$n_lowdiv), n = c(ref$n_obs, vec$n_obs), correct = FALSE, alternative="greater")
  vec$ztestgreat <- test$p.value
  test <- prop.test(x = c(ref$n_lowdiv, vec$n_lowdiv), n = c(ref$n_obs, vec$n_obs))
  vec$ztesttwosided <- test$p.value
  df_result<-rbind(df_result,vec)
}
df_result
df_result[df_result$ztestgreat>1e-10 & df_result$ztestless>1e-10,]
write.table(df_result,file="all_counts.Jan3.with_test.tsv",quote=FALSE,sep="\t",row.names=FALSE)
## Middle panel: same with counts ignoring singletons
df_lowdiv <- read.delim("input_lowdiv_barcharts.nosingleton.Jan3.tsv", stringsAsFactors = T)
summary(df_lowdiv)
df_lowdiv_eco <- df_lowdiv %>%
  filter(category=="ecosystem") %>%
  filter(!is.na(average_sub)) %>%
  mutate(group=factor(group,ordered=T,levels=rev(levels_ecosystem)))
df_lowdiv_tax <- df_lowdiv %>%
  filter(category=="taxon") %>%
  arrange(desc(average_lowdiv)) %>%
  filter(row_number() <=5 | row_number() >= (max(row_number()) - 5)) %>%
  mutate(group=factor(group,ordered=T,levels=rev(unique(group))))


# Prepare scale for class
expanded_purples<-colorRampPalette(brewer.pal(9,"Purples")) # Prep scale
# all
p1<-ggplot(df_lowdiv[df_lowdiv$category=="all",]) + geom_bar(aes(x=group,y=average_sub),col="black",fill="gray",stat="identity") + geom_errorbar(aes(x=group,ymin=average_sub-stdev_sub,ymax=average_sub+stdev_sub),width=0) + scale_y_continuous(labels=scales::percent,lim=c(0,0.6)) + blank_theme + coord_flip() + xlab("") + ylab("Average percentage of low-diversity sets") + theme(panel.grid.major.x=element_line(colour="lightgrey",linewidth=0.3))
# By ecosystem
p2<-ggplot(df_lowdiv_eco) + geom_bar(aes(x=group,y=average_sub,fill=group),col="black",stat="identity") + geom_errorbar(aes(x=group,ymin=average_sub-stdev_sub,ymax=average_sub+stdev_sub),width=0) + scale_y_continuous(labels=scales::percent,lim=c(0,0.6)) + blank_theme + coord_flip() + xlab("") + ylab("Average percentage of low-diversity sets") + theme(legend.position="none",panel.grid.major.x=element_line(colour="lightgrey",linewidth=0.3)) + scale_fill_manual(values=scale_ecosystem)
# By taxon
p3<-ggplot(df_lowdiv_tax) + geom_bar(aes(x=group,y=average_sub,fill=group),col="black",stat="identity") + geom_errorbar(aes(x=group,ymin=average_sub-stdev_sub,ymax=average_sub+stdev_sub),width=0) + scale_y_continuous(labels=scales::percent,lim=c(0,0.6)) + blank_theme + coord_flip() + xlab("") + ylab("Average percentage of low-diversity sets") + theme(legend.position="none",panel.grid.major.x=element_line(colour="lightgrey",linewidth=0.3))  + scale_fill_manual(values=expanded_purples(length(unique(df_lowdiv_tax$group))))
## Plot both atop of one another in a single pdf
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp1$widths<-gp3$widths
gp2$widths<-gp3$widths
pdf("Supp_mat/figsS9_mid.pdf",width=5,height=7)
grid.arrange(gp1,gp2,gp3,nrow=3,heights=c(1.25,5,5.5))
dev.off()

## Statistical test
df_lowdiv_fortest<-read.delim("all_counts.nosingleton.Jan3.tsv")
df_lowdiv_fortest <- df_lowdiv_fortest %>%
  filter(!is.na(n_lowdiv)) %>%
  filter(sampled=="complete")
ref <- df_lowdiv_fortest[1,]
df_result <- data.frame(ref)
df_result$ztestless <- NA
df_result$ztestgreat <- NA
df_result$ztesttwosided <- NA
for (i in 2:nrow(df_lowdiv_fortest)){
  print(i)
  vec <- df_lowdiv_fortest[i,]
  print(vec)
  test <- prop.test(x = c(ref$n_lowdiv, vec$n_lowdiv), n = c(ref$n_obs, vec$n_obs), correct = FALSE, alternative="less")
  vec$ztestless <- test$p.value
  test <- prop.test(x = c(ref$n_lowdiv, vec$n_lowdiv), n = c(ref$n_obs, vec$n_obs), correct = FALSE, alternative="greater")
  vec$ztestgreat <- test$p.value
  test <- prop.test(x = c(ref$n_lowdiv, vec$n_lowdiv), n = c(ref$n_obs, vec$n_obs))
  vec$ztesttwosided <- test$p.value
  df_result<-rbind(df_result,vec)
}
df_result
df_result[df_result$ztestgreat>1e-10 & df_result$ztestless>1e-10,]


## Finally, right panel is linking taxa and array types to low diversity spacer sets
## Panel C: Boxplot of the frequency of low-diversity sets when at least 1 low-diversity set 
df_betalowdiv_whole <- read.delim("spacer_sets_frozen.tsv", stringsAsFactors = T)
summary(df_betalowdiv_whole)
df_betalowdiv <- df_betalowdiv_whole %>%
  filter(low_diversity>=1) %>% ## Only take cases with at least 1 low diversity
  mutate(total_qualifying=(low_diversity+standard)) %>%
  filter(total_qualifying>=10) %>% ## At least 10 samples in which the set was not low coverage
  mutate(ratio_low_div=low_diversity/total_qualifying,ratio_standard=standard/total_qualifying) ## Ratio of low diversity among qualifying
df_betalowdiv_taxon <- df_betalowdiv %>%
  filter(!is.na(taxo)) %>%
  separate(taxo,sep=";",c("domain","phylum","class"),extra="drop",remove=FALSE) %>%
  filter(class!="c__unclassified") %>%
  filter(array_type!="Unknown") %>%
  unite("lca_class", domain:class, remove=TRUE, sep=";")
df_plotlowdiv <- df_betalowdiv_taxon %>%
  mutate(code=paste(lca_class,array_type))
# We ask for at least 10 unique repeats in each category to draw a boxplot
select <- df_plotlowdiv %>%
  group_by(code) %>%
  summarise(n=n()) %>%
  filter(n>=10) %>%
  pull(code)
df_plotlowdiv <- df_plotlowdiv %>%
  filter(code %in% select)
ggplot(df_plotlowdiv) + geom_boxplot(aes(x=code,y=ratio_low_div,fill=array_type)) + xlab("") + ylab("Percentage of samples displaying a low diversity of spacers (per array)") + scale_y_continuous(labels=scales::percent) + geom_hline(yintercept=0.25,col="black",linetype="dashed") + coord_flip() + blank_theme + scale_fill_manual(values=scale_type)
ggsave("FigsS9_right.pdf",width=5.3,height=4)