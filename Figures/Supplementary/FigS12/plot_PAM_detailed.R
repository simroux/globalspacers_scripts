library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))
zero_theme<-theme(axis.line=element_blank(),panel.background=element_blank(),legend.position="none",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))

## Load data: summarize by type, then by spacer
df_pam_type <- read.delim(paste("../../Main/Fig_3/Summary_PAM_detection_by_type.tsv",sep=""),stringsAsFactors = T)
summary(df_pam_type)
df_pam_freq <- read.delim(paste("../../Main/Fig_3/Summary_PAM_detection_by_spacer.tsv",sep=""),stringsAsFactors = T)
summary(df_pam_freq)

## Note: a repeat is associated to a PAM if it is found in > 50% of its neighborhoods (empirical cutoff determined by looking at the global distribution of the number of spacer with a match to the predicted PAM for a given repeat)
### Next, we select which CRISPR type we will use
total <- df_pam_type %>%
  group_by(Type) %>%
  summarise(total=sum(Count)) %>%
  filter(total>=20)  ## Showing only the CRISPR types for which we have at least 30 repeats, otherwise rare types make this impossible to visualize in a nice figure
## We do the sum to get nice bars -> note that we use all counts now
df_pam_type_sum <- df_pam_type %>%
  group_by(Lvl_1_type,Type,Expected) %>%
  summarise(Count=sum(Count)) %>%
  filter(Type %in% total$Type)

df_pam_type_sum$Type<-factor(df_pam_type_sum$Type,ordered=T,levels=rev(c(as.character(unique(total[total$Type!="Unknown",]$Type)),"Unknown")))
df_pam_type_sum$Expected<-factor(df_pam_type_sum$Expected,ordered=T,levels=rev(c("primary","other","no_motif")))

## First plot: how many arrays do we have a PAM for
ggplot(df_pam_type_sum) + geom_bar(aes(x=Type,y=Count,fill=Lvl_1_type,alpha=Expected),col="black",stat="identity",position="fill") + scale_x_discrete(drop=F) + scale_y_continuous(labels=scales::percent, expand=c(0.01,0.01)) + xlab("CRISPR array (predicted) type") + ylab("Percentage of arrays") + coord_flip() + scale_fill_brewer(palette="Set3") + scale_alpha_manual(values=c(0,0.3,1)) + blank_theme

## Second plot: how many spacers have a match to the predicted pam
df_pam_freq_selec <- df_pam_freq %>%
  filter(type %in% total$Type) %>%
  filter(count_neigh >= 10) %>% 
  filter(!is.na(refined_motif))
df_pam_freq_selec$type<-factor(df_pam_freq_selec$type,ordered=T,levels=levels(df_pam_type_sum$Type))

## Should never be below 50% ? To check -- ok
ggplot(df_pam_freq_selec) + geom_boxplot(aes(x=type,y=per_neigh/100,fill=lvl_1_type)) + scale_y_continuous(labels=scales::percent, expand=c(0.01,0.01)) + xlab("CRISPR array (predicted) type") + ylab("Percentage of neighborhoods") + coord_flip()  + scale_fill_brewer(palette="Set3") + blank_theme

## Third plot: show the dominant motif in each type
primary_motif <- df_pam_type %>%
  filter(Expected=="primary") %>%
  group_by(Type) %>%
  summarise(motifs = paste(Motif, collapse=", ")) %>%
  filter(Type %in% total$Type)
primary_motif$Type<-factor(primary_motif$Type,ordered=T,levels=rev(c(as.character(unique(total[total$Type!="Unknown",]$Type)),"Unknown")))
ggplot(primary_motif) + geom_text(aes(x=Type,y=1,label=motifs)) + xlab("CRISPR array (predicted) type") + scale_x_discrete(drop=F) + coord_flip() + scale_fill_brewer(palette="Set3") + scale_alpha_manual(values=c(0,0.3,1)) + zero_theme + theme(axis.ticks.x=element_blank(),axis.line.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank())

## Fourth plot: show the top 3 motifs for each array, if more than 10 observations, 
most_common_motifs <- df_pam_type %>%
  arrange(desc(Count)) %>%
  filter(!is.na(Motif)) %>%
  filter(Count_min_10>=10) %>%
  group_by(Type) %>%
  top_n(3) %>%
  summarise(motifs = paste(Motif, collapse=", ")) %>%
  filter(Type %in% total$Type)
most_common_motifs$Type<-factor(most_common_motifs$Type,ordered=T,levels=rev(c(as.character(unique(total[total$Type!="Unknown",]$Type)),"Unknown")))

## Arrange in grid:
p1<-ggplot(df_pam_type_sum) + geom_bar(aes(x=Type,y=Count,fill=Lvl_1_type,alpha=Expected),col="black",stat="identity",position="fill") + scale_x_discrete(drop=F) + scale_y_continuous(labels=scales::percent, expand=c(0.01,0.01)) + xlab("CRISPR array (predicted) type") + ylab("Percentage of arrays") + coord_flip() + scale_fill_brewer(palette="Set3") + scale_alpha_manual(values=c(0,0.3,1)) + blank_theme + theme(legend.position="bottom")
p2<-ggplot(df_pam_freq_selec) + geom_boxplot(aes(x=type,y=per_neigh/100,fill=lvl_1_type)) + scale_x_discrete(drop=F) + scale_y_continuous(labels=scales::percent, expand=c(0.01,0.01), breaks=c(0.5,0.75,1)) + xlab("CRISPR array (predicted) type") + ylab("Percentage of neighborhoods") + coord_flip()  + scale_fill_brewer(palette="Set3") + blank_theme + theme(legend.position = "none")
p3<-ggplot(primary_motif) + geom_text(aes(x=Type,y=1,label=motifs)) + xlab("CRISPR array (predicted) type") + scale_x_discrete(drop=F) + coord_flip() + scale_fill_brewer(palette="Set3") + scale_alpha_manual(values=c(0,0.3,1)) + zero_theme + theme(axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.line.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank())
p4<-ggplot(most_common_motifs) + geom_text(aes(x=Type,y=1,label=motifs)) + xlab("CRISPR array (predicted) type") + scale_x_discrete(drop=F) + coord_flip() + scale_fill_brewer(palette="Set3") + scale_alpha_manual(values=c(0,0.3,1)) + zero_theme + theme(axis.ticks.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_blank(),axis.line.x=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank())

gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp4 <- ggplot_gtable(ggplot_build(p4))
## Alternative way:
gp2$heights<-gp1$heights
gp3$heights<-gp1$heights
gp4$heights<-gp1$heights
pdf(paste("Fig_PAM_detailed.pdf",sep=""),width=8,height=5)
grid.arrange(gp1, gp3, gp4, gp2, nrow=1,ncol=4, widths=c(1.8,0.5,0.5,1))
dev.off()


