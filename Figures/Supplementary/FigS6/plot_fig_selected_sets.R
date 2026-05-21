library(ggplot2)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
source("../../color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3 Medium"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3 Medium"),legend.key.size = unit(1,"line"))
## Load data
df_betadiv_desc<-read.delim("input_description_betadiv.tsv",stringsAsFactors = T)
## Original summary
summary(df_betadiv_desc)
## Ordering taxa
df_betadiv_desc$taxonomy<-as.character(df_betadiv_desc$taxonomy)
df_betadiv_desc[is.na(df_betadiv_desc$taxonomy),]$taxonomy<-"unknown"
df_betadiv_desc$taxonomy<-factor(df_betadiv_desc$taxonomy)
df_betadiv_desc$taxonomy<-factor(df_betadiv_desc$taxonomy,ordered=T,levels=rev(levels(df_betadiv_desc$taxonomy)))
## Selecting types to display
df_betadiv_desc <- df_betadiv_desc %>%
  separate_wider_delim(type, "-", names = c("type_simple"), too_many = "drop", cols_remove=F) 
df_betadiv_desc$type_simple <- factor(df_betadiv_desc$type_simple,ordered=T,levels=rev(c("Unknown","VI","V","IV","III","II","I")))
# Displaying sub-types when they represent > 10% of the repeats, type otherwise
selected <- df_betadiv_desc %>%
  group_by(type) %>%
  summarise(count=n()) %>%
  mutate(pcent=count/sum(count)) %>%
  filter(pcent>=0.1) %>%
  select(type)
df_betadiv_desc <- df_betadiv_desc %>%
  mutate(type_for_plot = case_when(type %in% selected$type ~ type, .default = paste(type_simple,sep="")))
unique(df_betadiv_desc$type_for_plot)
# Looks like it only adds I-C and I-E to the display
df_betadiv_desc[df_betadiv_desc$type_for_plot=="I",]$type_for_plot<-"I-other"
df_betadiv_desc$type_for_plot <- factor(df_betadiv_desc$type_for_plot,ordered=T,levels=rev(c("Unknown","VI","V","IV","III","II","I-other","I-E","I-C")))
## Ordering ecosystem
df_betadiv_desc$ecosystem<-factor(df_betadiv_desc$ecosystem,ordered=T,levels=levels_ecosystem)
## Final summary
summary(df_betadiv_desc)

## Panel 0 : just the y axis
tmp <- df_betadiv_desc %>%
  group_by(taxonomy) %>%
  summarise(n_repeats=length(unique(repeat.)),n_samples=sum(n_samples),n_spacers=sum(n_spacers))
p0 <- ggplot(tmp) + geom_point(aes(x=taxonomy,y=0),stat="identity") + coord_flip() + blank_theme + xlab("")
p0 
gp0 <- ggplot_gtable(ggplot_build(p0))

## Panel 1: bar graph of number of repeats per taxon
p1<-ggplot(tmp) + geom_bar(aes(x=taxonomy,y=n_repeats),stat="identity") + scale_y_log10(labels=scales::comma) + coord_flip() + blank_theme + theme(panel.grid.major.x=element_line(colour="grey"),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.text.y.left=element_blank()) + ylab("number of repeats") + xlab("")
p1
gp1 <- ggplot_gtable(ggplot_build(p1))
# ggplot(tmp) + geom_bar(aes(x=taxonomy,y=n_spacers),stat="identity") + coord_flip()
## Panel 3: number of samples
p3<-ggplot(tmp) + geom_bar(aes(x=taxonomy,y=n_samples),stat="identity") + scale_y_log10(breaks=c(1000,10000,100000),labels=scales::comma) + coord_flip() + blank_theme + theme(panel.grid.major.x=element_line(colour="grey"),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.text.y.left=element_blank()) + ylab("number of samples") + xlab("")
p3
gp3 <- ggplot_gtable(ggplot_build(p3))
## Panel 5: number of spacers
p5<-ggplot(tmp) + geom_bar(aes(x=taxonomy,y=n_spacers),stat="identity") + scale_y_log10(breaks=c(100000,1000000,10000000),labels=scales::comma) + coord_flip() + blank_theme + theme(panel.grid.major.x=element_line(colour="grey"),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.text.y.left=element_blank()) + ylab("number of spacers") + xlab("")
p5
gp5 <- ggplot_gtable(ggplot_build(p5))

## Panel 2: distribution across type
tmp <- df_betadiv_desc %>%
  group_by(taxonomy,type_for_plot) %>%
  summarise(n_repeats=length(unique(repeat.)),n_spacers=sum(n_spacers),n_samples=sum(n_samples)) %>%
  mutate(r_repeats=n_repeats/sum(n_repeats)*100)
p2<-ggplot(tmp) + geom_tile(aes(x=type_for_plot,y=taxonomy,fill=r_repeats)) + blank_theme + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.text.y.left=element_blank()) + scale_fill_distiller(palette="GnBu",direction=1,name="Frequency in taxon (%)",na.value="white") + scale_linewidth(range=c(0,1)) + guides(linewidth = "none") + theme(legend.position="bottom") + xlab("CRISPR repeat - predicted type") + ylab("Repeat taxonomic assignment") + ylab("")
p2
gp2 <- ggplot_gtable(ggplot_build(p2))

## Panel 4: distribution across ecosystems
tmp <- df_betadiv_desc %>%
  group_by(taxonomy,ecosystem) %>%
  summarise(n_spacers=sum(n_spacers),n_samples=sum(n_samples)) %>%
  mutate(r_samples=n_samples/sum(n_samples)*100)
p4<-ggplot(tmp) + geom_tile(aes(x=ecosystem,y=taxonomy,fill=r_samples)) + blank_theme + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),axis.text.y.left=element_blank()) + scale_fill_distiller(palette="YlOrRd",direction=1,name="Frequency in taxon (%)",na.value="white") + scale_linewidth(range=c(0,1)) + guides(linewidth = "none") + theme(legend.position="bottom") + xlab("Sample ecosystem") + ylab("Repeat taxonomic assignment") + ylab("")
p4
gp4 <- ggplot_gtable(ggplot_build(p4))

gp4$heights<-gp2$heights
gp0$heights<-gp2$heights
gp1$heights<-gp2$heights
gp3$heights<-gp2$heights
gp5$heights<-gp2$heights

# pdf("Betadiv/test.pdf",width=12,height=5.6)
pdf("Fig_S6_bottom_draft.pdf",width=12,height=5.6)
grid.arrange(gp0,gp1,gp2,gp3,gp4,gp5,nrow=1,ncol=6,widths=c(2,0.3,0.8,0.35,0.8,0.4))
dev.off()

