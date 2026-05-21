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
df_setsize<-read.delim("../../Main/Fig_2/Spacer_set_size_info.tsv",stringsAsFactors = T)
df_setsize$ecosystem <- factor(df_setsize$ecosystem,ordered=T,levels=levels_ecosystem)
df_setsize$type <- factor(df_setsize$type,ordered=T,levels=levels_type)


## Histogram of set size
p1 <- ggplot(df_setsize) + geom_histogram(aes(x=spacer_set_size),bins=50,col="black",fill="gray") + geom_vline(aes(xintercept=median(spacer_set_size)),linetype="dashed",col="red") + scale_x_log10(labels=scales::comma) + scale_y_continuous(labels=scales::comma) + blank_theme + xlab("spacer set size (number of spacers for a given repeat in a given sample)") + theme(panel.grid.major.x=element_line(colour="grey"),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank()) + ylab("number of observations")
## Histogram of set size with only common spacers
p2 <- ggplot(df_setsize) + geom_histogram(aes(x=common_spacer_set_size),bins=15,col="black",fill="gray") + geom_vline(aes(xintercept=median(common_spacer_set_size)),linetype="dashed",col="red") + scale_x_log10(labels=scales::comma) + scale_y_continuous(labels=scales::comma) + blank_theme + xlab("common spacer set size (number of 'common' spacers for a given repeat in a given sample)") + theme(panel.grid.major.x=element_line(colour="grey"),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank()) + ylab("number of observations")

## Correlation between set size and max coverage in set for a subselection of:
# - 12 arrays
# - 100 sets per arrays
# Boxplot of the distribution of correlation values for individual arrays (with > XX sets) 
array_rand <- df_setsize %>%
  group_by(array) %>%
  summarise(array=first(array),count=n()) %>%
  filter(count>=100) %>%
  ungroup() %>%
  slice_sample(n=12)
forplotxy <- df_setsize %>%
  filter(array %in% array_rand$array) %>%
  group_by(array,sample) %>%
  slice_sample(n=100)
p3 <- ggplot(forplotxy) + geom_point(aes(x=max_cover,y=spacer_set_size,col=array),alpha=0.3) + geom_smooth(aes(x=max_cover,y=spacer_set_size,col=array,fill=array),method="lm",alpha=0.2) + scale_x_log10() + scale_y_log10() + xlab("maximum spacer coverage in set") + ylab("spacer set size (number of spacers)") + scale_color_brewer(palette="Set3") + scale_fill_brewer(palette="Set3") + blank_theme + theme(legend.position="bottom",panel.grid.major.x=element_line(colour="grey"),panel.grid.major.y=element_line(colour="grey"),panel.grid.minor=element_blank())

## Plot correlations, which are pre-computed in python (see Readme.md)
df_setsize_corr<-read.delim("correlation_results.tsv",stringsAsFactors = T)
df_setsize_corr <- df_setsize_corr %>%
  filter(!is.na(p_spacer_maxcover)) %>%
  mutate(is_sig = case_when(p_spacer_maxcover<0.05 ~ "yes", .default="no")) %>%
  mutate(is_sig = factor(is_sig))
summary(df_setsize_corr)

## Plot the distribution of correlation coefficients
p4 <- ggplot(df_setsize_corr) + geom_boxplot(aes(x=is_sig,y=r2_spacer_maxcover,fill=is_sig)) + blank_theme + theme(legend.position="none")+ ylab("correlation coefficient") + xlab("statistically significant correlation")
p5 <- ggplot(df_setsize_corr) + geom_bar(aes(x=1,fill=is_sig),col="black",position="fill") + scale_y_continuous(labels=scales::percent) + theme(axis.ticks.x = element_blank(),axis.text.x = element_blank()) + xlab("") + ylab("percentage of repeats") + blank_theme 

gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp4 <- ggplot_gtable(ggplot_build(p4))
gp5 <- ggplot_gtable(ggplot_build(p5))
gp2$widths <- gp1$widths

pdf("Supp_mat/Fig_S8_top_A_to_C.pdf",width=11,height=4)
grid.arrange(
  gp1,gp2,gp3,gp4,gp5,
  widths = c(2, 1.4, 0.75, 0.25),
  layout_matrix = rbind(c(1, 3, 5, NA),
                        c(2, 3, 4, 4))
)
dev.off()

## 
df_setsize$bins_maxcover <- cut(x=df_setsize$max_cover,c(seq(1,10,1),seq(20,100,10),seq(200,1000,100),1E+100),include.lowest = TRUE, right = FALSE)

## Look at distribution of spacer set size across coverage bins for different metadata
# Calculate the relative size of each array within these bins
df_setsize <- df_setsize %>%
  group_by(bins_maxcover) %>%
  mutate(fc_median = spacer_set_size / median(spacer_set_size), pc_rank = percent_rank(spacer_set_size))
pdf("Fig_S8D_left.pdf",width=6,height=3)
ggplot(df_setsize) + geom_boxplot(aes(x=type,y=pc_rank,fill=type)) + coord_flip() + xlab("predicted CRISPR type") + ylab("set size percentile within individual coverage bin") + scale_fill_manual(values=scale_type) + blank_theme + theme(legend.position="bottom") + scale_y_continuous(labels=scales::percent)
dev.off()
pdf("Fig_S8D_right.pdf",width=6,height=4)
ggplot(df_setsize) + geom_boxplot(aes(x=ecosystem,y=pc_rank,fill=ecosystem)) + coord_flip() + xlab("predicted CRISPR type") + ylab("set size percentile within individual coverage bin") + scale_fill_manual(values=scale_ecosystem) + blank_theme + theme(legend.position="bottom") + scale_y_continuous(labels=scales::percent)
dev.off()
#
## Prepare a plot that shows where the selected genera are across the coverage bins
# Look for the genera that are very often / always in the top of the coverage bin in terms of spacer set size
df_genus <- df_setsize %>%
  filter(!is.na(genus)) %>%
  filter (!grepl("g__unclassified$", genus)) %>%
  filter (!grepl("g__$", genus)) %>%
  group_by(genus) %>%
  summarise(n_obs_array = n(), p_high=sum(pc_rank>=0.8)/n()) %>%
  filter(n_obs_array>=200) %>%
  arrange(desc(p_high))
# select top 10, and show the distribution of percentile ranks
selected_genus <- df_genus %>%
  top_n(10) %>%
  pull(genus)
for (i in 1:5){
  ggplot() + geom_boxplot(data=df_setsize,aes(x=bins_maxcover,y=spacer_set_size),outliers = FALSE) + geom_jitter(data=df_setsize[df_setsize$genus %in% selected_genus[((i*2)-1):(i*2)],],aes(x=bins_maxcover,y=spacer_set_size,col=genus),alpha=0.3) + scale_y_log10() + xlab("maximum spacer coverage in set (binned)") + ylab("spacer set size (number of spacers)") + blank_theme + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)) + theme(legend.position="bottom")
  ggsave(paste("Fig_S8E",i,".pdf",sep=""),width=6,height=3)
}
