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
df_setsize<-read.delim("Figures/Main/Fig_2/Spacer_set_size_info.tsv",stringsAsFactors = T)
df_setsize$ecosystem <- factor(df_setsize$ecosystem,ordered=T,levels=levels_ecosystem)
df_setsize$type <- factor(df_setsize$type,ordered=T,levels=levels_type)

## x-y plot showing relationship between maximum coverage (x-axis) and spacer set size (y-axis)
forplot <- df_setsize %>% slice_sample(n=500000)
pdf("Fig_S8A.pdf",width=6,height=6)
ggplot(forplot) + geom_point(aes(x=max_cover,y=spacer_set_size),alpha=0.1) + geom_smooth(aes(x=max_cover,y=spacer_set_size),method="lm") + scale_x_log10() + scale_y_log10() + xlab("maximum spacer coverage in set") + ylab("spacer set size (number of spacers)")
dev.off()

## Similar data but for all (i.e. with 18.1 million spacers) and with the max cover binned
df_setsize$bins_maxcover <- cut(x=df_setsize$max_cover,c(seq(1,10,1),seq(20,100,10),seq(200,1000,100),1E+100),include.lowest = TRUE, right = FALSE)
pdf("Fig_S6B.pdf",width=6,height=3)
ggplot(df_setsize) + geom_boxplot(aes(x=bins_maxcover,y=spacer_set_size),outliers = FALSE) + scale_y_log10() + xlab("maximum spacer coverage in set (binned)") + ylab("spacer set size (number of spacers)") + blank_theme + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
dev.off()

## Look at distribution of spacer set size across coverage bins for different metadata
# Calculate the relative size of each array within these bins
df_setsize <- df_setsize %>%
  group_by(bins_maxcover) %>%
  mutate(fc_median = spacer_set_size / median(spacer_set_size), pc_rank = percent_rank(spacer_set_size))
pdf("Fig_S6C.pdf",width=6,height=3)
ggplot(df_setsize) + geom_boxplot(aes(x=type,y=pc_rank,fill=type)) + coord_flip() + xlab("predicted CRISPR type") + ylab("set size percentile within individual coverage bin") + scale_fill_manual(values=scale_type) + blank_theme + theme(legend.position="bottom") + scale_y_continuous(labels=scales::percent)
dev.off()
pdf("Fig_S6D.pdf",width=6,height=4)
ggplot(df_setsize) + geom_boxplot(aes(x=ecosystem,y=pc_rank,fill=ecosystem)) + coord_flip() + xlab("predicted CRISPR type") + ylab("set size percentile within individual coverage bin") + scale_fill_manual(values=scale_ecosystem) + blank_theme + theme(legend.position="bottom") + scale_y_continuous(labels=scales::percent)
dev.off()
## ANOVA:
res_aov <- aov(pc_rank ~ ecosystem, data = df_setsize)
summary(res_aov)
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
