library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))

out_file<-"Fig_S15_Repeats_classes_codetection_cotargeting.pdf"

## Load relevant data
df_tar_det<-read.delim("cotargeting_vs_codetection_min100.tsv",stringsAsFactors = T)
## Transform the co-targeting column into a factor
df_tar_det$cotargeting <- factor(df_tar_det$cotargeting,ordered=T,levels=c("no","yes"))
summary(df_tar_det)

## Make a column with a pseudo-count of 0.1 for repeat pairs with 0 co-detections
df_tar_det$common_obs_pcount  <- df_tar_det$common_obs
df_tar_det[df_tar_det$common_obs==0,]$common_obs_pcount <- 0.1
summary(df_tar_det)

## Subsample pairs with 200,000 yes and 200,000 no to balance the dataset
df_tar_det_sub <- df_tar_det %>%
  group_by(cotargeting) %>%
  slice_sample(n=350000)

## Boxplot of the number of co-detection for co-targeting vs non-co-targeting
plot_1<-ggplot(df_tar_det_sub) + geom_boxplot(aes(x=cotargeting,y=common_obs_pcount,fill=cotargeting)) + scale_y_log10(breaks=c(0.1,1,100,10000),labels=c("0","1","100","10,000")) + xlab("Co-targeting detection") + ylab("Number of co-detection\n(i.e. number of samples in which both repeats were detected)") + blank_theme  + theme(panel.grid.major.y=element_line(color="grey",linewidth = 0.2)) + ggtitle("Number of co-detection for co-targeting and not co-targeting repeats") + scale_fill_brewer(palette="Pastel1") + theme(legend.position="bottom")
plot_1
ggsave(out_file,width=4,height=4)

### Stats on complete data and on subsample
wilcox.test(df_tar_det[df_tar_det$cotargeting=="no",]$common_obs,df_tar_det[df_tar_det$cotargeting=="yes",]$common_obs, correct=F)
wilcox.test(df_tar_det_sub[df_tar_det_sub$cotargeting=="no",]$common_obs,df_tar_det_sub[df_tar_det_sub$cotargeting=="yes",]$common_obs, correct=F)
