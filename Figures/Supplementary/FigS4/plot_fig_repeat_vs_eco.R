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
df_sample_detect <-read.delim("Figures/Main/Fig_1/sample_to_array_counts.tsv",stringsAsFactors = T)
df_sample_detect$sra_run <- df_sample_detect$sra_run_1
df_sample_detect <- df_sample_detect %>%  
  mutate(total_cluster = replace_na(total_cluster, 0))
df_sample_detect$ecosystem_sum <- factor(df_sample_detect$ecosystem_sum,ordered=T,levels=levels_ecosystem)

### First- Illustration of the relationship between sequencing and number of arrays detected => x-y plot
## Split into 3 sub-panels: animal-associated, aquatic, and other ecosystems, so that it's more readable
df_animal <- df_sample_detect %>%
  filter(ecosystem_sum=="Human-associated_Digestive-system" | ecosystem_sum=="Human-associated_Other" | ecosystem_sum=="Animal-associated")  
## Need to down-sample to 5k
df_animal_fordots <- df_animal %>%
  group_by(ecosystem_sum) %>%
  sample_n(1000) %>%
  ungroup()
p1<-ggplot() + geom_point(data=df_animal_fordots,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),pch=19,size=1.5,alpha=0.3) + geom_smooth(data=df_animal,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),method="lm",formula=y~x+0)+ geom_rect(aes(xmin=0,xmax=2E+10,ymin=0,ymax=300)) + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3)) + scale_fill_manual(values=scale_ecosystem) + scale_color_manual(values=scale_ecosystem) + scale_y_continuous(lim=c(0,3100),labels=scales::comma) + xlab("Number of bases in metagenome") + ylab("Number of distinct CRISPR arrays in final dataset") + theme(legend.position="none") 
p1_inset<-ggplot() + geom_point(data=df_animal_fordots,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),pch=19,size=1.5,alpha=0.7) + geom_smooth(data=df_animal,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),method="lm",formula=y~x+0) + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3)) + scale_fill_manual(values=scale_ecosystem) + scale_color_manual(values=scale_ecosystem) + scale_y_continuous(lim=c(0,3100),labels=scales::comma) + xlab("Number of bases in metagenome") + ylab("Number of distinct CRISPR arrays in final dataset") + theme(legend.position="none") + coord_cartesian(xlim=c(0,2E+10),ylim=c(0,300))
df_aquatic <- df_sample_detect %>%
  filter(ecosystem_sum=="Aquatic_Marine" | ecosystem_sum=="Aquatic_Thermal-springs" | ecosystem_sum=="Aquatic_Sediment" | ecosystem_sum=="Aquatic_Freshwater" | ecosystem_sum=="Aquatic_Other")
df_aquatic_fordots <- df_aquatic %>%
  group_by(ecosystem_sum) %>%
  sample_n(1000) %>%
  ungroup()
p2<-ggplot() + geom_point(data=df_aquatic_fordots,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),pch=19,size=1.5,alpha=0.3) + geom_smooth(data=df_aquatic,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),method="lm",formula=y~x+0) + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3)) + scale_fill_manual(values=scale_ecosystem) + scale_color_manual(values=scale_ecosystem) + scale_y_continuous(lim=c(0,3100),labels=scales::comma) + xlab("Number of bases in metagenome") + ylab("Number of distinct CRISPR arrays in final dataset") + theme(legend.position="none")
p2_inset<-ggplot() + geom_point(data=df_aquatic_fordots,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),pch=19,size=1.5,alpha=0.7) + geom_smooth(data=df_aquatic,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),method="lm",formula=y~x+0) + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3)) + scale_fill_manual(values=scale_ecosystem) + scale_color_manual(values=scale_ecosystem) + scale_y_continuous(lim=c(0,3100),labels=scales::comma) + xlab("Number of bases in metagenome") + ylab("Number of distinct CRISPR arrays in final dataset") + theme(legend.position="none") + coord_cartesian(xlim=c(0,2E+10),ylim=c(0,300))

df_other <- df_sample_detect %>%
  filter(ecosystem_sum=="Terrestrial_Soil_and_plants" | ecosystem_sum=="Engineered" | ecosystem_sum=="Other")
df_other_fordots <- df_other %>%
  group_by(ecosystem_sum) %>%
  sample_n(1000) %>%
  ungroup()
p3<-ggplot() + geom_point(data=df_other_fordots,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),pch=19,size=1.5,alpha=0.3) + geom_smooth(data=df_other,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),method="lm",formula=y~x+0)  + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3)) + scale_fill_manual(values=scale_ecosystem) + scale_color_manual(values=scale_ecosystem) + scale_y_continuous(lim=c(0,3100),labels=scales::comma) + xlab("Number of bases in metagenome") + ylab("Number of distinct CRISPR arrays in final dataset") + theme(legend.position="none")
p3_inset<-ggplot() + geom_point(data=df_other_fordots,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),pch=19,size=1.5,alpha=0.7) + geom_smooth(data=df_other,aes(x=n_base,y=total_cluster,fill=ecosystem_sum,col=ecosystem_sum),method="lm",formula=y~x+0)  + blank_theme + theme(panel.grid.major=element_line(colour="lightgrey",linewidth=0.3)) + scale_fill_manual(values=scale_ecosystem) + scale_color_manual(values=scale_ecosystem) + scale_y_continuous(lim=c(0,3100),labels=scales::comma) + xlab("Number of bases in metagenome") + ylab("Number of distinct CRISPR arrays in final dataset") + theme(legend.position="none") + coord_cartesian(xlim=c(0,2E+10),ylim=c(0,300))
## One for the legend
pl <- ggplot(df_sample_detect,aes(x=n_base,y=total_cluster,col=ecosystem_sum)) + geom_smooth(method="lm",formula=y~x+0) + blank_theme + scale_color_manual(values=scale_ecosystem)
# Prep the grid
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gpl <- ggplot_gtable(ggplot_build(pl))
gp1i <- ggplot_gtable(ggplot_build(p1_inset))
gp2i <- ggplot_gtable(ggplot_build(p2_inset))
gp3i <- ggplot_gtable(ggplot_build(p3_inset))
##
gp2$widths <- gp1$widths
gp3$widths <- gp1$widths
gp2$heights <- gp1$heights
gp3$heights <- gp1$heights
##
gp2i$widths <- gp1i$widths
gp3i$widths <- gp1i$widths
gp2i$heights <- gp1i$heights
gp3i$heights <- gp1i$heights
##
grid.arrange(gp1,gp2,gp3,gpl,ncol=3,nrow=2)
pdf("Fig_S4_top_panel.pdf",width=9,height=6)
grid.arrange(gp1,gp2,gp3,gpl,ncol=3,nrow=2)
dev.off()

pdf("Fig_S4_insets.pdf",width=9,height=6)
grid.arrange(gp1i,gp2i,gp3i,gpl,ncol=3,nrow=2)
dev.off()



### Then - Two bar charts, with detailed ecosystem:
## one with the total number of samples in this ecosystem type, the other with the estimated average number of CRISPR per 1Gb
mean_w_err_eco2<-data.frame(eco=character(),eco_3=character(),total_obs=numeric(),total_mean=numeric(),total_lower=numeric(),total_upper=numeric(),rsq=numeric(),pval=numeric())
new.dat<-data.frame(n_base=1E+09) ## This will be the number for which we will predict this (i.e. 1Gb)
for (eco in levels(df_sample_detect$ecosystem_base)){
  print(eco)
  test<-df_sample_detect[df_sample_detect$ecosystem_base==eco,]
  total_obs <- nrow(test)
  eco_3 <- unique(test$ecosystem_sum)
  if (total_obs >= 50){
    lm.model<-lm(total_cluster ~ n_base + 0, data=test) ## Linear model, forcing it through zero
    toto <- summary(lm.model)
    titi <- predict(lm.model, newdata=new.dat, interval='confidence',level=0.99) ## This will give us a predicted value for 1E+09 with a confidence interval for this predicted value. We use "confidence" here, i.e. we have the confidence interval of the mean, not for any individual predicted value (see https://rpubs.com/aaronsc32/regression-confidence-prediction-intervals)
    mean_w_err_eco2<-rbind(mean_w_err_eco2,data.frame(eco,eco_3,total_obs,titi[1],titi[2],titi[3],toto$r.squared,toto$coefficients[4]))
  }
}
colnames(mean_w_err_eco2) <- c("eco","eco_color","total","mean","lower","upper","rsquared","pvalue")
mean_w_err_eco2 %>%
  arrange(mean)
export_tab <- df_sample_detect %>%
  group_by(ecosystem_sum,ecosystem_base) %>%
  summarise(n_sample = n()) %>%
  left_join(mean_w_err_eco2, by = join_by(ecosystem_base==eco)) %>%
  arrange(desc(mean)) %>%
  print(n=100) 
write.csv(export_tab,file="Supp_mat/Sup_table_number_of_array_per_sample_ecodetailed.csv",quote=F)

mean_w_err_eco2<-mean_w_err_eco2 %>%
  arrange(eco_color,eco) %>% 
  mutate(eco = factor(eco, ordered=T, levels=unique(eco)))
plot_1<-ggplot(mean_w_err_eco2,aes(x=eco,y=total)) + geom_bar(fill="grey",stat="identity",alpha=0.8,width=0.75,col="black") + scale_y_continuous(labels=scales::comma) + blank_theme + xlab("Ecosystem type") + ylab("Total number of samples") + coord_flip()
plot_2<-ggplot(mean_w_err_eco2,aes(x=eco,y=mean)) + geom_bar(aes(fill=eco_color),stat="identity",alpha=0.8,width=0.75,col="black") + geom_linerange(aes(ymin=lower, ymax=upper), col="gray30", alpha=0.8) + blank_theme  + scale_fill_manual(values=scale_ecosystem) + xlab("Ecosystem type") + ylab("Average number of CRISPR array (per 1Gb)") + coord_flip()

gp1 <- ggplot_gtable(ggplot_build(plot_1))
gp2 <- ggplot_gtable(ggplot_build(plot_2))
## One way
gp1$widths <- gp2$widths
gp1$heights <- gp2$heights
##
grid.arrange(gp1,gp2,ncol=2)
pdf("Fig_S4_bottom_panel.pdf",width=15,height=6)
grid.arrange(gp1,gp2,ncol=2)
dev.off()

