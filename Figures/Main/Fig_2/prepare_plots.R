library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
source("color_scales.R")
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))
## Panel A - Global distribution of spacer set size
## Load data
df_setsize<-read.delim("Spacer_set_size_info.tsv",stringsAsFactors = T)
df_setsize$ecosystem <- factor(df_setsize$ecosystem,ordered=T,levels=levels_ecosystem)
df_setsize$type <- factor(df_setsize$type,ordered=T,levels=levels_type)
## Plots
p1 <- ggplot(df_setsize) + geom_histogram(aes(x=spacer_set_size),bins=60) + xlab("spacer set size (number of spacers)") + ylab("number of observations") + blank_theme + scale_x_log10(labels=scales::comma,breaks=c(5,10,20,100,1000,10000,100000),lim=c(min(df_setsize$spacer_set_size),max(df_setsize$spacer_set_size)),expand=expansion(add=0.05)) + scale_y_continuous(labels=scales::comma,expand=expansion(mult=0.02))
# Make categories for maximum set coverage and set size, then show the boxplot
df_setsize$cat_maxcover <-  cut(df_setsize$max_cover,c(0,10,50,100,1000,1E100))
df_setsize$cat_maxcover <- df_setsize$cat_maxcover %>%
  recode("(0,10]" = "0-10", "(10,50]" = "11-50", "(50,100]" = "51-100", "(100,1e+03]" = "101-1000", "(1e+03,1e+100]" = "> 1000")
# Add boxplot for how maximum set coverage is linked to set size
p2 <- ggplot(df_setsize) + geom_boxplot(aes(x=cat_maxcover,y=spacer_set_size,fill=cat_maxcover),outliers=FALSE) + xlab("maximum spacer coverage in set") + ylab("spacer set size (number of spacers)") + blank_theme + scale_y_log10(labels=scales::comma,breaks=c(5,10,20,100,1000,10000,100000),lim=c(min(df_setsize$spacer_set_size),max(df_setsize$spacer_set_size)),expand=expansion(add=0.05)) + coord_flip() + scale_fill_brewer(palette="Reds") + theme(legend.position="none")
## Plot both atop of one another in a single pdf
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp1$widths<-gp2$widths
pdf("Fig2A.pdf",width=3,height=2.5)
grid.arrange(gp1,gp2,nrow=2,heights=c(2,1.3))
dev.off()

### Panel B: typical spacer coverage curve
df_adiv_curve <- read.delim("spacer_sets_alphadiv_forplot_subsamples.tsv", stringsAsFactors = T)
df_adiv_curve$ecosystem<-factor(df_adiv_curve$ecosystem,ordered=T,levels=levels_ecosystem)
df_adiv_curve$code<-paste(df_adiv_curve$crispr_array,df_adiv_curve$sra_run,sep="_")
summary(df_adiv_curve)
## List cases for which we have a max >=20x
relevant_codes <- df_adiv_curve %>%
  group_by(code) %>%
  summarise(max_coverage=max(total_coverage)) %>%
  filter(max_coverage>=20) 
df_adiv_curve_filt <- df_adiv_curve %>%
  filter(code %in% relevant_codes$code)
df_adiv_curve_summary <- df_adiv_curve_filt %>%
  group_by(rank) %>%
  summarise(mean_coverage = mean(total_coverage), mean_relative_coverage = mean(relative_coverage), median_coverage = median(total_coverage), median_relative_coverage = median(relative_coverage))
df_adiv_curve_points <- df_adiv_curve_filt %>%
  filter(total_coverage > 0) %>%
  slice_sample(n = 10000) %>%
  mutate(category = cut(relative_coverage,c(0,0.05,0.5,2),include.lowest=TRUE,right=FALSE))
ggplot() + geom_point(data=df_adiv_curve_points,aes(x=rank,y=relative_coverage,col=category),alpha=0.15) + geom_smooth(data=df_adiv_curve_summary,aes(x=rank,y=median_relative_coverage),method="loess",span=0.05,alpha=1,col="black") + scale_x_log10(labels=scales::comma,expand=expansion(mult = 0.01)) + scale_y_continuous(expand=expansion(mult = 0.01),labels=scales::percent,breaks=c(0,0.05,0.25,0.5,0.75,1)) + geom_hline(yintercept=0.05,linetype="dashed",col="blue") + geom_hline(yintercept=0.5001,linetype="dashed",col="red") + blank_theme  + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "inches"),legend.position="none") + ylab("Relative coverage within spacer set") + xlab("Rank within spacer set") + scale_color_manual(values=c("blue","gray","red"))
ggsave("Fig_2B_alphadiv_summary.pdf",width=5,height=3.5)


### Panel C - singletons percentage at different clustering level
sing_col <- c("#1c9099","#feb24c","#e34a33")
df_singleton<-read.delim("input_figure_singletons.tsv",stringsAsFactors = T)
## For simplicity, we only show up to 200 spacers, we will have the boxplots in supplementary figures
## we also calculate smoothed lines for the median, as well as max and min of the ribbon
df_singplot <- df_singleton %>%
  filter(total_spacer<=200) %>%
  group_by(id) %>%
  mutate(smoothed_median = (loess(median ~ total_spacer, span=0.25))$fitted) %>%
  mutate(smoothed_min = (loess(stderr_min ~ total_spacer, span=0.25))$fitted) %>%
  mutate(smoothed_max = (loess(stderr_max ~ total_spacer, span=0.25))$fitted) %>%
  ungroup()
ggplot(df_singplot) + geom_ribbon(aes(x=total_spacer,ymin=smoothed_min,ymax=smoothed_max,fill=id),alpha=1) + geom_line(aes(x=total_spacer,y=smoothed_median,col=id),linewidth=0.5) + scale_x_continuous(lim=c(5,200),breaks=c(5,50,100,150,200),expand = expansion(add=0)) + scale_y_continuous(lim=c(0,1),labels=scales::percent) + xlab("Spacer set size (number of spacers)") + ylab("Percentage of singletons") + blank_theme + theme(panel.grid.major.y = element_line(colour="gray")) + scale_colour_manual(values=sing_col) + scale_fill_manual(values=sing_col) + theme(legend.position="bottom")
ggsave("Fig_2C_singletons.pdf",width=3,height=2.25)

################## Panel D
df_curvebeta <- read.delim("spacer_sets_betadiv_forplot_subsamples-summarized.tsv",stringsAsFactors = T)
df_curvebeta$ecosystem<-factor(df_curvebeta$ecosystem,ordered=T,levels=c(levels_ecosystem,"all"))
df_curvebeta_all <- df_curvebeta %>%
  filter(ecosystem=="all")
df_curvebeta_all_medians <- df_curvebeta_all %>%
  group_by(n_samples) %>%
  summarise(median_percentage_all=median(percentage_all),median_percentage_sgton=median(percentage_sgton),median_percentage_common=median(percentage_common))
p1 <- ggplot() + geom_point(data=df_curvebeta_all,aes(x=n_samples,y=percentage_all),col="gray",alpha=0.15) + geom_line(data=df_curvebeta_all_medians,aes(x=n_samples,y=median_percentage_all),col="black",linewidth=1.5) + scale_x_log10(labels=scales::comma,expand=expansion(mult = 0.03)) + scale_y_continuous(expand=expansion(mult = 0.03),labels=scales::percent,lim=c(0,1)) + geom_vline(xintercept=1.5,linetype="dashed",col="blue") + blank_theme  + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "inches")) + ylab("Percentage of spacer in set") + xlab("Number of samples")
p1
p2 <- ggplot() + geom_point(data=df_curvebeta_all,aes(x=n_samples,y=percentage_common),col="gray",alpha=0.15) + geom_line(data=df_curvebeta_all_medians,aes(x=n_samples,y=median_percentage_common),col="black",linewidth=1.5) + scale_x_log10(labels=scales::comma,expand=expansion(mult = 0.03)) + scale_y_continuous(expand=expansion(mult = 0.03),labels=scales::percent,lim=c(0,1)) + geom_vline(xintercept=1.5,linetype="dashed",col="blue") + blank_theme  + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "inches")) + ylab("Percentage of spacer in set") + xlab("Number of samples")
p2
p3 <- ggplot() + geom_point(data=df_curvebeta_all,aes(x=n_samples,y=percentage_sgton),col="gray",alpha=0.15) + geom_line(data=df_curvebeta_all_medians,aes(x=n_samples,y=median_percentage_sgton),col="black",linewidth=1.5) + scale_x_log10(labels=scales::comma,expand=expansion(mult = 0.03)) + scale_y_continuous(expand=expansion(mult = 0.03),labels=scales::percent,lim=c(0,1)) + geom_vline(xintercept=1.5,linetype="dashed",col="blue") + blank_theme  + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "inches")) + ylab("Percentage of spacer in set") + xlab("Number of samples")
p3

### Same panel, but within subject
df_betacurve_ind <- read.delim("spacer_sets_betadiv_forplot_within-individuals.tsv",stringsAsFactors = T)
df_betacurve_ind_medians <- df_betacurve_ind %>%
  group_by(n_samples) %>%
  summarise(median_percentage_all=median(percentage_all),median_percentage_sgton=median(percentage_sgton),median_percentage_common=median(percentage_common))
pi1 <- ggplot() + geom_point(data=df_betacurve_ind,aes(x=n_samples,y=percentage_all),col="gray",alpha=0.15) + geom_line(data=df_betacurve_ind_medians,aes(x=n_samples,y=median_percentage_all),col="black",linewidth=1.5) + scale_x_log10(labels=scales::comma,expand=expansion(mult = 0.03)) + scale_y_continuous(expand=expansion(mult = 0.03),labels=scales::percent,lim=c(0,1)) + geom_vline(xintercept=1.5,linetype="dashed",col="blue") + blank_theme  + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "inches")) + ylab("Percentage of spacer in set") + xlab("Number of samples")
pi1
pi2 <- ggplot() + geom_point(data=df_betacurve_ind,aes(x=n_samples,y=percentage_common),col="gray",alpha=0.15) + geom_line(data=df_betacurve_ind_medians,aes(x=n_samples,y=median_percentage_common),col="black",linewidth=1.5) + scale_x_log10(labels=scales::comma,expand=expansion(mult = 0.03)) + scale_y_continuous(expand=expansion(mult = 0.03),labels=scales::percent,lim=c(0,1)) + geom_vline(xintercept=1.5,linetype="dashed",col="blue") + blank_theme  + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "inches")) + ylab("Percentage of spacer in set") + xlab("Number of samples")
pi2
pi3 <- ggplot() + geom_point(data=df_betacurve_ind,aes(x=n_samples,y=percentage_sgton),col="gray",alpha=0.15) + geom_line(data=df_betacurve_ind_medians,aes(x=n_samples,y=median_percentage_sgton),col="black",linewidth=1.5) + scale_x_log10(labels=scales::comma,expand=expansion(mult = 0.03)) + scale_y_continuous(expand=expansion(mult = 0.03),labels=scales::percent,lim=c(0,1)) + geom_vline(xintercept=1.5,linetype="dashed",col="blue") + blank_theme  + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "inches")) + ylab("Percentage of spacer in set") + xlab("Number of samples")
pi3

## Plot both atop of one another in a single pdf
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp4 <- ggplot_gtable(ggplot_build(pi1))
gp5 <- ggplot_gtable(ggplot_build(pi2))
gp6 <- ggplot_gtable(ggplot_build(pi3))
gp2$widths<-gp1$widths
gp3$widths<-gp1$widths
gp4$widths<-gp1$widths
gp5$widths<-gp1$widths
gp6$widths<-gp1$widths
pdf("fig2D.pdf",width=9,height=5)
grid.arrange(gp1,gp2,gp3,gp4,gp5,gp6,ncol=3,nrow=2)
dev.off()
