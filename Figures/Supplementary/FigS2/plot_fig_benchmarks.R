library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
loadfonts()
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(colour="grey"),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))

df_bench<-read.delim("Benchmark_method/Compare_spacer_recovery_sub200_metrics.tsv",stringsAsFactors = T)
df_bench$tool<-factor(df_bench$tool,ordered=T,levels=c("crass","metacrast","spacer_extractor","spacer_extractor_nosing"))
df_bench$coverage<-factor(df_bench$coverage,ordered=T,levels=c("0.8","10","100"))
df_bench$error_profile<-factor(df_bench$error_profile,ordered=T,levels=c("HS25qs3","HS25"))
summary(df_bench)

scale_tools<-c("crass"="#80B1D3","metacrast"="#8DD3C7","spacer_extractor"="#FDB462","spacer_extractor_nosing"="#FFFFB3")
p1 <- ggplot(df_bench,aes(x=coverage,y=recall_avg,fill=tool)) + geom_bar(stat="identity",position="dodge",col="black") + geom_errorbar(aes(x=coverage,ymin=recall_avg-recall_stderr,ymax=recall_avg+recall_stderr),width=0,position=position_dodge(.9)) + facet_grid(cols = vars(error_profile)) + scale_fill_manual(values=scale_tools) + scale_y_continuous(lim=c(0,1),labels=scales::percent) + blank_theme + xlab("Genome coverage (x)") + ylab("Recall")
p2 <- ggplot(df_bench,aes(x=coverage,y=fdr_avg,fill=tool)) + geom_bar(stat="identity",position="dodge",col="black") + geom_errorbar(aes(x=coverage,ymin=fdr_avg-fdr_stderr,ymax=fdr_avg+fdr_stderr),width=0,position=position_dodge(.9)) + facet_grid(cols = vars(error_profile)) + scale_fill_manual(values=scale_tools) + scale_y_continuous(lim=c(0,1),labels=scales::percent) + blank_theme + xlab("Genome coverage (x)") + ylab("FDR")
p3 <- ggplot(df_bench,aes(x=coverage,y=f1_avg,fill=tool)) + geom_bar(stat="identity",position="dodge",col="black") + geom_errorbar(aes(x=coverage,ymin=f1_avg-f1_stderr,ymax=f1_avg+f1_stderr),width=0,position=position_dodge(.9)) + facet_grid(cols = vars(error_profile)) + scale_fill_manual(values=scale_tools) + scale_y_continuous(lim=c(0,1),labels=scales::comma) + blank_theme + xlab("Genome coverage (x)") + ylab("F1-measure")

gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp2$widths<-gp1$widths
gp3$widths<-gp1$widths
pdf("Fig_S2_benchmarking.pdf",width=5,height=7)
grid.arrange(gp1,gp2,gp3,nrow=3)
dev.off()

