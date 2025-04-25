library(ggplot2)
library(ggpubfigs)
library(extrafont)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gggenomes)
source("../../color_scales.R")
custom_theme_hits <- theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_line(linewidth=0.3,colour="grey"),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.ticks.y=element_blank(),axis.text=element_text(color="black",size=8,family="Source Sans 3"),axis.text.y=element_blank(),legend.key.size = unit(1,"line"),strip.background=element_blank(),strip.text.y = element_text(color="black",size=8,family="Source Sans 3",angle = 0))

## Select which input file to work with (panel_c for Fig. S15, panel_d for Fig. S16]
panel <- "panel_c"
# panel <- "panel_d"

in_file <- paste(panel,"/coverage_for_fig.tsv",sep="")
in_df <- read.delim(in_file,stringsAsFactors = T)

## Prepare what's needed for the plots
in_df$avg_x<-(in_df$start+in_df$end)/2
in_df$width_x<-in_df$end-in_df$start+1
in_df$n_mis_f<-factor(in_df$n_mis,ordered=T,levels=c(0,1,2,3))
summary(in_df)

library(ggrastr)
list_viruses <- as.character(levels(in_df$virus))
## Prepare a plot for each virus, and print in a separate pdf
for (i in list_viruses){
  df_filt <- in_df %>%
    filter(virus==i)
  p1<-ggplot(df_filt) + rasterise(geom_tile(aes(x=avg_x,width=width_x,y=n_mis,fill=n_mis_f),linetype=0,alpha=0.8,height=0.5),dpi=600) + facet_grid(rows = vars(taxon)) + facet_grid(rows = vars(taxon)) + scale_fill_manual(values=c("#3c65e7","#4d9221","#fc8d59","#d73027")) + xlab("Coordinate on virus genome") + ylab("Number of mismatches") + theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="left",panel.grid.major=element_line(linewidth=0.3,colour="grey"),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.ticks.y=element_blank(),axis.text=element_text(color="black",size=8,family="Source Sans 3"),strip.text.y=element_blank(),legend.key.size = unit(1,"line"),strip.background=element_blank()) + ggtitle(i)
  p2<-ggplot(df_filt) + geom_bar(aes(x=n_mis,fill=n_mis_f),col="black") + coord_flip() + facet_grid(rows = vars(taxon)) + scale_fill_manual(values=c("#3c65e7","#4d9221","#fc8d59","#d73027")) + xlab("") + ylab("Number of spacer hits") + theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="none",panel.grid.major=element_line(linewidth=0.3,colour="grey"),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.ticks.y=element_blank(),axis.text=element_text(color="black",size=8,family="Source Sans 3"),axis.text.y=element_blank(),legend.key.size = unit(1,"line"),strip.background=element_blank(),strip.text.y = element_text(color="black",size=8,family="Source Sans 3",angle = 0))
  ## Plot both atop of one another in a single pdf
  gp1 <- ggplot_gtable(ggplot_build(p1))
  gp2 <- ggplot_gtable(ggplot_build(p2))
  gp1$heights<-gp2$heights
  grid.arrange(gp1,gp2,ncol=2,widths=c(3,1.3))
  n_hosts <- length(unique(df_filt$taxon))
  height <- 0.65*n_hosts
  out_file <- paste(panel,"/Split_figs_",i,"_panel.pdf",sep="")
  pdf(out_file,width=9.5,height=height)
  grid.arrange(gp1,gp2,ncol=2,widths=c(3,1.3))
  dev.off()
}
