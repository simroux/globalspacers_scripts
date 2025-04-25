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
blank_theme<-theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="right",panel.grid.major=element_blank(),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.text=element_text(color="black",size=8,family="Source Sans 3"),legend.key.size = unit(1,"line"))

### Choice of UViG
# uvig_id <- "IMGVR_UViG_3300045988_178767"
uvig_id <- "IMGVR_UViG_3300007356_000002"
###

df_exone<-read.delim(paste(uvig_id,"_spacer_hits_selected.tsv",sep=""),stringsAsFactors = T)
df_exone$avg_x<-(df_exone$hit_start+df_exone$hit_end)/2
df_exone$width_x<-df_exone$hit_end-df_exone$hit_start+1
df_exone$is_unique<-"no"
test <- nrow(df_exone[df_exone$n_votu==1,])
if (test>0){
  df_exone[df_exone$n_votu==1,]$is_unique<-"yes"  
}
df_exone$is_unique<-factor(df_exone$is_unique)
# Select libraries
list_lib <- df_exone %>% 
  group_by(library) %>%
  summarise(individual_sample_cover=first(individual_sample_cover)) %>%
  arrange(individual_sample_cover) %>%
  tail(n=10) ## Show only top 10 libraries, should not do anything if there are more than 10
list_lib
df_exone$library<-factor(df_exone$library,ordered=T,levels=list_lib$library)
summary(df_exone)
# Get the genome map too
map_genome <- read.delim(paste(uvig_id,"_genome.tsv",sep=""))
map_genes <- read.delim(paste(uvig_id,"_genes.tsv",sep=""))
## Make a known / hypothetical protein column
map_genes <- map_genes %>%
  mutate(simp_product=case_when(sum_product == "hypothetical protein" ~ "hypothetical protein", .default="annotated"))
# And also get the histogram of the array size
## while determining targeting level too
df_tmp <- df_exone %>% 
  group_by(library) %>%
  filter(!is.na(library)) %>%
  summarise(hits_target=max(hits_in_sample),hits_offtarget=max((n_spacers-hits_in_sample)),n_hits=n(),individual_sample_cover=first(individual_sample_cover)) %>%
  mutate(targeting_level = case_when(n_hits < 10 | individual_sample_cover < 200 ~ "low", individual_sample_cover >= 1000 ~ "high", .default="medium")) %>%
  pivot_longer(names_to="hit_type",values_to="count",cols=c("hits_target","hits_offtarget")) %>%
  mutate(library_w_level = paste(library,targeting_level,sep=" - "))
df_tmp$library<-factor(df_tmp$library,ordered=T,levels=list_lib$library)
summary(df_tmp)

### Figure
p1 <- gggenomes(seqs=map_genome,genes=map_genes) + geom_seq() + geom_gene(aes(fill = simp_product), position = "strand") + blank_theme + theme(legend.position="left")  + scale_fill_manual(values=c("#595858","#e5e5e5"))
p2 <- ggplot(df_exone) + geom_tile(aes(x=avg_x,width=width_x,y="1",fill=is_common,alpha=is_common),height=0.5) + scale_fill_manual(values=c("#2c7fb8","#fc4e2a")) + scale_alpha_manual(values=c(0.2,1)) + scale_x_continuous(lim=c(0,map_genome$length[1])) + blank_theme + theme(legend.position="left") + xlab("Virus genome") + ylab("Samples")
p3 <- ggplot(df_exone[!is.na(df_exone$library),]) + geom_tile(aes(x=avg_x,width=width_x,y=library,fill=is_common,col=is_common),height=0.5,alpha=0.5) + scale_fill_manual(values=c("#2c7fb8","#fc4e2a")) + scale_color_manual(values=c("#2c7fb8","#fc4e2a")) + scale_x_continuous(lim=c(0,map_genome$length[1])) + blank_theme + theme(legend.position="left") + xlab("Virus genome") + ylab("Samples")
p4<-ggplot(df_tmp) + geom_bar(aes(x=library,y=count,fill=hit_type),stat="identity",col="black",width=0.7) + scale_x_discrete(position="top",breaks=df_tmp$library,labels=df_tmp$library_w_level) + coord_flip() + blank_theme + scale_fill_manual(values=c("white","black")) + theme(legend.position="bottom")
# same width for plots
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp4 <- ggplot_gtable(ggplot_build(p4))

## Alternative way:
gp2$widths<-gp3$widths
gp1$widths<-gp3$widths
gp4$heights<-gp3$heights
blank <- grid.rect(gp=gpar(col="white"))
## Export
pdf(paste(uvig_id,"_spacer_hits_selected.pdf",sep=""),width=9,height=4)
grid.arrange(gp1, blank, gp2, blank, gp3, gp4, nrow=3, ncol=2,heights=c(2,1,3),widths=c(5,1.5))
dev.off()
