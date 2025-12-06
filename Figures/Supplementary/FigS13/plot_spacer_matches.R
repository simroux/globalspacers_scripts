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
# uvig_id <- "IMGVR_UViG_638275558_000001"
uvig_id <- "IMGVR_UViG_2911423748_000003"
### Choice of hosts
# host_one <- "d__Bacteria;p__Bacteroidota;c__Rhodothermia;o__Rhodothermales;f__Rhodothermaceae;g__Rhodothermus"
# host_two <- "d__Bacteria;p__Bacillota_B;c__Desulfotomaculia;o__Ammonifexales;f__Ammonificaceae;g__Ammonifex"
host_one <- "d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Agathobacter"
host_two <- "d__Bacteria;p__Bacillota_A;c__Clostridia;o__Lachnospirales;f__CAG-274;g__CAG-274"

###
df_known_ill <- read.delim(paste(uvig_id,"_spacer_hits_full.tsv",sep=""),stringsAsFactors = T)
# Get the genome map too
map_genome <- read.delim(paste(uvig_id,"_genome.tsv",sep=""))
map_genes <- read.delim(paste(uvig_id,"_genes.tsv",sep=""))
## Make a known / hypothetical protein column
map_genes <- map_genes %>%
  mutate(simp_product=case_when(sum_product == "hypothetical protein" ~ "hypothetical protein", .default="annotated"))
p1 <- gggenomes(seqs=map_genome,genes=map_genes) + geom_seq() + geom_gene(aes(fill = simp_product), position = "strand") + blank_theme + theme(legend.position="left")  + scale_fill_manual(values=c("#e8bc52","#595858"))+ scale_x_continuous(expand=expansion(mult = 0.01))
### Prep hits visualization
df_filt <- df_known_ill %>%
  filter(lca_origin=="Genome_high-confidence" | lca_origin=="Genome_medium-confidence") %>%
  filter(lca_genus=="d__Bacteria;p__Bacteroidota;c__Rhodothermia;o__Rhodothermales;f__Rhodothermaceae;g__Rhodothermus" | lca_genus=="d__Bacteria;p__Bacillota_B;c__Desulfotomaculia;o__Ammonifexales;f__Ammonificaceae;g__Ammonifex") %>%
  arrange(n_mismatches) %>% # To take the "best" hit for each spacer
  group_by(cluster_id) %>%
  summarise(n_mismatches=first(n_mismatches),hit_start=first(hit_start),hit_end=first(hit_end),lca_genus=first(lca_genus)) %>%
  mutate(avg_x=(hit_start+hit_end)/2) %>%
  mutate(width_x=(hit_end-hit_start+1)) %>%
  mutate(n_mis_f=factor(n_mismatches,ordered=T,levels=c("0","1","2","3")))
## Set the width at 200 to be visible
p2 <- ggplot(df_filt) + geom_tile(aes(x=avg_x,width=200,y=n_mismatches,fill=n_mis_f),linetype=0,alpha=0.8,height=0.5) + facet_grid(rows = vars(lca_genus)) + scale_fill_manual(values=c("#3c65e7","#4d9221","#fc8d59","#d73027")) + xlab("Coordinate on virus genome") + ylab("Number of mismatches") + theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="left",panel.grid.major=element_line(linewidth=0.3,colour="grey"),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.ticks.y=element_blank(),axis.text=element_text(color="black",size=8,family="Source Sans 3"),strip.text.y=element_blank(),legend.key.size = unit(1,"line"),strip.background=element_blank()) + ggtitle(i) + scale_x_continuous(expand=expansion(mult = 0.01))
p3 <- ggplot(df_filt) + geom_bar(aes(x=n_mismatches,fill=n_mis_f),col="black") + coord_flip() + facet_grid(rows = vars(lca_genus)) + scale_fill_manual(values=c("#3c65e7","#4d9221","#fc8d59","#d73027")) + xlab("") + ylab("Number of spacer hits") + theme(axis.line=element_blank(),panel.background=element_rect(linewidth=1,colour="black",fill="white"),legend.position="none",panel.grid.major=element_line(linewidth=0.3,colour="grey"),panel.grid.minor=element_blank(),text=element_text(color="black",size=8,family="Source Sans 3"),axis.ticks=element_line(color="black"),axis.ticks.y=element_blank(),axis.text=element_text(color="black",size=8,family="Source Sans 3"),axis.text.y=element_blank(),legend.key.size = unit(1,"line"),strip.background=element_blank(),strip.text.y = element_text(color="black",size=8,family="Source Sans 3",angle = 90))
# same width and height for plots when needed
gp1 <- ggplot_gtable(ggplot_build(p1))
gp2 <- ggplot_gtable(ggplot_build(p2))
gp3 <- ggplot_gtable(ggplot_build(p3))
gp1$widths<-gp2$widths
gp3$heights<-gp2$heights
# grid.arrange(gp2,gp3,gp1,ncol=2,nrow=2,widths=c(6,1))
pdf(paste(uvig_id,"_spacer_hits_summary.pdf",sep=""),width=9,height=4)
grid.arrange(gp2,gp3,gp1,ncol=2,nrow=2,widths=c(6,1))
dev.off()
