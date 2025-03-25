## Source of taxo info
scale_taxo<-c("Genome_high-confidence"="#00441b","Genome_medium-confidence"="#41ab5d","Genome_low-confidence"="#c7e9c0","Contig"="#fed976","No-information_long"="#737373","No-information_short"="#f0f0f0")
levels_taxo<-rev(c("Genome_high-confidence","Genome_medium-confidence","Genome_low-confidence","Contig","No-information_long","No-information_short"))

## CRISPR Types
# scale_type<-c("#f0f0f0",brewer.pal("Set3",n=6))
scale_type<-c("Unknown"="#f0f0f0","VI"="#8DD3C7","V"="#FFFFB3","IV"="#BEBADA","III"="#FB8072","II"="#80B1D3","I"="#FDB462")
levels_type<-c("Unknown","VI","V","IV","III","II","I")

## Phyla
scale_phylum<-c("d__Viruses_p__Uroviricota"="#fa9fb5","d__Archaea;Other"="#ffffe5","d__Archaea_p__Thermoproteota"="#fec44f","d__Archaea_p__Halobacteriota"="#fd8d3c","d__Bacteria_p__unclassified"="#f0f0f0","d__Bacteria;Other"="#737373","d__Bacteria_p__Verrucomicrobiota"="#bcbddc","d__Bacteria_p__Pseudomonadota"="#081d58","d__Bacteria_p__Planctomycetota"="#d9f0a3","d__Bacteria_p__Desulfobacterota"="#41b6c4","d__Bacteria_p__Cyanobacteriota"="#006837","d__Bacteria_p__Chloroflexota"="#f7fcb9","d__Bacteria_p__Campylobacterota"="#fc9272","d__Bacteria_p__Bacteroidota"="#d7301f","d__Bacteria_p__Bacillota_C"="#67000d","d__Bacteria_p__Bacillota_A"="#ef6548","d__Bacteria_p__Bacillota"="#fdbb84","d__Bacteria_p__Actinomycetota"="#662506","d__Bacteria_p__Acidobacteriota"="#cc4c02")
levels_phylum<-c("d__Viruses_p__Uroviricota","d__Archaea;Other","d__Archaea_p__Thermoproteota","d__Archaea_p__Halobacteriota","d__Bacteria_p__unclassified","d__Bacteria;Other","d__Bacteria_p__Verrucomicrobiota","d__Bacteria_p__Pseudomonadota","d__Bacteria_p__Planctomycetota","d__Bacteria_p__Desulfobacterota","d__Bacteria_p__Cyanobacteriota","d__Bacteria_p__Chloroflexota","d__Bacteria_p__Campylobacterota","d__Bacteria_p__Bacteroidota","d__Bacteria_p__Bacillota_C","d__Bacteria_p__Bacillota_A","d__Bacteria_p__Bacillota","d__Bacteria_p__Actinomycetota","d__Bacteria_p__Acidobacteriota")

## Ecosystems
scale_ecosystem <- c("Human-associated_Digestive-system"="#ef3b2c","Human-associated_Other"="#f768a1","Animal-associated"="#a24963","Engineered"="#9e9ac8","Aquatic_Freshwater"="#7bccc4","Aquatic_Marine"="#2171b5","Aquatic_Thermal-springs"="#fec44f","Aquatic_Sediment"="#9bc14a","Aquatic_Other"="#c6dbef","Terrestrial_Soil_and_plants"="#41ab5d","Other"="#969696","Unknown"="#d9d9d9")
levels_ecosystem<-c("Human-associated_Digestive-system","Human-associated_Other","Animal-associated","Engineered","Aquatic_Freshwater","Aquatic_Marine","Aquatic_Thermal-springs","Aquatic_Sediment","Aquatic_Other","Terrestrial_Soil_and_plants","Other","Unknown")

## Taxo levels
scale_taxolvl<-c("none"="#bdbdbd","domain"="#c7e9b4","phylum"="#7fcdbb","class"="#41b6c4","order"="#1d91c0","family"="#225ea8","genus"="#253494","species"="#081d58")
levels_taxolvl<-c("none","domain","phylum","class","order","family","genus","species")

## Size of sets
levels_sizeset<-c("[1]","[5-19]","[20-99]","[100p]") ## [1] should never ever happen in theory
scale_sizeset<-c("[1]"="#edf8b1","[5-19]"="#7fcdbb","[20-99]"="#1d91c0","[100p]"="#0c2c84")

## Max coverage of set 
level_sizecover<-c("[1]","[2-4]","[5-9]","[10-19]","[20-49]","[50-99]","[100p]")
scale_sizecover<-c("[1]"="#fee8c8","[2-4]"="#fdd49e","[5-9]"="#fdbb84","[10-19]"="#fc8d59","[20-49]"="#ef6548","[50-99]"="#d7301f","[100p]"="#b30000")

## Type of spacer
level_spacertypes<-rev(c("all","common","rare"))
scale_spacertypes<-c("rare"="#fec44f","common"="#3690c0","all"="#969696")

## Taxonomy ranks
scale_taxorank<-c("domain"="#e9a3c9","phylum"="#d73027","class"="#fc8d59","order"="#fee090","family"="#4d9221","genus"="#91bfdb","species"="#4575b4")
