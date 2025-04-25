# How-to - reconstruct panels of Fig. S13
This supplementary figure requires a few external files, but can be (relatively) easily adjusted to visualize the targeting for any virus-repeat pair.  
First, gff files need to be extracted from IMG. For the two viruses used in Fig. S13, the gff files are available in "Data/Additional/Data". Next, you can generate the input file needed for gggenomes using `get_input_gggenomes.pl`.  

The second step is to prepare the file for the spacer hits. First, use `get_single_virus_table.py` to extract all spacers matching the virus of interest. Then, `filter_hit_table.pl` can be used to crate the input file for the plot for a given repeat, and the plots can be generated with `plot_spacer_matches.R`.  

For Fig. S13, the following commands were used:
```
./get_input_gggenomes.pl -i IMGVR_UViG_3300007356_000002
python get_single_virus_table.py -i IMGVR_UViG_3300007356_000002
./filter_hit_table.pl -i IMGVR_UViG_3300007356_000002_spacer_hits_full.tsv -a Ac_12820
```
and
```
./get_input_gggenomes.pl -i IMGVR_UViG_3300045988_178767 
python get_single_virus_table.py -i IMGVR_UViG_3300045988_178767
./filter_hit_table.pl -i IMGVR_UViG_3300045988_178767_spacer_hits_full.tsv -a Ac_16888
```
then use `plot_spacer_matches.R` to prepare the plot panels.