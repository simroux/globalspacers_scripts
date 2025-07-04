# globalspacers_scripts
## Content
This repository includes the main scripts used in the global spacer project. **To access the data directly and perform new analyses, please see the [companion website](https://spacers.jgi.doe.gov)**. The scripts in this repository are instead shared for documentation purposes only, and are organized in two main folders. Files in "Analyses" include scripts generating files used at multiple steps and/or for multiple figures, for instance calculating spacer alpha diversity, matching spacers to IMG/VR and IMG/PR sequences, or building the duckdb database used to prepare some of the figures. Files in "Figures" include scripts used to generate the final input files as well as the R scripts that were used to generate the different visualizations.  

To retrace all the steps included, a folder called "Data" is included, and meant to host raw input files. Some of these raw input files are already included as part of this repository, in "Additional_data". The other raw input files are available as a tar archive named "Spacerdb_raw_files.tar.gz", which needs to be downloaded from ([https://portal.nersc.gov/dna/microbial/prokpubs/spacer_database_resources/](https://portal.nersc.gov/dna/microbial/prokpubs/spacer_database_resources/)), uncompressed, and moved to the "Data/Spacer_db" folder.  

## Usage
The data and scripts presented here are meant to provide documentation for reproducing the analysis and figures presented in the [global spacer manuscript (now on bioRxiv)](https://doi.org/10.1101/2025.06.12.659409). **Optimized versions of both scripts and databases are now available, so we encourage researchers to use these rather than the non-optimized data and scripts available here**. These new versions include the [SpacerExtractor tool](https://code.jgi.doe.gov/SRoux/spacerextractor) as well as [a companion website](https://spacers.jgi.doe.gov) illustrating how to access and interrogate the spacer database. 
Each figure and supplementary figure has its own folder in "Figures/", and each folder includes the relevant scripts and a Readme file. If the scripts rely on files previously generated as part of a broader analysis, this will be indicated in the Readme, and the corresponding scripts will be found in the "Analysis/" folders.  

## Dependencies
The scripts are a mix of Perl, Python, R, bash, and SQL queries executed in DuckDB. The main dependencies are:
* [DuckDB](https://duckdb.org). Version v0.10.1 4a89d97db8 was used in the analysis.
* Perl 5, with modules autodie, Cwd, File, Getopt, and strict (typically included by default), as well as Parallel::ForkManager and Text::Unidecode
* Python 3, with packages argparse, csv, os, and sys (typically included by default), as well as BioPython (SeqIO, Seq, and SeqUtils), duckdb, kcounter, pandas, and polars (all available from conda and/or PyPI)
* R 4, with libraries dplyr, extrafont, gggenomes, ggplot2, ggpubfigs, ggrastr, grid, gridExtra, RColorBrewer, stingr, sysfonts, and tidyr.
* For some analyses and/or figures, other external programs are used including minced, MetaCRAST, Crass, Bowtie1, blast, samtools, and SpacerExtractor.  

## Connected resources
* The tool [SpacerExtractor](https://code.jgi.doe.gov/SRoux/spacerextractor) allows users to mine their own metagenome (short reads) for CRISPR spacers, as well as matching (large) spacer databases to (large) databases of potential targets.
* A [companion website](https://spacers.jgi.doe.gov) provides access to a streamlined version of the spacer databases, and detailed examples for how to interrogate the database and extract e.g. information on spacers associated with a specific taxon or matching a specific target. **This is the recommended way to interact with the data for any new analysis or query**.

## Questions/Issues
* Found a bug, or have a question ? Please open an [issue](https://github.com/simroux/globalspacers_scripts/issues) in this repository ! 
