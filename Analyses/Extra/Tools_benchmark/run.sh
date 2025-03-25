## List_genomes_sub200.tsv has the list of genomes we selected
## Spacer_Table_sub200.tsv has the list of predicted spacers and repeats on these genomes, based on CRISPR-Cas Typer
## Run art to simulate reads
conda activate art
mkdir Reads_sub200/
./generate_reads.pl run
./clean_and_import_reads.pl run

### Run SpacerExtractor
conda activate spacerextractor_dev
./run_SE_base.pl run
## Run Crass base
conda activate crass
./run_crass_base.pl run
## Run MetaCRAST base
conda activate metaCRAST
./run_metacrast_base.pl run

### Compile results
./compare_all_tools.pl run
