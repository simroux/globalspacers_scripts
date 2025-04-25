# Documentation/Archive on metadata collection for SRA runs used in the study
This is based on a dump of all "Run Info" from SRA for shotgun metagenomes, and an xml dump of the corresponding BioSample metadata, both read as tsv files. THIS IS FOR DOCUMENTATION PURPOSES ONLY, AND THERE ARE PROBABLY MORE EFFICIENT WAYS OF PERFORMING THE SAME TASKS NOW. These scripts are shared to show the logic of "semi-automated data curation" we used, but should probably not re-used as is.

## Step 1 - try to assign BioSamples to our ecosystem classification
`biosample_to_ecosystem.pl`
This script uses a dump of BioSample metadata ("All_biosample_info_NCBIAPI.tsv", transformed from an xml in a tsv), which looks like
```
Query	Attribute	Information
SAMN00272922	DATA_Organism	Homo sapiens
SAMN00272922	DATA_Title	RNA sample from blood of a human male participant in the dbGaP study "NHLBI Framingham SNP Health Association Resource (SHARe)"
SAMN00272922	DATA_analyte_type	RNA
SAMN00272922	DATA_biospecimen_repository	Framingham_SHARe
SAMN00272922	DATA_biospecimen_repository_sample_id	11334

```
along with a dump of SRA Run Information ("SraRunInfo_MetaG_not_Amplicon_Public-All-Combined-Nov_18_2023.csv") as a csv files, and some rules for translating keywords into our ecosystem classification ("Acceptable_combination.txt","Final_keyword_to_category.tsv","Ordered_category_scores.txt", provided here).

This script then tries to consider different information across the different BioSample attributes to assign each BioSample to an ecosystem type.

## Step 2 - link biosample information to SRA Runs
`get_srr_to_ecosystem.pl`
At this step, we filter some runs based on their sequencing characteristics (e.g. low-diversity metagenomes, amplified metagenomes, etc), and we also import additional ecosystem information for metagenomes found in IMG and MGNify. Eventually, we get a table linking individual runs to key run information (number of spots, number of bases) and key sample information (title, ecosystem type).


## Step 3 - prepare the database file
`prep_for_db.pl`
We put the final touches to our sample table, including (i) verifying that we have sample information for all spacers we will have in the database, (ii) excluding a few more samples (the ones with ecosystem "Exclude"), and (iii) adding the summarized ecosystem type in addition to the "full" ecosystem type

# Additional analyses
## Identify potential time series, i.e. sets of samples obtained from a single individual (for host-associated microbiomes) or the same location (for environmental microbiomes)
`identify_best_time_series.pl`
Here we identify the series of metagenomes that may represent time series, and we create the folder structure that will help us create separate files for each "time series".
