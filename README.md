# oral_comparison
Data and R scripts for comparison of oral microbiome datasets

<b>Data files/</b>
* `metadata_all_samples.tsv` : Metadata for all samples used in the analysis.
* `tax_lineage_all_samples.tsv` : Taxonomic analysis results for all samples.
* `subgingival_conditions.txt` : Disease/health metadata for the subgingival plaque samples.
* `tax_lineage_subgingival.tsv` : Taxonomic analysis resilts for subgingival plaque samples.
* `subg_ml_with_metadata.csv` : Taxonomic analysis results (clr) and metadata for subgingival plaque samples for machine learning task.

<b>R/</b>
* `all_plaque_pca.R` : Converts taxonomic analysis results from raw counts to clr and performcs PCA and ANOSIM.
* `effect_plot.R` : Converts taxonomic analysis results from raw counts to clr and produces effect plot.
* `machine_learning.R` : Trains and tests machine learning classifiers.
