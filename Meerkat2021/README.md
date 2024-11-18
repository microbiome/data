
## Diurnal Oscillations in Gut Bacterial Load and Composition

This dataset, titled **"Diurnal oscillations in gut bacterial load and composition eclipse seasonal and lifetime dynamics in wild meerkats, Suricata suricatta,"** contains data and code associated with a study published in *Nature Communications* (2021). The research investigates the dynamics of gut microbiomes in wild meerkats, focusing on the significant diurnal variations that surpass seasonal and lifetime changes.

## Data Origin

The core of this dataset is the `processed_data_phyloseq.RDS` object, which is a phyloseq object containing processed data for 1109 samples. This object includes:

- **ASV table**: Amplicon Sequence Variants quantified across samples.
- **taxonomic classifications**: Taxonomy assigned to each ASV.
- **phylogenetic tree**: A tree structure representing the evolutionary relationships among ASVs.
- **sample metadata**: Relevant metadata for each sample used in the analysis.

To facilitate further analysis and integration, the `processed_data_phyloseq.RDS` object has been converted into a `TreeSummarizedExperiment` object; `Meerkat2021.rds` using [mia](https://microbiome.github.io/mia/) tools. 

## Accessing the Data

For additional details, raw data, and code, please refer to the original zenodo repository:

[Original data and code for Diurnal oscillations in gut bacterial load and composition](https://zenodo.org/records/5337076)

### References 

- Alice Risely, Kerstin Wilhelm, Marta B. Manser, Tim Clutton-Brock, & Simone Sommer. (2021). Data and code for: Diurnal oscillations in gut bacterial load and composition eclipse seasonal and lifetime dynamics in wild meerkats, Suricata suricatta (Version 1.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5337076
