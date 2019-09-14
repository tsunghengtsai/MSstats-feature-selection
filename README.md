# MSstats-feature-selection

This repository hosts a set of R scripts for the MSstats analyses and for evaluating and summarizing the analysis results presented in the manuscript: 

	    T-H Tsai, M Choi, B Banfai, Y Liu, T Dunkley, and O Vitek (2019),
	    "Selection of features with consistent profiles improves relative protein quantification in mass spectrometry experiments."

## Datasets

### Benchmarks

- DIA Bruderer (MassIVE ID: MSV000081828)
- DIA Navarro (MassIVE ID: MSV000081024)
- DDA Choi (MassIVE ID: MSV000084181)
- DDA iPRG (MassIVE ID: MSV000079843)

### Biological investigations

- DIA Selevsek (MassIVE ID: MSV000081677)

## Step 1: Performing MSstats analyses

MSstats analyses are performed with the options of using 

- all features (default)
- informative features (`featureSubset = "highQuality"`, `remove_uninformative_feature_outlier = TRUE`)
- top-n features (`featureSubset = "topN"`)

### DIA and DDA Benchmarks

Run `run-msstats-benchmark.R` to perform MSstats analyses of the datasets of the benchmark controlled mixtures.

### Biological investigations

Run `run-msstats-selevsek.R` to perform MSstats analyses of the datasets from the DIA Selevsek investigation. 

## Step 2: Evaluation

### DIA benchmarks

Run `eval-DIAbenchmark.R` to evaluate the analysis results of the DIA benchmarks. The evaluation is based on estimation accuracy, detection of differential abundancem and reproducibility across data processing tools. The ground-truth defined in 

- `gold_dia_bruderer.R`

### DDA benchmarks

Run `eval-DDAbenchmark.R` to evaluate the analysis results of the DDA benchmarks. The evaluation is based on estimation accuracy, detection of differential abundancem and reproducibility across data processing tools. The ground-truth defined in 

- `gold_dda_choi.R`
- `gold_dda_iprg.R`

### Biological investigations

Run `eval-selevsek.R` to evaluate the analysis results of the Selevsek's datasets. The evaluation focuses on impact of data processing options.
