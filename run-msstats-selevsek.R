#######################################################################
## MSstats analyses with different options on the biological experiment
##   - DIA Selevsek (MassIVE ID: MSV000081677)
#######################################################################

# Datasets from two Skyline analyses and two q-filtering strategies -------

library(MSstats)
library(dplyr)

if (!dir.exists("data")) dir.create("data")

# Full analysis (Reanalysis ID on MassIVE: RMSV000000251.2)
raw <- read.csv("data/Selevsek2015-MSstats-input-90nodup-i.csv")
annotation <- read.csv("data/Selevsek2015_DIA_Skyline_all_annotation.csv")

# Full analysis with sparse q-filter
quant <- SkylinetoMSstatsFormat(
    raw,
    annotation = annotation,
    filter_with_Qvalue = TRUE,   # same as default
    qvalue_cutoff = 0.01,        # same as default
    fewMeasurements = "remove",  # same as default
    removeProtein_with1Feature = TRUE
)
save(quant, file = "data/input.dia.selevsek.full_sparse.rda")

# Full analysis with 50% q-filter
raw <- raw[, -2]
raw <- raw %>%
    filter(StandardType != "iRT") %>%
    mutate(DetectionQValue = as.numeric(as.character(DetectionQValue))) %>%
    group_by(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge) %>%
    filter(sum(DetectionQValue < 0.01) > 9) %>%
    ungroup()
raw <- as.data.frame(raw)
quant <- SkylinetoMSstatsFormat(
    raw,
    annotation = annotation,
    filter_with_Qvalue = FALSE,  # disabled, as 50% q-filter has been applied
    fewMeasurements = "remove",  # same as default
    removeProtein_with1Feature = TRUE
)
save(quant, file = "data/input.dia.selevsek.full_50pct.rda")

# LowCV analysis (Reanalysis ID on MassIVE: RMSV000000251.1)
raw <- read.csv("data/Selevsek2015-MSstats-input-90lowcv-i.csv")
annotation <- read.csv("data/Selevsek2015_DIA_Skyline_lowcv_annotation.csv")

# LowCV analysis with sparse q-filter
quant <- SkylinetoMSstatsFormat(
    raw,
    annotation = annotation,
    filter_with_Qvalue = TRUE,   # same as default
    qvalue_cutoff = 0.01,        # same as default
    fewMeasurements = "remove",  # same as default
    removeProtein_with1Feature = TRUE
)
save(quant, file = "data/input.dia.selevsek.lowcv_sparse.rda")

# LowCV analysis with 50% q-filter
raw <- raw[,-2]
raw <- raw %>%
    filter(StandardType != "iRT") %>%
    mutate(DetectionQValue = as.numeric(as.character(DetectionQValue))) %>%
    group_by(ProteinName, PeptideModifiedSequence, PrecursorCharge, FragmentIon, ProductCharge) %>%
    filter(sum(DetectionQValue < 0.01) > 9) %>%
    ungroup()
raw <- as.data.frame(raw)
quant <- SkylinetoMSstatsFormat(
    raw,
    annotation = annotation,
    filter_with_Qvalue = FALSE,  # disabled, as 50% q-filter has been applied
    fewMeasurements = "remove",  # same as default
    removeProtein_with1Feature = TRUE
)
save(quant, file = "data/input.dia.selevsek.lowcv_50pct.rda")

# Running MSstats with different options ----------------------------------

library(MSstats)

if (!dir.exists("processed")) dir.create("processed")
if (!dir.exists("test")) dir.create("test")

dset <- c("dia.selevsek.full_sparse", "dia.selevsek.lowcv_sparse",
          "dia.selevsek.full_50pct", "dia.selevsek.lowcv_50pct")

for (i in seq_along(dset)) {
    # Loading data 
    load(paste0("data/input.", dset[i], ".rda"))
    
    # Analysis with all features
    processed <- dataProcess(
        quant,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = "0",
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )
    # Testing
    comparison1 <- matrix(c(-1, 1, 0, 0, 0, 0), nrow = 1)
    comparison2 <- matrix(c(-1, 0, 1, 0, 0, 0), nrow = 1)
    comparison3 <- matrix(c(-1, 0, 0, 1, 0, 0), nrow = 1)
    comparison4 <- matrix(c(-1, 0, 0, 0, 1, 0), nrow = 1)
    comparison5 <- matrix(c(-1, 0, 0, 0, 0, 1), nrow = 1)
    comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5)
    row.names(comparison) <- c("T1-T0", "T2-T0", "T3-T0", "T4-T0", "T5-T0")
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult
    
    # Saving intermediate/testing results
    processed_name <- paste0("processed/processed.", dset[i], ".all.rda")
    save(processed, file = processed_name)
    test_name <- paste0("test/test.", dset[i], ".all.rda")
    save(test, file = test_name)
    
    # Analysis with selected informative features
    processed <- dataProcess(
        quant,
        featureSubset = "highQuality",
        remove_uninformative_feature_outlier = TRUE,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = "0", 
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult
    processed_name <- paste0("processed/processed.", dset[i], ".inf.rda")
    save(processed, file = processed_name)
    test_name <- paste0("test/test.", dset[i], ".inf.rda")
    save(test, file = test_name)
}
