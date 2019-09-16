#######################################################################
## MSstats analysis with feature selection on the biological experiment
##   - SRM Dunkley (Panorama link: https://panoramaweb.org/neuroSRM.url)
#######################################################################


# Running MSstats with feature selection option ---------------------------

library(MSstats)

if (!dir.exists("processed")) dir.create("processed")

raw <- read.csv("data/SRM_Dunkley/MSstats Input.csv")

# Convert Skyline output to MSstats - iRT standard and proteins with only one 
# feature are removed
quant <- SkylinetoMSstatsFormat(
    raw,
    removeiRT = TRUE,
    filter_with_Qvalue = FALSE,
    removeProtein_with1Feature = TRUE,
    fewMeasurements = "remove"
)

# Data processing with feature selection to flag uninformative features/outliers
processed <- dataProcess(
    quant,
    featureSubset = "highQuality",
    summaryMethod = "TMP",
    cutoffCensored = "minFeature",
    censoredInt = "0", 
    MBimpute = TRUE,
    maxQuantileforCensored = 0.999
)
save(processed, file = "processed/processed.srm.dunkley.sl.inf.rda")
