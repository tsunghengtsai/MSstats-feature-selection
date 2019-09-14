######################################################################
## MSstats analyses with different options on the following benchmarks
##   - DDA Choi (MassIVE ID: MSV000084181)
##   - DDA iPRG (MassIVE ID: MSV000079843)
##   - DIA Bruderer (MassIVE ID: MSV000081828)
##   - DIA Navarro (MassIVE ID: MSV000081024)
######################################################################

library(MSstats)

if (!dir.exists("processed")) dir.create("processed")
if (!dir.exists("summarized")) dir.create("summarized")
if (!dir.exists("test")) dir.create("test")

# Running MSstats for DDA Choi with different options ---------------------

dset <- c("dda.choi.mq", "dda.choi.pg", "dda.choi.pd", "dda.choi.sl")
cval <- c("NA", "0", "NA", "0")

for (i in seq_along(dset)) {
    # Loading data 
    dset_name <- paste0("input.", dset[i])
    load(paste0("data/", dset_name, ".rda"))
    eval(parse(text = paste0("quant <- ", dset_name)))
    
    # Analysis with all features
    processed <- dataProcess(
        quant,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = cval[i],
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )
    
    # Testing
    comparison1 <- matrix(c(1, -1, 0, 0, 0), nrow = 1)
    comparison2 <- matrix(c(1, 0, -1, 0, 0), nrow = 1)
    comparison3 <- matrix(c(1, 0, 0, -1, 0), nrow = 1)
    comparison4 <- matrix(c(1, 0, 0, 0, -1), nrow = 1)
    comparison5 <- matrix(c(0, 1, -1, 0, 0), nrow = 1)
    comparison6 <- matrix(c(0, 1, 0, -1, 0), nrow = 1)
    comparison7 <- matrix(c(0, 1, 0, 0, -1), nrow = 1)
    comparison8 <- matrix(c(0, 0, 1, -1, 0), nrow = 1)
    comparison9 <- matrix(c(0, 0, 1, 0, -1), nrow = 1)
    comparison10 <- matrix(c(0, 0, 0, 1, -1), nrow = 1)
    comparison <- rbind(comparison1, comparison2, comparison3, comparison4, 
                        comparison5, comparison6, comparison7, comparison8, 
                        comparison9, comparison10)
    row.names(comparison) <- c("M1-M2", "M1-M3", "M1-M4", "M1-M5", "M2-M3", 
                               "M2-M4", "M2-M5", "M3-M4", "M3-M5", "M4-M5")
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult

    # Saving intermediate/testing results
    summarized_name <- paste0("summarized/summarized.", dset[i], ".all.rda")
    eval(parse(text = paste0("summarized.", dset[i], ".all <- processed$RunlevelData")))
    eval(parse(text = paste0("save(summarized.", dset[i], ".all, file=summarized_name)")))
    test_name <- paste0("test/test.", dset[i], ".all.rda")
    eval(parse(text = paste0("test.", dset[i], ".all <- test")))
    eval(parse(text = paste0("save(test.", dset[i], ".all, file=test_name)")))
    
    # Analysis with the proposed approach
    processed <- dataProcess(
        quant,
        featureSubset = "highQuality",
        remove_uninformative_feature_outlier = TRUE,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = cval[i], 
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult

    # Saving intermediate/testing results
    processed_name <- paste0("processed/processed.", dset[i], ".inf.rda")
    eval(parse(text = paste0("processed.", dset[i], ".inf <- processed")))
    eval(parse(text = paste0("save(processed.", dset[i], ".inf, file=processed_name)")))
    test_name <- paste0("test/test.", dset[i], ".inf.rda")
    eval(parse(text = paste0("test.", dset[i], ".inf <- test")))
    eval(parse(text = paste0("save(test.", dset[i], ".inf, file=test_name)")))
    
    # Analysis with top-n features
    for (n in c(3, 5)) {
        processed <- dataProcess(
            quant,
            featureSubset = "topN",
            n_top_feature = n,
            summaryMethod = "TMP",
            cutoffCensored = "minFeature",
            censoredInt = cval[i],
            MBimpute = TRUE,
            maxQuantileforCensored = 0.999
        )
        test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult
        summarized_name <- paste0("summarized/summarized.", dset[i], ".top", n, ".rda")
        eval(parse(text = paste0("summarized.", dset[i], ".top", n, " <- processed$RunlevelData")))
        eval(parse(text = paste0("save(summarized.", dset[i], ".top", n, ", file=summarized_name)")))
        test_name <- paste0("test/test.", dset[i], ".top", n, ".rda")
        eval(parse(text = paste0("test.", dset[i], ".top", n, " <- test")))
        eval(parse(text = paste0("save(test.", dset[i], ".top", n, ", file=test_name)")))
    }
}

# Running MSstats for DDA iPRG with different options ---------------------

dset <- c("dda.iprg.mq", "dda.iprg.pg", "dda.iprg.sl")
cval <- c("NA", "0", "0")

processed <- vector("list", length(dset))
test <- vector("list", length(dset))
for (i in seq_along(dset)) {
    # Loading data 
    dset_name <- paste0("input.", dset[i])
    load(paste0("data/", dset_name, ".rda"))
    eval(parse(text = paste0("quant <- ", dset_name)))
    
    # Analysis with all features
    processed <- dataProcess(
        quant,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = cval[i],
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )

    # Testing
    comparison1 <- matrix(c(-1,1,0,0), nrow=1)
    comparison2 <- matrix(c(-1,0,1,0), nrow=1)
    comparison3 <- matrix(c(-1,0,0,1), nrow=1)
    comparison4 <- matrix(c(0,-1,1,0), nrow=1)
    comparison5 <- matrix(c(0,-1,0,1), nrow=1)
    comparison6 <- matrix(c(0,0,-1,1), nrow=1)
    comparison <- rbind(comparison1, comparison2, comparison3, comparison4, 
                        comparison5, comparison6)
    row.names(comparison) <- c("C2-C1", "C3-C1", "C4-C1", "C3-C2", "C4-C2", "C4-C3")
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult

    # Saving intermediate/testing results
    summarized_name <- paste0("summarized/summarized.", dset[i], ".all.rda")
    eval(parse(text = paste0("summarized.", dset[i], ".all <- processed$RunlevelData")))
    eval(parse(text = paste0("save(summarized.", dset[i], ".all, file=summarized_name)")))
    test_name <- paste0("test/test.", dset[i], ".all.rda")
    eval(parse(text = paste0("test.", dset[i], ".all <- test")))
    eval(parse(text = paste0("save(test.", dset[i], ".all, file=test_name)")))
    
    # Analysis with the proposed approach
    processed <- dataProcess(
        quant,
        featureSubset = "highQuality",
        remove_uninformative_feature_outlier = TRUE,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = cval[i], 
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult

    # Saving intermediate/testing results
    processed_name <- paste0("processed/processed.", dset[i], ".inf.rda")
    eval(parse(text = paste0("processed.", dset[i], ".inf <- processed")))
    eval(parse(text = paste0("save(processed.", dset[i], ".inf, file=processed_name)")))
    test_name <- paste0("test/test.", dset[i], ".inf.rda")
    eval(parse(text = paste0("test.", dset[i], ".inf <- test")))
    eval(parse(text = paste0("save(test.", dset[i], ".inf, file=test_name)")))
    
    # Analysis with top-n features
    for (n in c(3, 5)) {
        processed <- dataProcess(
            quant,
            featureSubset = "topN",
            n_top_feature = n,
            summaryMethod = "TMP",
            cutoffCensored = "minFeature",
            censoredInt = cval[i],
            MBimpute = TRUE,
            maxQuantileforCensored = 0.999
        )
        test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult
        summarized_name <- paste0("summarized/summarized.", dset[i], ".top", n, ".rda")
        eval(parse(text = paste0("summarized.", dset[i], ".top", n, " <- processed$RunlevelData")))
        eval(parse(text = paste0("save(summarized.", dset[i], ".top", n, ", file=summarized_name)")))
        test_name <- paste0("test/test.", dset[i], ".top", n, ".rda")
        eval(parse(text = paste0("test.", dset[i], ".top", n, " <- test")))
        eval(parse(text = paste0("save(test.", dset[i], ".top", n, ", file=test_name)")))
    }
}

# Running MSstats for DIA Bruderer with different options -----------------

dset <- c("dia.bruderer.du", "dia.bruderer.sl", "dia.bruderer.sn")
cval <- c("NA", "0", "0")

processed <- vector("list", length(dset))
test <- vector("list", length(dset))
for (i in seq_along(dset)) {
    # Loading data 
    dset_name <- paste0("input.", dset[i])
    load(paste0("data/", dset_name, ".rda"))
    eval(parse(text = paste0("quant <- ", dset_name)))
    
    # Analysis with all features
    processed <- dataProcess(
        quant,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = cval[i],
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )

    # Testing
    comparison1 <- matrix(c(-1,1,0,0,0,0,0), nrow=1)
    comparison2 <- matrix(c(-1,0,1,0,0,0,0), nrow=1)
    comparison3 <- matrix(c(-1,0,0,1,0,0,0), nrow=1)
    comparison4 <- matrix(c(-1,0,0,0,1,0,0), nrow=1)
    comparison5 <- matrix(c(-1,0,0,0,0,1,0), nrow=1)
    comparison6 <- matrix(c(-1,0,0,0,0,0,1), nrow=1)
    comparison8 <- matrix(c(0,-1,1,0,0,0,0), nrow=1)
    comparison9 <- matrix(c(0,-1,0,1,0,0,0), nrow=1)
    comparison10 <- matrix(c(0,-1,0,0,1,0,0), nrow=1)
    comparison11 <- matrix(c(0,-1,0,0,0,1,0), nrow=1)
    comparison12 <- matrix(c(0,-1,0,0,0,0,1), nrow=1)
    comparison14 <- matrix(c(0,0,-1,1,0,0,0), nrow=1)
    comparison15 <- matrix(c(0,0,-1,0,1,0,0), nrow=1)
    comparison16 <- matrix(c(0,0,-1,0,0,1,0), nrow=1)
    comparison17 <- matrix(c(0,0,-1,0,0,0,1), nrow=1)
    comparison19 <- matrix(c(0,0,0,-1,1,0,0), nrow=1)
    comparison20 <- matrix(c(0,0,0,-1,0,1,0), nrow=1)
    comparison21 <- matrix(c(0,0,0,-1,0,0,1), nrow=1)
    comparison23 <- matrix(c(0,0,0,0,-1,1,0), nrow=1)
    comparison24 <- matrix(c(0,0,0,0,-1,0,1), nrow=1)
    comparison26 <- matrix(c(0,0,0,0,0,-1,1), nrow=1)
    comparison <- rbind(comparison1, comparison2, comparison3, comparison4,
                        comparison5, comparison6, comparison8, comparison9, 
                        comparison10, comparison11, comparison12, comparison14,
                        comparison15, comparison16, comparison17, comparison19, 
                        comparison20, comparison21, comparison23, comparison24,
                        comparison26)
    row.names(comparison) <- c("S2-S1", "S3-S1", "S4-S1", "S5-S1", "S6-S1",
                               "S7-S1", "S3-S2", "S4-S2", "S5-S2", "S6-S2",
                               "S7-S2", "S4-S3", "S5-S3", "S6-S3", "S7-S3",
                               "S5-S4", "S6-S4", "S7-S4", "S6-S5", "S7-S5",
                               "S7-S6")
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult

    # Saving intermediate/testing results
    summarized_name <- paste0("summarized/summarized.", dset[i], ".all.rda")
    eval(parse(text = paste0("summarized.", dset[i], ".all <- processed$RunlevelData")))
    eval(parse(text = paste0("save(summarized.", dset[i], ".all, file=summarized_name)")))
    test_name <- paste0("test/test.", dset[i], ".all.rda")
    eval(parse(text = paste0("test.", dset[i], ".all <- test")))
    eval(parse(text = paste0("save(test.", dset[i], ".all, file=test_name)")))
    
    # Analysis with the proposed approach
    processed <- dataProcess(
        quant,
        featureSubset = "highQuality",
        remove_uninformative_feature_outlier = TRUE,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = cval[i], 
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult

    # Saving intermediate/testing results
    processed_name <- paste0("processed/processed.", dset[i], ".inf.rda")
    eval(parse(text = paste0("processed.", dset[i], ".inf <- processed")))
    eval(parse(text = paste0("save(processed.", dset[i], ".inf, file=processed_name)")))
    test_name <- paste0("test/test.", dset[i], ".inf.rda")
    eval(parse(text = paste0("test.", dset[i], ".inf <- test")))
    eval(parse(text = paste0("save(test.", dset[i], ".inf, file=test_name)")))
    
    # Analysis with top-n features
    for (n in c(3, 5, 10)) {
        processed <- dataProcess(
            quant,
            featureSubset = "topN",
            n_top_feature = n,
            summaryMethod = "TMP",
            cutoffCensored = "minFeature",
            censoredInt = cval[i],
            MBimpute = TRUE,
            maxQuantileforCensored = 0.999
        )
        test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult
        summarized_name <- paste0("summarized/summarized.", dset[i], ".top", n, ".rda")
        eval(parse(text = paste0("summarized.", dset[i], ".top", n, " <- processed$RunlevelData")))
        eval(parse(text = paste0("save(summarized.", dset[i], ".top", n, ", file=summarized_name)")))
        test_name <- paste0("test/test.", dset[i], ".top", n, ".rda")
        eval(parse(text = paste0("test.", dset[i], ".top", n, " <- test")))
        eval(parse(text = paste0("save(test.", dset[i], ".top", n, ", file=test_name)")))
    }
}

# Running MSstats for DIA Navarro with different options ------------------

dset <- c("dia.navarro.du", "dia.navarro.os", "dia.navarro.sl", "dia.navarro.sn")
cval <- c("NA", "0", "0", "0")

processed <- vector("list", length(dset))
test <- vector("list", length(dset))
for (i in seq_along(dset)) {
    # Loading data 
    dset_name <- paste0("input.", dset[i])
    load(paste0("data/", dset_name, ".rda"))
    eval(parse(text = paste0("quant <- ", dset_name)))
    
    # Analysis with all features
    processed <- dataProcess(
        quant,
        normalization = FALSE,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = cval[i],
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )

    # Testing
    comparison <- matrix(c(1,-1), nrow=1)
    row.names(comparison) <- c("A-B")
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult

    # Saving intermediate/testing results
    summarized_name <- paste0("summarized/summarized.", dset[i], ".all.rda")
    eval(parse(text = paste0("summarized.", dset[i], ".all <- processed$RunlevelData")))
    eval(parse(text = paste0("save(summarized.", dset[i], ".all, file=summarized_name)")))
    test_name <- paste0("test/test.", dset[i], ".all.rda")
    eval(parse(text = paste0("test.", dset[i], ".all <- test")))
    eval(parse(text = paste0("save(test.", dset[i], ".all, file=test_name)")))
    
    # Analysis with the proposed approach
    processed <- dataProcess(
        quant,
        normalization = FALSE,
        featureSubset = "highQuality",
        remove_uninformative_feature_outlier = TRUE,
        summaryMethod = "TMP",
        cutoffCensored = "minFeature",
        censoredInt = cval[i], 
        MBimpute = TRUE,
        maxQuantileforCensored = 0.999
    )
    test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult

    # Saving intermediate/testing results
    processed_name <- paste0("processed/processed.", dset[i], ".inf.rda")
    eval(parse(text = paste0("processed.", dset[i], ".inf <- processed")))
    eval(parse(text = paste0("save(processed.", dset[i], ".inf, file=processed_name)")))
    test_name <- paste0("test/test.", dset[i], ".inf.rda")
    eval(parse(text = paste0("test.", dset[i], ".inf <- test")))
    eval(parse(text = paste0("save(test.", dset[i], ".inf, file=test_name)")))
    
    # Analysis with top-n features
    for (n in c(3, 5, 10)) {
        processed <- dataProcess(
            quant,
            normalization = FALSE,
            featureSubset = "topN",
            n_top_feature = n,
            summaryMethod = "TMP",
            cutoffCensored = "minFeature",
            censoredInt = cval[i],
            MBimpute = TRUE,
            maxQuantileforCensored = 0.999
        )
        test <- groupComparison(contrast.matrix = comparison, data = processed)$ComparisonResult
        summarized_name <- paste0("summarized/summarized.", dset[i], ".top", n, ".rda")
        eval(parse(text = paste0("summarized.", dset[i], ".top", n, " <- processed$RunlevelData")))
        eval(parse(text = paste0("save(summarized.", dset[i], ".top", n, ", file=summarized_name)")))
        test_name <- paste0("test/test.", dset[i], ".top", n, ".rda")
        eval(parse(text = paste0("test.", dset[i], ".top", n, " <- test")))
        eval(parse(text = paste0("save(test.", dset[i], ".top", n, ", file=test_name)")))
    }
}
