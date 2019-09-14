######################################################################
## Evaluation of the MSstats analyses for the Selevsek experiment
##   - DIA Selevsek (MassIVE ID: MSV000081677)
## 
## Please run first the MSstats analyses with run-msstats-selevsek.R
## 
## The results are summarized and presented in the manuscript:
##   T-H Tsai, M Choi, B Banfai, Y Liu, T Dunkley, and O Vitek (2019),
##   "Selection of features with consistent profiles improves relative
##   protein quantification in mass spectrometry experiments."
######################################################################

library(tidyverse)

meth <- c("All", "Proposed")
meth_code <- c("all", "inf")

dset_names <- c("dia.selevsek.full_sparse", "dia.selevsek.lowcv_sparse", 
                "dia.selevsek.full_50pct", "dia.selevsek.lowcv_50pct")
test_all <- vector("list", length(meth_code))
for (m in seq_along(meth_code)) {
    test_onemeth <- vector("list", length(dset_names))
    for (ii in seq_along(dset_names)) {
        dset_name <- dset_names[ii]
        test_name <- paste0("test/test.", dset_name, ".", meth_code[m])
        load(paste0(test_name, ".rda"))
        
        test_onemeth[[ii]] <- as_tibble(test) %>% 
            select(protein = Protein, label = Label, log2FC, pvalue, adj.pvalue) %>% 
            mutate_if(is.factor, as.character) %>% 
            mutate(dset = dset_name)
    }
    test_all[[m]] <- bind_rows(test_onemeth) %>% 
        mutate(method = meth[m])
}
test_all <- bind_rows(test_all)

# Number of statistically significant changes (Fig. 6a)
test_all %>% 
    filter(method == "All") %>% 
    filter(!is.na(pvalue), adj.pvalue < 0.05) %>% 
    count(dset, label)

# Agreement across data processing options --------------------------------

library(eulerr)
library(grid)
library(gridExtra)

# Euler diagram for all-features and proposed methods (Fig. 6e, Fig. 6f)
uniq_ctrx <- unique(test_all$label)
for (i in seq_along(uniq_ctrx)) {
    test_tmp <- test_all %>% 
        filter(label == uniq_ctrx[i]) %>% 
        filter(!is.na(pvalue)) %>% 
        filter(adj.pvalue < 0.05) %>% 
        mutate(input = ifelse(str_detect(dset, "full"), "full", "lowcv")) %>% 
        mutate(Qfilter = ifelse(str_detect(dset, "sparse"), "sparse", "half")) %>% 
        mutate(de = TRUE)
    
    # Testing results across processing options
    test_xcv <- bind_rows(
        test_tmp %>% 
            filter(method == "All", Qfilter == "half") %>% 
            select(protein, input, de) %>% 
            spread(input, de, fill = FALSE) %>% 
            mutate(method = "All", Qfilter = "half"), 
        test_tmp %>% 
            filter(method == "Proposed", Qfilter == "half") %>% 
            select(protein, input, de) %>% 
            spread(input, de, fill = FALSE) %>% 
            mutate(method = "Proposed", Qfilter = "half"), 
        test_tmp %>% 
            filter(method == "All", Qfilter == "sparse") %>% 
            select(protein, input, de) %>% 
            spread(input, de, fill = FALSE) %>% 
            mutate(method = "All", Qfilter = "sparse"), 
        test_tmp %>% 
            filter(method == "Proposed", Qfilter == "sparse") %>% 
            select(protein, input, de) %>% 
            spread(input, de, fill = FALSE) %>% 
            mutate(method = "Proposed", Qfilter = "sparse")
    )
    
    # Euler diagram 
    for (j in seq_along(meth)) {
        test_tmp <- test_xcv %>% 
            filter(method == meth[j])
        test_tmp_half <- test_tmp %>% 
            filter(Qfilter == "half") %>% 
            rename(`Full_50%` = full, `LowCV_50%` = lowcv) %>% 
            select(-Qfilter)
        test_tmp_sparse <- test_tmp %>% 
            filter(Qfilter == "sparse") %>% 
            rename(Full_Sparse = full, LowCV_Sparse = lowcv) %>% 
            select(-Qfilter)
        test_tmp2 <- full_join(test_tmp_half, test_tmp_sparse) %>% 
            select(-protein, -method)
        test_tmp2[is.na(test_tmp2)] <- FALSE
        
        fit2 <- euler(as.matrix(test_tmp2))
        obj_ep <- plot(fit2, quantities = list(fontsize = 28, col = "#2F4F4F"), 
                       lty = 1:4, labels = list(fontsize = 28))
        dset_meth <- paste0(uniq_ctrx[i], " (", meth[j], ")")
        pngname <- paste0("repsigcv_", str_replace(uniq_ctrx[i], "-", ""), "_", 
                          str_to_lower(meth[j]), ".png")
        png(file = pngname, width=800, height=600)
        print(grid.arrange(
            grobs = list(obj_ep),
            top = textGrob(dset_meth, x = 0.02, y = 0.2, just = "left", gp=gpar(fontsize=36,font=2))
        ))
        dev.off()
    }
}

# Example profile plot with testing result --------------------------------

plot_runsum <- function(one_processed, one_runsum, text_title) {
    one_processed %>% 
        ggplot(aes(run, log2inty, group = feature)) + 
        geom_point(aes(shape = is_outlier), color = "gray") + 
        geom_line(aes(linetype = feature_quality, alpha = feature_quality), color = "gray", size = 0.5) + 
        geom_point(data = one_runsum, aes(run, log2inty, color = method), size = 2.5) + 
        geom_line(data = one_runsum, aes(run, log2inty, group = method, color = method), size = 1.5) + 
        geom_vline(xintercept = c(3.5, 6.5, 9.5, 12.5, 15.5), linetype = "dashed") + 
        annotate("text", x = c(2, 5, 8, 11, 14, 17), y = 21.5, size = 5, 
                 label = paste0("T", 0:5, "\n(", c(0, 15, 30, 60, 90, 120)," min)")) + 
        coord_cartesian(ylim = c(0, 22.5)) +
        scale_x_discrete(breaks = c(3, 6, 9, 12, 15, 18)) +
        scale_color_manual(values = c("#666666", "#E69F00"), breaks=c("Proposed", "All")) + 
        scale_shape_manual(values = c(16, 1), guide = FALSE) + 
        scale_linetype_discrete(guide = FALSE) + 
        scale_alpha_manual(values = c(1, 0.5), guide = FALSE) + 
        labs(x = "Run", y = "Log2-intensity", color = "") + 
        ggtitle(text_title) + 
        theme_bw() + 
        theme(text = element_text(size = 18)) + 
        theme(legend.position = c(0.8, -0.117), legend.direction = "horizontal") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

# Profile plot for the example protein (Fig. 6b, Fig. 6c, Fig. 6d)
dset_names <- c("dia.selevsek.full_sparse", "dia.selevsek.lowcv_sparse", 
                "dia.selevsek.full_50pct", "dia.selevsek.lowcv_50pct")
dset_text <- c("Full_Sparse", "LowCV_Sparse", "Full_50%", "LowCV_50%")

ii <- 1  # (ii in 1:3)

# Processed data with all the features
load(paste0("processed/processed.", dset_names[ii], ".all.rda"))
processed.all <- processed

# Processed data with the selected informative features
load(paste0("processed/processed.", dset_names[ii], ".inf.rda"))

# Peak intensity
pdata <- as_tibble(processed$ProcessedData) %>% 
    select(
        protein = PROTEIN, peptide = PEPTIDE, feature = FEATURE, group = GROUP, 
        run = RUN, log2inty = ABUNDANCE, censored, feature_quality, is_outlier
    ) %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate(run = factor(run, levels = as.character(1:n_distinct(processed$ProcessedData$RUN))))

# Protein-level summarization
runsum <- bind_rows(
    as_tibble(processed.all$RunlevelData) %>% 
        mutate_if(is.factor, as.character) %>% 
        mutate(method = "All"), 
    as_tibble(processed$RunlevelData) %>% 
        mutate_if(is.factor, as.character) %>% 
        mutate(method = "Proposed")
) %>% 
    select(protein = Protein, run = RUN, log2inty = LogIntensities, method) %>% 
    mutate(run = factor(run, levels = as.character(1:n_distinct(processed$ProcessedData$RUN))))

# Example protein
prot_name <- "YPL117C"

# Processed and protein-level data for the example protein
one_processed <- pdata %>% 
    filter(protein == prot_name)
one_runsum <- runsum %>% 
    filter(protein == prot_name) %>% 
    mutate(feature = method)

# Profile plot
plot_runsum(one_processed, one_runsum, paste0(dset_text[ii], ", ", prot_name))

# Summary of relative protein quantification
test_all %>%
    filter(dset == dset_names[ii]) %>%
    filter(protein == prot_name, label == "T2-T0") %>%
    mutate(method = factor(method, levels = c("Proposed", "All"))) %>%
    arrange(method)
