######################################################################
## Evaluation of the MSstats analyses for the DIA benchmarks 
##   - DIA Bruderer (MassIVE ID: MSV000081828)
##   - DIA Navarro (MassIVE ID: MSV000081024)
## 
## The evaluation is based on the ground-truth defined in 
##   - gold_dia_bruderer.R (for DIA Bruderer)
## 
## Please run first the MSstats analyses with run-msstats-benchmark.R
## 
## The results are summarized and presented in the manuscript:
##   T-H Tsai, M Choi, B Banfai, Y Liu, B MacLean, T Dunkley, and O Vitek (2019),
##   "Selection of features with consistent profiles improves relative
##   protein quantification in mass spectrometry experiments."
######################################################################

library(tidyverse)

list_dset <- list(
    c("dia.bruderer.du", "dia.bruderer.sl", "dia.bruderer.sn"),
    c("dia.navarro.du", "dia.navarro.os", "dia.navarro.sl", "dia.navarro.sn")
)

# Number of informative features ------------------------------------------

nb_ftr_all <- vector("list", length(list_dset))
for (ds in seq_along(list_dset)) {
    dset_names <- list_dset[[ds]]
    nb_ftr <- vector("list", length(dset_names))
    for (ii in seq_along(dset_names)) {
        dset_name <- dset_names[ii]
        proc_name <- paste0("processed.", dset_name, ".inf")
        load(paste0("processed/", proc_name, ".rda"))
        eval(parse(text = paste0("processed <- ", proc_name)))
        eval(parse(text = paste0("rm(", proc_name, ")")))
        
        nb_ftr[[ii]] <- processed$ProcessedData %>% 
            as_tibble() %>% 
            select(protein = PROTEIN, feature = FEATURE, feature_quality) %>% 
            mutate_if(is.factor, as.character) %>% 
            distinct(protein, feature, feature_quality) %>% 
            group_by(protein) %>% 
            summarise(All = n(), Proposed = sum(feature_quality == "Informative")) %>% 
            mutate(tool = str_sub(dset_name, start = -2, end = -1))
    }
    nb_ftr_all[[ds]] <- bind_rows(nb_ftr) %>% 
        mutate(expt = str_sub(dset_names[1], start = 1, end = -4))
}
nb_ftr_all <- bind_rows(nb_ftr_all) %>% 
    mutate(expt = expt %>% str_replace("\\.", " ") %>% str_to_title()) %>%
    mutate(expt = expt %>% str_replace_all(c("Dia" = "DIA"))) %>%
    mutate(tool = tool %>% str_replace_all(
        c("du" = "DIA Umpire", "os" = "OpenSWATH", "sl" = "Skyline", "sn" = "Spectronaut")
    ))

# Number of features before/after feature selection (Fig. 2a)
nb_ftr_all %>% 
    gather("type", "nb_feature", All, Proposed) %>% 
    mutate(type = factor(type, levels = c("Proposed", "All"))) %>% 
    ggplot(aes(tool, nb_feature)) + 
    geom_boxplot(aes(fill = type)) + 
    geom_hline(yintercept = 0, linetype = 1) +
    scale_fill_manual(values = c("#E69F00", "#999999")) +
    labs(x = "", y = "# features", fill = "") +
    facet_grid(. ~ expt, scales = "free_x", space = "free_x") + 
    coord_cartesian(ylim = c(0, 70)) + 
    theme_bw() + 
    theme(text = element_text(size = 16)) + 
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold"))
# ggsave("boxnb-dia.png", width = 12, height = 4)

# Load all group comparison results ---------------------------------------

meth <- c("All", "Proposed", "Top10", "Top5", "Top3")
meth_code <- c("all", "inf", "top10", "top5", "top3")

sensPrecCurves <- function(test_results) {
    meths <- unique(test_results$method)
    test_meth = vector("list", length = length(meths))
    for (i in seq_along(meths)) {
        test_sub <- test_results %>% 
            filter(method == meths[i], !is.na(pvalue)) %>% 
            arrange(adj.pvalue)
        ordered_de <- test_sub$de == "With change"
        test_meth[[i]] <- tibble(
            sens = cumsum(ordered_de) / sum(ordered_de), 
            fdr = cumsum(1 - ordered_de) / seq_along(ordered_de),
            fpr = cumsum(!ordered_de) / sum(!ordered_de)
        ) %>% mutate(method = meths[i])
    }
    
    return(bind_rows(test_meth))
}

test_res_all <- vector("list", length(list_dset))
sp_curves_all <- vector("list", length(list_dset))
for (ds in seq_along(list_dset)) {
    dset_names <- list_dset[[ds]]
    test_res <- vector("list", length(dset_names))
    for (ii in seq_along(dset_names)) {
        dset_name <- dset_names[ii]
        test_sub <- vector("list", length(meth_code))
        
        input_name <- paste0("input.", dset_name)
        load(paste0("data/", input_name, ".rda"))
        eval(parse(text = paste0("tmp <- ", input_name)))
        tmp <- as_tibble(tmp) %>% 
            rename_at(vars(contains("Peptide")), funs(sub("Modified", "", .))) %>% 
            rename(Protein = ProteinName) %>% 
            mutate_if(is.factor, as.character)
        tmp2 <- tmp %>% 
            group_by(Protein) %>% 
            summarise(nb_ftr = n_distinct(PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge))
        eval(parse(text = paste0("rm(", input_name, ")")))

        for (m in seq_along(meth_code)) {
            test_name <- paste0("test.", dset_name, ".", meth_code[m])
            load(paste0("test/", test_name, ".rda"))
            eval(parse(text = paste0("test_sub[[m]] <- ", test_name)))
            test_sub[[m]] <- test_sub[[m]] %>% 
                as_tibble() %>% 
                mutate_if(is.factor, as.character) %>% 
                mutate(method = meth[m], tool = str_sub(dset_name, start = -2, end = -1))
            eval(parse(text = paste0("rm(", test_name, ")")))
        }
        test_res[[ii]] <- bind_rows(test_sub) %>% 
            left_join(tmp2)
    }
    test_res <- bind_rows(test_res) %>% 
        mutate(expt = str_sub(dset_names[1], start = 1, end = -4))
    
    # Gold standard
    if (str_detect(dset_name, "navarro")) {
        gold <- test_res %>% 
            distinct(Protein) %>% 
            mutate(
                trueFC = ifelse(str_detect(Protein, "HUMAN"), 1, ifelse(str_detect(Protein, "ECOLI"), 0.25, 2)), 
                is_BG = ifelse(str_detect(Protein, "HUMAN"), TRUE, FALSE)
            )
        test_res <- test_res %>% 
            left_join(gold) %>% 
            mutate(de = ifelse(trueFC == 1, "No change", "With change")) %>% 
            mutate(trueLog2FC = log2(trueFC))
    } else if (str_detect(dset_name, "bruderer")) {
        source("gold_dia_bruderer.R")
        test_res <- bind_rows(
            test_res %>% filter(Protein %in% unlist(mix)) %>% left_join(gold) %>% mutate(is_BG = FALSE), 
            test_res %>% filter(!(Protein %in% unlist(mix))) %>% mutate(trueFC = 1, is_BG = TRUE)
        )
        test_res <- test_res %>% 
            mutate(de = ifelse(trueFC == 1, "No change", "With change")) %>% 
            mutate(trueLog2FC = log2(trueFC))
    }
    test_res_all[[ds]] <- test_res
    
    # Detection of differential abundance
    sp_curves <- vector("list", length(dset_names))
    for (ii in seq_along(dset_names)) {
        sp_curves[[ii]] <- sensPrecCurves(
            test_res %>% filter(paste0(expt, ".", tool) == dset_names[ii])
        ) %>% mutate(dset = dset_names[ii])
    }
    sp_curves_all[[ds]] <- bind_rows(sp_curves)
}
test_res_all <- bind_rows(test_res_all) %>% 
    mutate(expt = expt %>% str_replace("\\.", " ") %>% str_to_title()) %>%
    mutate(expt = expt %>% str_replace_all(c("Dia" = "DIA"))) %>%
    mutate(tool = tool %>% str_replace_all(
        c("du" = "DIA Umpire", "os" = "OpenSWATH", "sl" = "Skyline", "sn" = "Spectronaut")
    ))

sp_curves_all <- bind_rows(sp_curves_all) %>% 
    mutate(expt = str_sub(dset, start = 1, end = -4)) %>% 
    mutate(tool = str_sub(dset, start = -2, end = -1)) %>% 
    mutate(expt = expt %>% str_replace("\\.", " ") %>% str_to_title()) %>%
    mutate(expt = expt %>% str_replace_all(c("Dia" = "DIA"))) %>%
    mutate(tool = tool %>% str_replace_all(
        c("du" = "DIA Umpire", "os" = "OpenSWATH", "sl" = "Skyline", "sn" = "Spectronaut")
    ))

# Example profile plot with testing result --------------------------------

plot_ftrpair <- function(one_processed, one_runsum, text_title) {
    cbPalette <- c("#E69F00", "#666666", "#56B4E9")
    pair_processed <- bind_rows(
        one_processed %>% 
            select(protein:log2inty) %>% 
            mutate(dset = "Before"), 
        one_processed %>% 
            mutate(log2inty = ifelse(is_outlier | feature_quality == "Noninformative", NA, log2inty)) %>% 
            select(protein:log2inty) %>% 
            mutate(dset = "After")
    ) %>% 
        mutate(dset = factor(dset, levels = c("Before", "After")))
    
    one_runsum <- one_runsum %>% 
        mutate(dset = ifelse(method == "Proposed", "After", "Before")) %>% 
        mutate(dset = factor(dset, levels = c("Before", "After")), 
               method = factor(method, levels = c("Proposed", "All", "Top3")))
    
    pair_processed %>% 
        ggplot(aes(run, log2inty, group = feature)) + 
        geom_point(color = "gray", size = 1.5) + 
        geom_line(color = "gray", color = "gray", size = 0.5) + 
        geom_point(data = one_runsum, aes(run, log2inty, color = method), size = 2.5) +
        scale_color_manual(values = cbPalette) + 
        geom_line(data = one_runsum, aes(run, log2inty, group = method, color = method), size = 1.5) +
        geom_vline(xintercept = c(3.5), linetype = "dashed") + 
        annotate("text", x = 2, y = 16, label = "Condition1", size = 6) + 
        annotate("text", x = 5, y = 16, label = "Condition2", size = 6) + 
        coord_cartesian(ylim = c(0, 17)) + 
        facet_grid(dset ~ .) + 
        labs(x = "Run", y = "Log2-intensity", color = "Method") + 
        ggtitle(text_title) +
        theme_bw() + 
        theme(text = element_text(size = 16)) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold")) + 
        theme(legend.position = "bottom")
}

# Uniprot accession
regex_uniprot_iso <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}")

# Example of Skyline (Fig. 2b)
dset_name <- "dia.navarro.sl"; dset_text <- "DIA Navarro (Skyline)"
prot_name <- "1/tr|C8ZIG9|C8ZIG9_YEAS8"

# Example of Spectronaut (Fig. 2c)
# dset_name <- "dia.navarro.sn"; dset_text <- "DIA Navarro (Spectronaut)"
# prot_name <- "1/tr|C8ZIG9|C8ZIG9_YEAS8"

# Example of DIA Umpire (Fig. S4a)
# dset_name <- "dia.navarro.du"; dset_text <- "DIA Navarro (DIA Umpire)"
# prot_name <- "tr|C8ZGU7|C8ZGU7_YEAS8"

# Example of OpenSWATH (Fig. S4b)
# dset_name <- "dia.navarro.os"; dset_text <- "DIA Navarro (OpenSWATH)"
# prot_name <- "1/sp|P75913|GHRA_ECOLI"

# Loading processed data with annotation of uninformative features
proc_name <- paste0("processed.", dset_name, ".inf")
load(paste0("processed/processed.", dset_name, ".inf.rda"))
eval(parse(text = paste0("processed <- processed.", dset_name, ".inf")))

pdata <- as_tibble(processed$ProcessedData) %>% 
    select(
        protein = PROTEIN, peptide = PEPTIDE, feature = FEATURE, group = GROUP, 
        run = RUN, log2inty = ABUNDANCE, censored, feature_quality, is_outlier
    ) %>% 
    mutate_if(is.factor, as.character)
pdata <- pdata %>% 
    mutate(run = factor(run, levels = as.character(1:n_distinct(pdata$run))))

# Loading protein-level data
load(paste0("summarized/summarized.", dset_name, ".all.rda"))
eval(parse(text = paste0("run_all <- summarized.", dset_name, ".all")))
load(paste0("summarized/summarized.", dset_name, ".top3.rda"))
eval(parse(text = paste0("run_top3 <- summarized.", dset_name, ".top3")))

runsum <- bind_rows(
    as_tibble(run_all) %>% mutate_if(is.factor, as.character) %>% mutate(method = "All"), 
    as_tibble(processed$RunlevelData) %>% mutate_if(is.factor, as.character) %>% mutate(method = "Proposed"), 
    as_tibble(run_top3) %>% mutate_if(is.factor, as.character) %>% mutate(method = "Top3")
) %>% 
    select(protein = Protein, run = RUN, log2inty = LogIntensities, method) 
runsum <- runsum %>% 
    mutate(run = factor(run, levels = as.character(1:n_distinct(runsum$run))))

# Processed and protein-level data for the example protein
one_processed <- pdata %>% 
    filter(protein == prot_name)
one_runsum <- runsum %>% 
    filter(protein == prot_name) %>% 
    mutate(feature = method)

# Profile plot before and after feature selection
plot_ftrpair(one_processed, one_runsum, 
             paste0(dset_text, ", ", str_extract(prot_name, regex_uniprot_iso)))

# Summary of relative protein quantification
test_res_all %>% 
    filter(Protein == prot_name, Label == "A-B") %>% 
    filter(paste0(expt, " (", tool, ")") == dset_text) %>% 
    filter(method %in% c("Proposed", "All", "Top3")) %>% 
    select(Protein, Label, log2FC, SE, adj.pvalue, method)

# Estimation of fold change -----------------------------------------------

myPalette <- c("#E69F00", "#999999", "#7285A5", "#56B4E9", "#daeefe")
colmeth <- c("Proposed", "All", "Top10", "Top5", "Top3")

# Absolute errors for background proteins (Fig. 3a)
df_dummy <- data.frame(
    tool = "DIA Umpire", log2FC = 0.31, trueLog2FC = 0, 
    method = factor(rep(colmeth, 2), levels = colmeth), 
    expt = rep(c("DIA Bruderer", "DIA Navarro"), each = 5)
)
test_res_all %>% 
    filter(is_BG) %>% 
    filter(!is.na(pvalue)) %>% 
    mutate(method = factor(method, levels = colmeth)) %>% 
    ggplot(aes(tool, abs(log2FC - trueLog2FC))) + 
    geom_boxplot(aes(fill = method)) + 
    geom_text(aes(label = method), nudge_x = c(-0.3, -0.15, 0, 0.15, 0.3) - 0.025, 
              vjust = 0, hjust = 1, angle = 90, size = 3, data = df_dummy) +
    geom_hline(yintercept = 0, linetype = 2) + 
    scale_fill_manual(values = myPalette) + 
    labs(x = "", y = "| Deviation of estimated \nlog2-fold change from truth |", fill = "Method") +
    coord_cartesian(ylim = c(0, 0.3)) +
    facet_grid(. ~ expt, scales = "free_x", space = "free_x") + 
    theme_bw() + 
    theme(text = element_text(size = 16)) + 
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold")) + 
    theme(legend.position = "bottom")
# ggsave("boxabserr-dia_bg.png", width = 10, height = 4.5)

# Standard errors for background proteins (Fig. 3b)
df_dummy <- data.frame(
    tool = "DIA Umpire", SE = 0.31, 
    method = factor(rep(colmeth, 2), levels = colmeth), 
    expt = rep(c("DIA Bruderer", "DIA Navarro"), each = 5)
)
test_res_all %>% 
    filter(is_BG) %>% 
    filter(!is.na(pvalue)) %>% 
    mutate(method = factor(method, levels = colmeth)) %>% 
    ggplot(aes(tool, SE)) + 
    geom_boxplot(aes(fill = method)) + 
    geom_text(aes(label = method), nudge_x = c(-0.3, -0.15, 0, 0.15, 0.3) - 0.025, 
              vjust = 0, hjust = 1, angle = 90, size = 3, data = df_dummy) +
    geom_hline(yintercept = 0, linetype = 2) + 
    scale_fill_manual(values = myPalette) + 
    labs(x = "", y = " \nStandard error", fill = "Method") +
    coord_cartesian(ylim = c(0, 0.3)) +
    facet_grid(. ~ expt, scales = "free_x", space = "free_x") + 
    theme_bw() + 
    theme(text = element_text(size = 16)) + 
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold")) + 
    theme(legend.position = "bottom")
# ggsave("boxse-dia_bg.png", width = 10, height = 4.5)

# Detection of differential abundance -------------------------------------

myPalette <- c("#E69F00", "#6b6b6b", "#0E4D92", "#008ECC", "#89CFF0")
colmeth <- c("Proposed", "All", "Top10", "Top5", "Top3")

# ROC curve for the Bruderer benchmark (Fig. S5)
sp_curves_all %>% 
    filter(expt == "DIA Bruderer") %>% 
    mutate(Method = factor(method, levels = colmeth)) %>% 
    ggplot(aes(x = fpr, y = sens, group = Method, color = Method, linetype = Method, size = Method)) + 
    geom_line() + 
    scale_linetype_manual(values = c(rep("solid", 3), "dotted", "dotted")) +
    scale_size_manual(values = c(rep(0.55, 3), 0.5, 0.5)) +
    scale_color_manual(values = myPalette) + 
    coord_cartesian(xlim = c(0, 0.15), ylim = c(0.55, 0.95)) + 
    facet_wrap(~ tool) + 
    labs(x = "False positive rate", y = "Sensitivity") + 
    theme_bw() + 
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold")) + 
    theme(legend.position = "bottom")
# ggsave("rocs-dia_bruderer_doc.png", width = 6, height = 3)

# ROC curve for the Navarro benchmark (Fig. S6)
sp_curves_all %>% 
    filter(expt == "DIA Navarro") %>% 
    mutate(Method = factor(method, levels = colmeth)) %>% 
    ggplot(aes(x = fpr, y = sens, group = Method, color = Method, linetype = Method, size = Method)) + 
    geom_line() + 
    scale_linetype_manual(values = c(rep("solid", 3), "dashed", "dotted")) +
    scale_size_manual(values = c(rep(0.55, 3), 0.5, 0.5)) +
    scale_color_manual(values = myPalette) + 
    coord_cartesian(xlim = c(0, 0.15), ylim = c(0.55, 0.95)) + 
    facet_wrap(~ tool) + 
    labs(x = "False positive rate", y = "Sensitivity") + 
    theme_bw() + 
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold")) + 
    theme(legend.position = "bottom")
# ggsave("rocs-dia_navarro_doc.png", width = 6, height = 6)

# Reproducibility of detecting true changes -------------------------------

library(eulerr)
library(grid)
library(gridExtra)

regex_uniprot_iso <- regex("([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})([-]\\d{1,}){0,1}")

test_res_tp <- test_res_all %>% 
    filter(!is.na(pvalue), adj.pvalue <= 0.05) %>% 
    filter(de == "With change") %>% 
    mutate(Protein = str_extract(Protein, pattern = regex_uniprot_iso))

# Euler diagram for the Bruderer benchmark (Fig. 4a)
bm <- "DIA Bruderer"
for (j in seq_along(meth)) {
    test_tmp <- test_res_tp %>% 
        filter(expt == bm) %>% 
        filter(method == meth[j]) %>% 
        select(Protein, Label, tool) %>% 
        mutate(de = TRUE) %>% 
        spread(tool, de, fill = FALSE)
    fit2 <- euler(as.matrix(test_tmp %>% select(-Protein, -Label)))
    obj_ep <- plot(fit2, quantities = list(fontsize = 28), lty = 1:4, labels = list(fontsize = 28))
    pngname <- paste0("reptp_", str_replace(str_to_lower(bm), " ", "_"), "_", str_to_lower(meth[j]), ".png")
    dset_meth <- paste0(bm, "\n(", meth[j], ")")
    png(file = pngname, width=800, height=600)
    print(grid.arrange(grobs = list(obj_ep),
                       top = textGrob(dset_meth, x = 0.02, y = 0.2, just = "left", gp=gpar(fontsize=36,font=2))))
    dev.off()
}

# Euler diagram for the Navarro benchmark (Fig. 4b)
bm <- "DIA Navarro"
for (j in seq_along(meth)) {
    test_tmp <- test_res_tp %>% 
        filter(expt == bm, tool != "DIA Umpire") %>% 
        filter(method == meth[j]) %>% 
        select(Protein, Label, tool) %>% 
        mutate(de = TRUE) %>% 
        spread(tool, de, fill = FALSE)
    fit2 <- euler(as.matrix(test_tmp %>% select(-Protein, -Label)))
    obj_ep <- plot(fit2, quantities = list(fontsize = 28), lty = 1:4, labels = list(fontsize = 28))
    pngname <- paste0("reptp_", str_replace(str_to_lower(bm), " ", "_"), "_", str_to_lower(meth[j]), ".png")
    dset_meth <- paste0(bm, "\n(", meth[j], ")")
    png(file = pngname, width=800, height=600)
    print(grid.arrange(grobs = list(obj_ep),
                       top = textGrob(dset_meth, x = 0.02, y = 0.2, just = "left", gp=gpar(fontsize=36,font=2))))
    dev.off()
}

# Estimation accuracy versus number of features ---------------------------

nb_ftr_lbnd <- 7
nb_ftr_hbnd <- 20

nb_ftr_lchr <- str_c("# features ", "1 - ", nb_ftr_lbnd - 1)
nb_ftr_hchr <- str_c("# features > ", nb_ftr_hbnd)
nb_ftr_mchr <- str_c("# features ", nb_ftr_lbnd, " - ", nb_ftr_hbnd)

test_res_all <- test_res_all %>% 
    mutate(nb_ftr_cat = ifelse(
        nb_ftr < nb_ftr_lbnd, 
        nb_ftr_lchr, 
        ifelse(nb_ftr > nb_ftr_hbnd, nb_ftr_hchr, nb_ftr_mchr)
    )) %>% 
    mutate(nb_ftr_cat = factor(nb_ftr_cat, levels = c(nb_ftr_lchr, nb_ftr_mchr, nb_ftr_hchr)))

est_nbftr <- test_res_all %>% 
    filter(!is.na(pvalue)) %>%
    group_by(expt, tool, method, nb_ftr_cat) %>%
    summarise(med_abserr = median(abs(log2FC - trueLog2FC), na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(dset = str_c(expt, tool, sep = ": "), acq = str_sub(expt, start = 1, end = 3))

# DIA benchmarks (Fig. S17a)
est_nbftr %>% 
    filter(!str_detect(method, "Top")) %>% 
    ggplot(aes(method, med_abserr, color = expt, shape = tool, group = dset)) + 
    geom_point(size = 3, alpha = 0.8) + 
    geom_line(aes(linetype = tool)) + 
    scale_color_manual(values = c("#56B4E9", "#CC79A7")) + 
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    facet_wrap(~ nb_ftr_cat) +
    guides(linetype = "none", color = guide_legend(order = 1), shape = guide_legend(order = 2)) + 
    labs(x = "Method", y = "Median absolute estimation error", 
         color = "Benchmark", shape = "Tool", title = "DIA benchmarks") + 
    theme_bw(base_size = 12) + 
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold")) + 
    theme(legend.position = "bottom")
# ggsave("lineerr-ftr-dia.png", width = 9, height = 4)
