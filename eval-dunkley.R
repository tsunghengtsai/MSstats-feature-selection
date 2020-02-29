######################################################################
## Evaluation of the MSstats analyses for the Dunkley experiment
##   - SRM Dunkley (Panorama link: https://panoramaweb.org/neuroSRM.url)
## 
## Please run first the MSstats analysis with run-msstats-dunkley.R
## 
## The results are summarized and presented in the manuscript:
##   T-H Tsai, M Choi, B Banfai, Y Liu, B MacLean, T Dunkley, and O Vitek (2019),
##   "Selection of features with consistent profiles improves relative
##   protein quantification in mass spectrometry experiments."
######################################################################

# Comparison of the manual curation and feature selection method ----------

library(tidyverse)

load("processed/processed.srm.dunkley.sl.inf.rda")

# Data processed with feature selection option
annotated_automatic <- as_tibble(processed$ProcessedData) %>%
    select(protein = PROTEIN, peptide = PEPTIDE, feature = FEATURE, run = RUN, 
           file_name = originalRUN, label = LABEL, log2inty = ABUNDANCE, 
           censored, feature_quality, is_outlier) %>% 
    mutate_if(is.factor, as.character) %>% 
    mutate(peptide = str_sub(peptide, 1, -3))

# Uninformative features detected by the feature selection method
flagged_automatic <- annotated_automatic %>% 
    filter(feature_quality == "Uninformative") %>% 
    distinct(protein, peptide, feature, label)

# Annotation by manual curation in 2 steps (available at https://panoramaweb.org/neuroSRM.url)
dunkley_full <- tbl_df(read.csv("data/SRM_Dunkley/transition_lists/Neuron SRM all features after peak review.csv", stringsAsFactors = F))
dunkley_step1 <- tbl_df(read.csv("data/SRM_Dunkley/transition_lists/Neuron SRM after manual peptide exclusion.csv", stringsAsFactors = F))
dunkley_step2 <- tbl_df(read.csv("data/SRM_Dunkley/transition_lists/Neuron SRM after manual peptide exclusion and quality filter.csv", stringsAsFactors = F))

dunkley_full <- dunkley_full %>% 
    mutate(feature = str_c(Peptide.Sequence, Precursor.Charge, Fragment.Ion, Product.Charge, sep = "_"), 
           label = ifelse(Isotope.Label.Type == "light", "L", "H")) %>% 
    select(protein = Protein.Name, peptide = Peptide.Sequence, feature, run = File.Name, 
           label, intensity = Area, Standard.Type)
dunkley_step1 <- dunkley_step1 %>% 
    mutate(feature = str_c(Peptide.Sequence, Precursor.Charge, Fragment.Ion, Product.Charge, sep = "_"), 
           label = ifelse(Isotope.Label.Type == "light", "L", "H")) %>% 
    select(protein = Protein.Name, peptide = Peptide.Sequence, feature, run = File.Name, 
           label, intensity = Area, Standard.Type, exclude, meets.quality.filter)
dunkley_step2 <- dunkley_step2 %>% 
    mutate(feature = str_c(Peptide.Sequence, Precursor.Charge, Fragment.Ion, Product.Charge, sep = "_"), 
           label = ifelse(Isotope.Label.Type == "light", "L", "H")) %>% 
    select(protein = Protein.Name, peptide = Peptide.Sequence, feature, run = File.Name, 
           label, intensity = Area, Standard.Type, exclude, meets.quality.filter)

# Features flagged as bad-quality in manual curation
flagged_step1 <- dunkley_full %>% 
    distinct(protein, peptide, feature, label) %>% 
    anti_join(dunkley_step1 %>% distinct(protein, peptide, feature, label))
flagged_step2 <- dunkley_step1 %>% 
    filter(meets.quality.filter == "NO") %>% 
    distinct(protein, peptide, feature, label)

# Exclude those proteins removed during conversion (iRT, proteins with one feature)
flagged_step1 <- flagged_step1 %>% 
    semi_join(annotated_automatic %>% distinct(protein, peptide, feature, label))
flagged_step2 <- flagged_step2 %>% 
    semi_join(annotated_automatic %>% distinct(protein, peptide, feature, label))

flagged_manual <- bind_rows(
    flagged_step1 %>% mutate(meth = "step1"), 
    flagged_step2 %>% mutate(meth = "step2")
)

# Proteins with inconsistent features detected by the proposed approach and the manual curation
flagged_automatic %>% 
    distinct(protein) %>% 
    semi_join(flagged_manual %>% distinct(protein))

# Proteins with inconsistent features detected by the proposed approach alone
flagged_automatic %>% 
    distinct(protein) %>% 
    anti_join(flagged_manual %>% distinct(protein))

# Proteins with inconsistent features detected by the manual curation alone
flagged_manual %>% 
    distinct(protein) %>% 
    anti_join(flagged_automatic %>% distinct(protein))

# Profile plots for example proteins --------------------------------------

annotated_automatic <- annotated_automatic %>% 
    mutate(run = factor(run, levels = as.character(1:n_distinct(annotated_automatic$run))))

# Color palette
cbp <- c("#009E73", "#0072B2", "#D55E00")

# Plot annotations
df_ann1 <- data.frame(
    run = rep(c("17", "46", "68", "87"), 2), 
    log2inty = rep(26, 8), 
    ann = rep(c("QC", "SA001 GE1", "SA001 GE2", "SA001"), 2), 
    label = rep(c("Reference", "Endogenous"), each = 4)
)
df_ann2 <- data.frame(
    run = rep(c("39", "45", "53", "61", "67", "75", "82", "87", "92"), 2), 
    log2inty = rep(24.5, 18), 
    ann = rep(c("D0", "D14", "D41"), 6), 
    label = rep(c("Reference", "Endogenous"), each = 9)
)

# Example protein with inconsistent features detected by both approaches (Fig. S12)
prot_name <- "sp|P63104|1433Z_HUMAN"
one_processed <- annotated_automatic %>% 
    filter(protein == prot_name)

# For annotation
df_ftr <- data.frame(
    run = c("12", "12"), 
    log2inty = c(7, 10), 
    ann = rep(c("Inconsistent\nfeatures"), 2), 
    label = c("Reference", "Endogenous")
)

one_processed %>% 
    mutate(label = ifelse(label == "H", "Reference", "Endogenous")) %>% 
    mutate(label = factor(label, levels = c("Reference", "Endogenous"))) %>% 
    ggplot(aes(run, log2inty)) + 
    geom_point(aes(shape = is_outlier, alpha = feature_quality, color = peptide)) + 
    geom_line(aes(group = feature, linetype = feature_quality, alpha = feature_quality, color = peptide)) + 
    scale_color_manual(values = cbp) + 
    scale_shape_manual(values = c(16, 1), guide = FALSE) + 
    scale_linetype_discrete(guide = FALSE) +
    scale_alpha_manual(values = c(1, 0.5), limits = c("Informative", "Uninformative"), guide = FALSE) +
    scale_x_discrete(breaks = c(35, 41, 48, 57, 63, 70, 79, 84, 88, 95)) +
    geom_vline(xintercept = c(35, 41, 48, 57, 63, 70, 79, 84, 88) + 0.5, linetype = "dashed") + 
    coord_cartesian(ylim = c(0, 26)) +
    geom_label(aes(label = ann), size = 4, label.size = NA, data = df_ann1) + 
    geom_label(aes(label = ann), size = 3.5, label.size = NA, data = df_ann2) + 
    geom_text(aes(label = ann), size = 5, color = "#D55E00", fontface = "italic", data = df_ftr) + 
    facet_grid(~ label) + 
    labs(x = "Run", y = "Log2-intensity", title = prot_name, color = "") + 
    theme_bw() + 
    theme(text = element_text(size = 16)) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold")) + 
    theme(legend.position = c(0.002, 0.005), legend.justification=c(0,0), legend.direction = "horizontal")
# ggsave("srm_both.png", width = 15, height = 5)

# Example protein with inconsistent features detected only by the proposed approach (Fig. S15a)
prot_name <- "sp|P51532|SMCA4_HUMAN"
one_processed <- annotated_automatic %>% 
    filter(protein == prot_name)

# For annotation
df_ftr <- data.frame(
    run = c("12", "12"), 
    log2inty = c(11, 8), 
    ann = rep(c("Inconsistent\nfeatures"), 2), 
    label = c("Reference", "Endogenous")
)

one_processed %>% 
    mutate(label = ifelse(label == "H", "Reference", "Endogenous")) %>% 
    mutate(label = factor(label, levels = c("Reference", "Endogenous"))) %>% 
    ggplot(aes(run, log2inty)) + 
    geom_point(aes(shape = is_outlier, alpha = feature_quality, color = peptide)) + 
    geom_line(aes(group = feature, linetype = feature_quality, alpha = feature_quality, color = peptide)) + 
    scale_color_manual(values = cbp) + 
    scale_shape_manual(values = c(16, 1), guide = FALSE) + 
    scale_linetype_discrete(guide = FALSE) +
    scale_alpha_manual(values = c(1, 0.5), limits = c("Informative", "Uninformative"), guide = FALSE) +
    scale_x_discrete(breaks = c(35, 41, 48, 57, 63, 70, 79, 84, 88, 95)) +
    geom_vline(xintercept = c(35, 41, 48, 57, 63, 70, 79, 84, 88) + 0.5, linetype = "dashed") + 
    coord_cartesian(ylim = c(0, 26)) +
    geom_label(aes(label = ann), size = 4, label.size = NA, data = df_ann1) + 
    geom_label(aes(label = ann), size = 3.5, label.size = NA, data = df_ann2) + 
    geom_text(aes(label = ann), size = 5, color = "#009E73", fontface = "italic", data = df_ftr) + 
    facet_grid(~ label) + 
    labs(x = "Run", y = "Log2-intensity", title = prot_name, color = "") + 
    theme_bw() + 
    theme(text = element_text(size = 16)) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold")) + 
    theme(legend.position = c(0.002, 0.005), legend.justification=c(0,0), legend.direction = "horizontal")
# ggsave("srm_stat.png", width = 15, height = 5)

# Example protein with features flagged as bad-quality by the manual curation (Fig. S15b)
prot_name <- "sp|Q9ULK0|GRID1_HUMAN"
one_processed <- annotated_automatic %>% 
    filter(protein == prot_name)

# For annotation
df_ftr <- data.frame(
    run = c("18", "18"), 
    log2inty = c(15, 4), 
    ann = rep(c("All flagged as bad quality\nin original study"), 2), 
    label = c("Reference", "Endogenous")
)
one_processed %>% 
    mutate(label = ifelse(label == "H", "Reference", "Endogenous")) %>% 
    mutate(label = factor(label, levels = c("Reference", "Endogenous"))) %>% 
    ggplot(aes(run, log2inty)) + 
    geom_point(aes(shape = is_outlier, alpha = feature_quality, color = peptide)) + 
    geom_line(aes(group = feature, linetype = feature_quality, alpha = feature_quality, color = peptide)) + 
    scale_color_manual(values = cbp) + 
    scale_shape_manual(values = c(16, 1), guide = FALSE) + 
    scale_linetype_discrete(guide = FALSE) +
    scale_alpha_manual(values = c(1, 0.5), limits = c("Informative", "Uninformative"), guide = FALSE) +
    scale_x_discrete(breaks = c(35, 41, 48, 57, 63, 70, 79, 84, 88, 95)) +
    geom_vline(xintercept = c(35, 41, 48, 57, 63, 70, 79, 84, 88) + 0.5, linetype = "dashed") + 
    coord_cartesian(ylim = c(0, 26)) +
    geom_label(aes(label = ann), size = 4, label.size = NA, data = df_ann1) + 
    geom_label(aes(label = ann), size = 3.5, label.size = NA, data = df_ann2) + 
    geom_text(aes(label = ann), size = 5, fontface = "italic", data = df_ftr) + 
    facet_grid(~ label) + 
    labs(x = "Run", y = "Log2-intensity", title = prot_name, color = "") + 
    theme_bw() + 
    theme(text = element_text(size = 16)) + 
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face="bold")) + 
    theme(legend.position = c(0.002, 0.005), legend.justification=c(0,0), legend.direction = "horizontal")
# ggsave("srm_man.png", width = 15, height = 5)
