######################################################################
## Ground-truth for the DDA iPRG benchmark
## 
## The ground-truth defined here is used for evaluating the iPRG 
## benchmark in eval-DDAbenchmark.R
## 
## The results are summarized and presented in the manuscript:
##   T-H Tsai, M Choi, B Banfai, Y Liu, B MacLean, T Dunkley, and O Vitek (2019),
##   "Selection of features with consistent profiles improves relative
##   protein quantification in mass spectrometry experiments."
######################################################################

# Spike-in mix
mix <- list(
    c("P44015", "sp|P44015|VAC2_YEAST"), 
    c("P55752", "sp|P55752|ISCB_YEAST"), 
    c("P44374", "sp|P44374|SFG2_YEAST"), 
    c("P44983", "sp|P44983|UTR6_YEAST"), 
    c("P44683", "sp|P44683|PGA4_YEAST"), 
    c("P55249", "sp|P55249|ZRT4_YEAST")
)

conc <- list(
    c(65, 55, 15, 2), 
    c(55, 15, 2, 65), 
    c(15, 2, 65, 55), 
    c(2, 65, 55, 15), 
    c(11, 0.6, 10, 500), 
    c(10, 500, 11, 0.6)
)

# Relative concentration across conditions of the mix
cond_code <- "C"
nb_cond <- length(conc[[1]])
lab <- rep(NA, (nb_cond ^ 2 - nb_cond) / 2)
fc <- vector("list", length(conc))
for (m in 1:length(fc)) {
    fc[[m]] <- rep(NA, (nb_cond ^ 2 - nb_cond) / 2)
}

idx <- 1
for (i in 1:(nb_cond - 1)) {
    for (j in (i + 1):nb_cond) {
        lab[idx] <- paste0(cond_code, j, "-", cond_code, i)
        for (m in 1:length(fc)) {
            fc[[m]][idx] <- conc[[m]][j] / conc[[m]][i]
        }
        idx <- idx + 1
    }
}

# Ground-truth for performance evaluation
gold <- tibble(
    Protein = rep(unlist(mix), each = length(lab)), 
    Label = rep(lab, sum(map_int(mix, length))), 
    trueFC = unlist(map2(fc, map(mix, length), rep))
)
