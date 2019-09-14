######################################################################
## Ground-truth for the DDA Choi benchmark
## 
## The ground-truth defined here is used for evaluating the Choi 
## benchmark in eval-DDAbenchmark.R
## 
## The results are summarized and presented in the manuscript:
##   T-H Tsai, M Choi, B Banfai, Y Liu, T Dunkley, and O Vitek (2019),
##   "Selection of features with consistent profiles improves relative
##   protein quantification in mass spectrometry experiments."
######################################################################

# Spike-in mix
mix <- list(
    c("SS1_P00432", "SS1_P00711", "SS1_P02701", "SS1_P02754", "SS1_P0CG53", 
      "SS1_P24627", "SS1_P68082", "SS1_P80025", "SS1_Q29443"), 
    c("SS2_P00915", "SS2_P00921", "SS2_P01008", "SS2_P01012", "SS2_P02662",
      "SS2_P02663", "SS2_P02666", "SS2_P02787", "SS2_P05307", "SS2_P61769"), 
    c("SS3_P00563", "SS3_P00698", "SS3_P02769", "SS3_Q3SX14", "SS3_Q58D62",
      "SS4_P00004", "SS4_P00442", "SS4_P01133", "SS4_P02753")
)

# Relative concentration across conditions of the mix
conc <- list(
    2 ^ c(0, -1, -2, -1, 0), 
    2 ^ c(0, 2, 1, 2, 1), 
    2 ^ c(0, -1, 1, 0, 1)
)

cond_code <- "M"
nb_cond <- length(conc[[1]])
lab <- rep(NA, (nb_cond ^ 2 - nb_cond) / 2)
fc <- vector("list", length(conc))
for (m in 1:length(fc)) {
    fc[[m]] <- rep(NA, (nb_cond ^ 2 - nb_cond) / 2)
}

idx <- 1
for (i in 1:(nb_cond - 1)) {
    for (j in (i + 1):nb_cond) {
        lab[idx] <- paste0(cond_code, i, "-", cond_code, j)
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
