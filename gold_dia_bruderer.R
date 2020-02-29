######################################################################
## Ground-truth for the DIA Bruderer benchmark
## 
## The ground-truth defined here is used for evaluating the Bruderer 
## benchmark in eval-DIAbenchmark.R
## 
## The results are summarized and presented in the manuscript:
##   T-H Tsai, M Choi, B Banfai, Y Liu, B MacLean, T Dunkley, and O Vitek (2019),
##   "Selection of features with consistent profiles improves relative
##   protein quantification in mass spectrometry experiments."
######################################################################

# Spike-in mix
mix <- list(
    c("P02754", "P00921", "P80025", "P02662", "P00366", "sp|P02754|LACB_BOVIN", 
      "sp|P00921|CAH2_BOVIN", "sp|P80025|PERL_BOVIN", "sp|P02662|CASA1_BOVIN", 
      "sp|P00366|DHE3_BOVIN"), 
    c("P12799", "P02672", "P02789", "P02676", "P61823", "sp|P12799|FIBG_BOVIN", 
      "sp|P02672|FIBA_BOVIN", "sp|P02789|TRFE_CHICK", "sp|P02676|FIBB_BOVIN", 
      "sp|P61823|RNAS1_BOVIN"), 
    c("P68082", "P02666", "sp|P68082|MYG_HORSE", "sp|P02666|CASB_BOVIN")
)

# Relative concentration across conditions of the mix
conc <- list(
    c(1.5, 1.65, 1.815, 1.995, 15, 16.515, 18.165), 
    c(100, 62.995, 39.685, 25, 2, 1.26, 0.795),
    c(0.05, 0.2, 0.8, 3.2, 12.8, 51.2, 204.8)
)

cond_code <- "S"
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
