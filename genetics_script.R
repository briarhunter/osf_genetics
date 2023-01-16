#### For Fst values ####

library(fst)
library(tidyfst)

# write_fst("Fst/fst_elk_brook_vs_morris_valley.weir.fst", compress = 50, uniform_encoding = TRUE)
# fst <- import_fst("Fst/fst_elk_brook_vs_morris_valley.weir.fst")
# 
# summary("fst_elk_brook_")
####### Cannot get above code to work because FST files are apparently in old format #######


getwd()
# read the output

het_imputed = read.table("imputed_het_without_OSF_315.het", header = T)

names(het_imputed)

# [1] "INDV"    "O.HOM."  "E.HOM."  "N_SITES" "F"

names(het_imputed)= c("names_ind", "nb_homO", "nb_homE", "nb_sites", "F")

het_imputed$nb_hetO = het_imputed$nb_sites -het_imputed$nb_homO

het_imputed$hetO = het_imputed$nb_hetO/het_imputed$nb_sites

view(het_imputed)
