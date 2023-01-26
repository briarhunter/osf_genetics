#### For Fst values ####

library(fst)
library(tidyfst)
library(tidyverse)

getwd()
# read the output

het_imputed = read.table("imputed_het_without_OSF_315.het", header = T)

names(het_imputed)

# [1] "INDV"    "O.HOM."  "E.HOM."  "N_SITES" "F"

names(het_imputed)= c("names_ind", "nb_homO", "nb_homE", "nb_sites", "F")

het_imputed$nb_hetO = het_imputed$nb_sites -het_imputed$nb_homO

het_imputed$hetO = het_imputed$nb_hetO/het_imputed$nb_sites

view(het_imputed)
summary(het_imputed) # 321

#### bring in frog info
frog_info = read.table("oregon_spotted_322_individuals.tsv", sep = "\t", header = T)
names(frog_info)
view(frog_info) # 322 individuals... need to remove OSF_315
### already changed OSF_273 and OSF_274 to be GVZoo - captive frogs (not wild / Maria)
frog_info <- frog_info %>% 
  filter(names_ind != 'OSFrog_315A')
summary(frog_info) # down to 321 individuals 

info = frog_info[,c(1,7:13)]
df <- merge(info, het_imputed, by = "names_ind")
view(df)

##### now have merged info for frogs along with calculated diversity indices
##### plot and analyze significance between populations

##############################################################################
#### 1. Compare F (inbreeding coefficients) between source populations #######
##############################################################################
library(ggplot2)
library(ggpubr)
hist(df$F)
ggqqplot(df$F) # looks normal
library(rstatix)
df %>% 
  group_by(source) %>% 
  shapiro_test(F) # VanAqua is only one not normal (p < 0.05)

# Null hypothesis of bartlett: variance is equal
boxplot(F ~ source, data = df) # this does not work because infinite xlim values
plot(df$F)
bartlett.test(F ~ source, data = df)
# p-value: 0.08804 - accept the null hypothesis. Variance is homogenous - proceed ANOVA
F.aov <- aov(F ~ source, data = df)
summary(F.aov) # p < 0.05; significant differences present 
TukeyHSD(F.aov) # post-hoc to compare all groups
# p < 0.05: Maria:Elk, Morris:Elk, VA:Elk, (almost ST:Morris)

## get some summary measurements
df %>%
  group_by(source) %>%
  summarise_at(vars(F), list(name = mean))

library(EnvStats) # for adding sample size
ggplot(data = df, aes(x = source, y = F, fill=source))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Population", y = "Inbreeding Coefficient (F)")+
  theme_classic()+
  scale_fill_brewer(palette = "Set2")+
  scale_x_discrete(labels = c("Elk Brook" = "Elk Brook", "Greater Vancouver Zoo" = "GVZoo",
                                 "Maria Slough" = "Maria Slough", "Morris Valley" = "Morris Valley",
                                 "Mountain Slough" = "Mountain", "Semmihault" = "Semmihault",
                                 "Toronto Zoo" = "TZoo", "Vancouver Aquarium" = "VanAqua"))+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = "Arial"))+
  theme(text = element_text(size = 12))+
  stat_n_text() #sample size - figure out how to put it in the actual box

### plot
plot <- ggboxplot(df, x = "source", y = "F")
compare_means(F ~ source, data = df)
my_comparisons <- list( c("Toronto Zoo", "Semmihault"), c("Toronto Zoo", "Elk Brook"),
                        c("Vancouver Aquarium", "Semmihault"), c("Vancouver Aquarium", "Elk Brook"),
                        c("Greater Vancouver Zoo", "Semmihault"), c("Greater Vancouver Zoo", "Elk Brook"),
                        c("Maria Slough", "Semmihault"), c("Maria Slough", "Elk Brook"), c("Morris Valley", "Semmihault"),
                        c("Morris Valley", "Elk Brook"))
ggboxplot(df, x = "source", y = "F", 
          color = "source", palette = "Set2")+
  stat_compare_means(comparisons = my_comparisons)

##############################################################################
#### 2. Compare Ho (observed heterozygosity) between source populations #######
##############################################################################
hist(df$hetO)
ggqqplot(df$hetO) # looks normal
df %>% 
  group_by(source) %>% 
  shapiro_test(hetO) # VanAqua is only one not normal (p < 0.05)

# Null hypothesis of bartlett: variance is equal
boxplot(hetO ~ source, data = df) # a couple potential outliers but none look major
bartlett.test(hetO ~ source, data = df)
# p-value: 0.08804 - accept the null hypothesis. Variance is homogenous - proceed ANOVA
Ho.aov <- aov(hetO ~ source, data = df)
summary(Ho.aov) # p < 0.05; significant differences present 
TukeyHSD(Ho.aov) # post-hoc to compare all groups
# p < 0.05: Elk - GVZ, MS, MV, (almost TZ), VA. (almost ST:Morris)

df %>%
  group_by(source) %>%
  summarise_at(vars(hetO), list(name = mean))

ggplot(data = df, aes(x = source, y = hetO, fill=source))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Population", y = "Heterozygosity (Ho)")+
  theme_classic()+
  scale_fill_brewer(palette = "Set2")+
  scale_x_discrete(labels = c("Elk Brook" = "Elk Brook", "Greater Vancouver Zoo" = "GVZoo",
                              "Maria Slough" = "Maria Slough", "Morris Valley" = "Morris Valley",
                              "Mountain Slough" = "Mountain", "Semmihault" = "Semmihault",
                              "Toronto Zoo" = "TZoo", "Vancouver Aquarium" = "VanAqua"))+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = "Arial", size = 12))+
  stat_n_text() #sample size - figure out how to put it in the actual box

## add zoo vs wild column and compare Ho / F
df <- df %>% 
  mutate(population_type = case_when(source == "Maria Slough" ~ "wild",
                                source == "Morris Valley" ~ "wild",
                                source == "Mountain Slough" ~ "wild",
                                source == "Semmihault" ~ "wild",
                                source == "Elk Brook" ~ "wild",
                                source == "Toronto Zoo" ~ "zoo",
                                source == "Vancouver Aquarium" ~ "zoo",
                                source == "Greater Vancouver Zoo" ~ "zoo")) %>% 
  relocate(population_type, .after = source)

ggplot(data = df, aes(x = population_type, y = hetO, colour=population_type))+
  geom_boxplot(show.legend = FALSE)+
  labs(x = "Population", y = "Heterozygosity (Ho)")+
  theme_classic()+
  theme(legend.title = element_blank())+
  theme(text = element_text(family = "Arial", size = 12))+
  stat_n_text()

## test significance
df %>% 
  group_by(population_type) %>% 
  shapiro_test(hetO) # zoo is not normal (because of VanAqua?) (p < 0.05)

# Null hypothesis of bartlett: variance is equal
boxplot(hetO ~ population_type, data = df) # 
bartlett.test(hetO ~ population_type, data = df)
# p<0.05, reject null. not homogenous. Cannot use ANOVA

kruskal.test(hetO ~ population_type, data = df)
# chi-squared = 0.568, p-value = 0.451 (accept the null) 
####### There are no significant differences in mean Ho between zoo and wild


#######################################################
##################### Fst #############################
#######################################################

fst.EMV = read.table("Fst/fst_elk_brook_vs_morris_valley.weir.fst", header = T, na = "NA")
view(fst.EMV)
str(fst.EMV)
summary(fst.EMV$WEIR_AND_COCKERHAM_FST)
fst.EMS = read.table("Fst/fst_elk_brook_vs_maria_slough_minus_three_ind.weir.fst", header = T)
view(fst.EMS)
summary(fst.EMS$WEIR_AND_COCKERHAM_FST) # updated after removing 3 individuals
fst.SV = read.table("Fst/fst_maria_slough_minus_three_ind_vs_morris_valley.weir.fst", header = T)
summary(fst.SV$WEIR_AND_COCKERHAM_FST)
fst.EG = read.table("Fst/fst_elk_brook_vs_great_vancouver_zoo.weir.fst", header = T)
summary(fst.EG$WEIR_AND_COCKERHAM_FST)
fst.ET = read.table("Fst/fst_elk_brook_vs_toronto_zoo.weir.fst", header = T)
summary(fst.ET$WEIR_AND_COCKERHAM_FST)
fst.EVA = read.table("Fst/fst_elk_brook_vs_vancouver_aquarium.weir.fst", header = T)
summary(fst.EVA$WEIR_AND_COCKERHAM_FST)
fst.MSG = read.table("Fst/fst_maria_slough_vs_great_vancouver_zoo.weir.fst", header = T)
summary(fst.MSG$WEIR_AND_COCKERHAM_FST)
fst.MST = read.table("Fst/fst_maria_slough_vs_toronto_zoo.weir.fst", header = T)
summary(fst.MST$WEIR_AND_COCKERHAM_FST)
fst.MSV = read.table("Fst/fst_maria_slough_vs_vancouver_aquarium.weir.fst", header = T)
summary(fst.MSV$WEIR_AND_COCKERHAM_FST)
fst.MVG = read.table("Fst/fst_morris_valley_vs_great_vancouver_zoo.weir.fst", header = T)
summary(fst.MVG$WEIR_AND_COCKERHAM_FST)
fst.MVT = read.table("Fst/fst_morris_valley_vs_toronto_zoo.weir.fst", header = T)
summary(fst.MVT$WEIR_AND_COCKERHAM_FST)
fst.MVV = read.table("Fst/fst_morris_valley_vs_vancouver_aquarium.weir.fst", header = T)
summary(fst.MVV$WEIR_AND_COCKERHAM_FST)
## need to find zoo vs zoo fst files
fst.GV = read.table("Fst/fst_greater_zoo_vancouver_vs_vancouver_aquarium.weir.fst", header = T)
summary(fst.GV$WEIR_AND_COCKERHAM_FST)
fst.TG = read.table("Fst/fst_toronto_zoo_vs_greater_vancouver_zoo.weir.fst", header = T)
summary(fst.TG$WEIR_AND_COCKERHAM_FST)
fst.TV = read.table("Fst/fst_toronto_zoo_vs_vancouver_aquarium.weir.fst", header = T)
summary(fst.TV$WEIR_AND_COCKERHAM_FST)
#######################################################################
####### Test for significant differences in these mean Fst values #####

### Need to combine all Fst values first
EG <- fst.EG %>% 
  rename(EK_GVZ = WEIR_AND_COCKERHAM_FST)
EMS <- fst.EMS %>% 
  rename(EK_MS = WEIR_AND_COCKERHAM_FST)
EMV <- fst.EMV %>% 
  rename(EK_MV = WEIR_AND_COCKERHAM_FST)
ET <- fst.ET %>% 
  rename(EK_TZ = WEIR_AND_COCKERHAM_FST)
EVA <- fst.EVA %>% 
  rename(EK_VA = WEIR_AND_COCKERHAM_FST)
GV <- fst.GV %>% 
  rename(GVZ_VA = WEIR_AND_COCKERHAM_FST)
MSG <- fst.MSG %>% 
  rename(MS_GVZ = WEIR_AND_COCKERHAM_FST)
MST <- fst.MST %>% 
  rename(MS_TZ = WEIR_AND_COCKERHAM_FST)
MSV <- fst.MSV %>% 
  rename(MS_VA = WEIR_AND_COCKERHAM_FST)
MVG <- fst.MVG %>% 
  rename(MV_GVZ = WEIR_AND_COCKERHAM_FST)
MVT <- fst.MVT %>% 
  rename(MV_TZ = WEIR_AND_COCKERHAM_FST)
MVV <- fst.MVV %>% 
  rename(MV_VA = WEIR_AND_COCKERHAM_FST)
TG <- fst.TG %>% 
  rename(TZ_GVZ = WEIR_AND_COCKERHAM_FST)
TV <- fst.TV %>% 
  rename(TZ_VA = WEIR_AND_COCKERHAM_FST)
SV <- fst.SV %>% 
  rename(MS_MV = WEIR_AND_COCKERHAM_FST)
FST <- merge(EG, EMS, EMV, ET, EVA, GV, MSG, MST, MSV, MVG, MVT, MVV, SV, TG, TV, by=c("CHROM", "POS"))
FST <- merge(EG, EMS, by=c("CHROM", "POS"))
FST <- merge(FST, EMV, by=c("CHROM", "POS"))
FST <- merge(FST, ET, by=c("CHROM", "POS"))
FST <- merge(FST, EVA, by=c("CHROM", "POS"))
FST <- merge(FST, GV, by=c("CHROM", "POS"))
FST <- merge(FST, MSG, by=c("CHROM", "POS"))
FST <- merge(FST, MST, by=c("CHROM", "POS"))
FST <- merge(FST, MSV, by=c("CHROM", "POS"))
FST <- merge(FST, MVG, by=c("CHROM", "POS"))
FST <- merge(FST, MVT, by=c("CHROM", "POS"))
FST <- merge(FST, MVV, by=c("CHROM", "POS"))
FST <- merge(FST, SV, by=c("CHROM", "POS"))
FST <- merge(FST, TG, by=c("CHROM", "POS"))
FST <- merge(FST, TV, by=c("CHROM", "POS"))
view(FST)
#### Now have all Fst values combined - make it long not wide
library(reshape2)
FST_l <- melt(setDT(FST), id = 1:2, value.name = "fst", variable.name = "pops", na.rm = F) # keep NaN
dim(FST_l) # 333,435 rows (22,229 x 15)
############################################
#### Test for significant differences in Fst between populations
# hist(FST_l$fst) # right tail
# ggqqplot(FST_l$fst) # definitely not normal
# 
# boxplot(fst ~ pops, data = FST_l)
# bartlett.test(fst ~ pops, data = FST_l) # p < 0.05, variance is NOT homogenous
# kruskal.test(fst ~ pops, data = FST_l) # p < 0.05, there are significant differences
# require(FSA)
# dunnTest(fst ~ pops, data = FST_l, method = "holm")

###############################################################################
###############################################################################
#### Cannot test for significance because measurements are not independent ####
###############################################################################
# Would have to use permutation test - but this is not common practice anyway #