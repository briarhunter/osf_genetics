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
library(ggpubr)
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
