#######Physiology#######
#physiology data and analysis for SW08 and WH
#SW08: Mecklenburg Bight (Western Baltic Sea), WH: Wilhelmshaven (North Sea)

library(tidyverse)
library(car)
library(arm)
library(MASS)
library(emmeans)
library(readxl)
library(lme4)
library(olsrr)

#setwd("~/Documents/Experiments/Physiology_SW08-WH/")
setwd("~/SalAdapt_transcriptomics/data/")

####acute survival####
dat_acute <- read.csv("Acute_sal.csv", sep = ";")

dat_acute$Time_factor <- as.factor(dat_acute$Time)
dat_acute$Station <- factor(dat_acute$Station, 
                              levels = c("Wilhelmshaven", "SW08"))
dat_acute$Treatment_factor <- as.factor(dat_acute$Treatment_sal)


#including Replicate as random factor does not improve model fit, no significance

#m <- bayesglm(Active ~ Station + Time + Treatment_sal, family = binomial(), data = dat_acute)

m <- glm(Active ~ Station*Time*Treatment_factor , family = binomial(), data = dat_acute)#Replicate not significant
summary(m)
anova(m, test = "Chisq")

#only include data below 15 PSU in plot
log_data_24 <- dat_acute %>%
  filter(Treatment_sal < 15)

#plot
ggplot(data = log_data_24, aes(x =Treatment_sal, y = Active))+
  theme_light(base_size = 14)+
  geom_point(aes(fill = Station), size = 2,shape = 21, alpha = 0.2,
             position = position_dodge(0.15))+
  geom_smooth(aes(group = interaction(Station, Time_factor), 
                  fill = Station, col = Station, linetype = Time_factor),
              method = "glm", method.args = list(family = "binomial"), alpha = 0.02)+
  scale_fill_manual(values=c("#DE3C22", "#83B3C2"))+
  scale_color_manual(values =c("#DE3C22", "#83B3C2"))+
  scale_linetype_manual(values = c("22", "solid")) +
  theme(legend.position = "bottom")+
  xlab("Salinity [PSU]")+ ylab("Proportion active")



#calculate LD50
get_LD50 = function(fit){
  data.frame(
    LD50 = dose.p(fit)[1],
    SE = attributes(dose.p(fit))$SE[,1]
  )
}

#LD50 <- dat_acute %>% group_by(Station, Original_sal, Time) %>% 
  #do(get_LD50(glm(Active ~ Treatment_sal, family = "binomial", data = .)))

LD50 <- dat_acute %>% group_by(Station, Original_sal, Time, Replicate) %>% 
  do(get_LD50(glm(Active ~ Treatment_sal, family = "binomial", data = .)))

LD50_min_max <- LD50 %>%
  group_by(Station, Time)%>%
  summarise(Min = min(LD50), Max = max(LD50))

LD50_mean <- LD50 %>%
  group_by(Station, Original_sal, Time)%>%
  summarise(Mean = mean(LD50), SD = sd(LD50))%>%
  left_join(LD50_min_max, by = c("Station", "Time"))

write.csv(LD50, "LD50.csv", row.names = FALSE)

LD50_mean$Time_factor <- as.factor(LD50_mean$Time)

ggplot(data = LD50_mean, aes(x = Time_factor, y = Mean)) +
  theme_light(base_size = 14) +
  geom_errorbar(aes(ymin = Min, ymax = Max, group = Station), 
                width = 0.05, linewidth = 0.7, 
                position = position_dodge(width = 0.5)) +  
  geom_point(aes(fill = Station), size = 12, shape = 21, 
             position = position_dodge(width = 0.5)) +   
  scale_fill_manual(values = c("#DE3C22", "#83B3C2")) +
  scale_y_reverse() +
  ylab(expression("LS"["50"]*" [PSU]")) +
  xlab("Time [h]")


####survival to adulthood#####
dat_survival <- read.csv("Survival.csv", sep = ";", header = TRUE)

#format and calculate variables
dat_survival$Fraction_Survived <- dat_survival$survived/dat_survival$Nauplii_t0
dat_survival$Fraction_Nauplii <- dat_survival$Nauplii/dat_survival$Nauplii_t0
dat_survival$Fraction_Copepodites_1.3 <- dat_survival$Copepodites_1.3/dat_survival$Nauplii_t0
dat_survival$Fraction_Copepodites_4.5 <- dat_survival$Copepodites_4.5/dat_survival$Nauplii_t0
dat_survival$Fraction_Adults <- dat_survival$Adults/dat_survival$Nauplii_t0
dat_survival$Fraction_Dead <- (20-(dat_survival$Nauplii+dat_survival$Copepodites_1.3+dat_survival$Copepodites_4.5+dat_survival$Adults))/dat_survival$Nauplii_t0
dat_survival$Day_factor <- as.factor(dat_survival$Day)
dat_survival$Treatment_factor <- as.factor(dat_survival$Treatment_sal)
dat_survival$GroupID <- paste(dat_survival$Station, dat_survival$Treatment_sal, dat_survival$Replicate, sep = "_")

#box plot
ggplot(dat_survival, aes(x = Day_factor, y = Fraction_Survived , fill = Station))+
  facet_wrap(~Treatment_factor, labeller = labeller(Treatment_factor = c("7" = "7 PSU", "15" = "15 PSU")))+
  geom_boxplot(alpha = 0.5)+
  theme_light(base_size = 14)+
  scale_fill_manual(values = c( "#83B3C2", "#DE3C22"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 12))+
  ylab("Fraction alive") + xlab("Time [days]")

#line plot
#mean and SE table
sum_survival <- dat_survival %>%
  filter(Day != "NA")%>%
  group_by(Day, Station, Treatment_factor) %>%
  dplyr::summarise(
    Mean_Survival = mean(Fraction_Survived),
    SD = sd(Fraction_Survived)
  )

#plot lines

ggplot(sum_survival, aes(x = Day, y = Mean_Survival)) +
  geom_point(aes(col = Station), size = 3)+
  geom_line(aes(col = Station, linetype = Treatment_factor)) +
  geom_errorbar(aes(ymin = Mean_Survival - SD, ymax = Mean_Survival + SD,
                    col = Station), width = 0.2) +
  labs(x = "Day", y = "Mean Survival") +
  theme_light()+
  scale_color_manual(values = c( "#83B3C2", "#DE3C22"))

#female : male

dat_survival %>%
  filter(Day == 21 & Treatment_factor == "15")%>%
  group_by(Day, Station, Treatment_factor) %>%
  dplyr::summarise(
    Females = mean(females),
    Males = mean(males),
    Ratio = Females/Males
  )


#stats
hist(dat_survival$Fraction_Survived, breaks = 15)#binomial distribution

stat_survival <- glmer(survived / Nauplii_t0 ~ Station + Treatment_factor * Day + (1 | GroupID), 
               data = dat_survival, family = binomial, weights = Nauplii_t0)
#this model has lower AIC than full interaction model
summary(stat_survival)


####egg production and hatching####
dat_hatch <- read.csv("Hatching.csv", sep = ";")%>%
  group_by(Station, Treatment_sal)

dat_hatch$Station <- factor(dat_hatch$Station, levels = c("WH", "SW08"))
dat_hatch$Treatment_factor <- as.factor(dat_hatch$Treatment_sal)

#add colum for total eggs
dat_hatch$Eggs_total <- dat_hatch$Eggs_unhatched + dat_hatch$Nauplii

dat_hatch %>%
  group_by(Station, Treatment_sal) %>%
  dplyr::summarise(
     Eggs = mean(Eggs_total),
     SD_Eggs = sd(Eggs_total),
     Nau = mean(Nauplii),
     SD_Nau = sd(Nauplii),
     per_hatch = Nau/Eggs*100)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(Station)%>%
  mutate(mean_hatch = mean(per_hatch),
         mean_Nau = mean(Nau),
         mean_SD_Nau = mean(SD_Nau),
         mean_Eggs = mean(Eggs),
         mean_SD_Eggs = mean(SD_Eggs))

#total eggs
ggplot(dat_hatch, aes(x = Station, y = Eggs_total, fill = Station, alpha = Treatment_factor))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("#DE3C22", "#83B3C2"))+
  theme_light()+
  scale_alpha_discrete(range=c(0.2, 0.85))

#hatched eggs (= nauplii)
ggplot(dat_hatch, aes(x = Station, y = Nauplii, fill = Station, alpha = Treatment_factor))+
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("#DE3C22", "#83B3C2"))+
  theme_light()+
  scale_alpha_discrete(range=c(0.2, 0.85))

#stats
hist(dat_hatch$Eggs_total)
stat_eggs <- lm(Eggs_total ~ Station*Treatment_sal, data = dat_hatch)

summary(stat_eggs)#not sig
ols_test_normality(stat_eggs)#Kolmogorov-Smirnov for observation > 50, not sig = normal distribution
ols_test_correlation(stat_eggs)#correlation between observed and expected residuals (under normality)

par(mfrow = c(2,2))
plot(stat_eggs)
par(mfrow = c(1,1))

durbinWatsonTest(stat_eggs)#no autocorrelation if p > 0.05

hist(dat_hatch$Nauplii)
stat_hatch <- lm(Nauplii ~ Station*Treatment_sal, data = dat_hatch)

summary(stat_hatch)#not sig
ols_test_normality(stat_hatch)#Kolmogorov-Smirnov for observation > 50, not sig = normal distribution
ols_test_correlation(stat_hatch)#correlation between observed and expected residuals (under normality)

par(mfrow = c(2,2))
plot(stat_hatch)
par(mfrow = c(1,1))

durbinWatsonTest(stat_hatch)#no autocorrelation if p > 0.05

####length - conversion factors ####
#dat_length <- read_excel("~/Downloads/Acartia_ImageJ_WH_SW08.xlsx")
dat_length <- read.csv("Length_all.csv", sep = ";")

#calculate mean length per animal 
dat_length<- dat_length %>%
  group_by(Animal_ID, Month) %>%
  dplyr::mutate(Mean_Length = mean(Length))

#remove unused columns and filter out unusable values (duplicates and bad qaulity)
dat_length_m <- dat_length %>%
  filter(Measurement_quality != "bad" |  is.na(Measurement_quality))

#keep only one row per animal  
dat_length_m <- dat_length_m %>%
  filter(File_name != "") 

ggplot(dat_length_m, aes(x = Station, y = Mean_Length, col = Station))+
  geom_boxplot()+
  scale_color_manual(values = my_colors)+
  geom_point(alpha = 0.2)+
  facet_wrap(~Treatment)+
  theme_light()

#calculate mean length per station and treatment to standardize
dat_length_m <- dat_length_m %>%
  group_by(Station, Treatment_PSU, Month) %>%
  mutate(Mean_Group = mean(Mean_Length))

# calculate dry weight
# Lugol factor 1.205 (to account for shrinkage in Lugol, Jaspers & Carstensen, 2009)
# W = 13.4 x Length^3 (Kiørboe et al. 1985)

dat_length_m$Weight <- (((dat_length_m$Mean_Group * 1.205)*0.001)^3)*13.4

#calculate conversion factors
conversion_factors <- dat_length_m %>%
  dplyr::summarise(Tank_ID, Weight)%>%
  dplyr::distinct(Tank_ID, .keep_all = TRUE) %>%  # Keep one row per Tank_ID
  dplyr::select(Tank_ID, Weight)  

conversion_factors$Month <- as.numeric(conversion_factors$Month)

conversion_factors <- conversion_factors %>%
  mutate(Group = paste(Station, Treatment_PSU, Month, sep = "_"))

####respiration####
#overnight acclimation
dat_resp_long <- read.csv("Respiration_long.csv", sep = ";")

#calculate variables
dat_resp_long$O2_final <- dat_resp_long$Oxygen_umol.h/dat_resp_long$Individuals
dat_resp_long$Salinity <- as.factor(dat_resp_long$Salinity)
dat_resp_long$Station <- factor(dat_resp_long$Station, levels = c("WH","SW08"))
dat_resp_long <- dat_resp_long %>%
  mutate(Group = paste(Station, Salinity, Month, sep = "_"))

dat_resp_long <- dat_resp_long %>%
  left_join(conversion_factors, by = "Group") %>%
  mutate(Oxygen_adjusted = (O2_final/Weight)*1000)#convert to mg body weight instead of ug

#per dry body weight (in ug)
#SW08 - 15: 13.08976
#SW08 - 5: 12.91851
#WH - 15: 13.29447
#WH - 5: 13.57943

#umol O2
ggplot(dat_resp_long, aes(x = factor(Station.x), y = O2_final, fill = factor(Station.x))) +
  theme_light()+
  facet_wrap(~Salinity, labeller = labeller(Salinity = c("5" = "5 PSU", "15" = "15 PSU")))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 12))+
  labs(x = "Station", y = "umol O2 h-1 ind-1") +
  geom_boxplot(alpha = 0.2, outlier.shape = NA)+
  geom_point(aes(color = factor(Station.x)), position = position_jitter(width = 0.15), size = 2)+
  scale_fill_manual(values = c("#DE3C22", "#83B3C2"), guide = "none") +  # Specify blue and orange colors for fill without legend
  scale_color_manual(values = c("#DE3C22", "#83B3C2"))

#adjusted for dry body weight
ggplot(dat_resp_long, aes(x = factor(Station.x), y = Oxygen_adjusted, fill = factor(Station.x))) +
  theme_light()+
  facet_wrap(~Salinity, labeller = labeller(Salinity = c("5" = "5 PSU", "15" = "15 PSU")))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 12))+
  labs(x = "Station", y = "umol O2 h-1 mg dry weight-1") +
  geom_boxplot(alpha = 0.2, outlier.shape = NA)+
  geom_point(aes(color = factor(Station.x)), position = position_jitter(width = 0.15), size = 2)+
  scale_fill_manual(values = c("#DE3C22", "#83B3C2"), guide = "none") +  # Specify blue and orange colors for fill without legend
  scale_color_manual(values = c("#DE3C22", "#83B3C2"))

#summary and stats
hist(dat_resp_long$Oxygen_adjusted, breaks = 15)
dat_resp_long %>%
  group_by(Salinity, Station.x) %>%
  summarise(Mean = mean(Oxygen_adjusted), SD = sd(Oxygen_adjusted)) 

stat_resp_long <- aov(Oxygen_adjusted ~ Station.x*Salinity, data = dat_resp_long)
summary(stat_resp_long)# sig diff for station not salinity

#check residuals
par(mfrow = c(2,2))
plot(stat_resp_long)
par(mfrow = c(1,1))

#short acclimation
dat_resp_short <- read.csv("Respiration_short.csv", sep = ";")

#calculate variables
dat_resp_short$O2_final <- dat_resp_short$Oxygen_umol.h/dat_resp_short$Individuals
dat_resp_short$Treatment_PSU <- as.factor(dat_resp_short$Treatment_PSU)
dat_resp_short$Station <- factor(dat_resp_short$Station, levels = c("WH","SW08"))
dat_resp_short <- dat_resp_short %>%
  mutate(Group = paste(Station, Treatment_PSU, Month, sep = "_"))

dat_resp_short <- dat_resp_short %>%
  left_join(conversion_factors, by = "Group") %>%
  dplyr::mutate(Oxygen_adjusted = (O2_final/Weight)*1000)

#umol O2
ggplot(dat_resp_short, aes(x = factor(Treatment_PSU.x), y = O2_final, fill = factor(Station.x))) +
  theme_light()+
  facet_wrap(~Station.x)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 12))+
  labs(x = "Salinity [PSU]", y = "umol o2 h-1 ind-1") +
  geom_boxplot(alpha = 0.2, outlier.shape = NA)+
  geom_point(aes(color = factor(Station.x)), position = position_jitter(width = 0.15), size = 2)+
  scale_fill_manual(values = c("#DE3C22", "#83B3C2"), guide = "none") +  # Specify blue and orange colors for fill without legend
  scale_color_manual(values = c("#DE3C22", "#83B3C2"))

#adjusted for dry body weight
ggplot(dat_resp_short, aes(x = factor(Treatment_PSU.x), y = Oxygen_adjusted, fill = factor(Station.x))) +
  theme_light()+
  facet_wrap(~Station.x)+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black', size = 12))+
  labs(x = "Salinity [PSU]", y = "umol o2 h-1 mg dry weight-1") +
  geom_boxplot(alpha = 0.2, outlier.shape = NA)+
  geom_point(aes(color = factor(Station.x)), position = position_jitter(width = 0.15), size = 2)+
  scale_fill_manual(values = c("#DE3C22", "#83B3C2"), guide = "none") +  # Specify blue and orange colors for fill without legend
  scale_color_manual(values = c("#DE3C22", "#83B3C2"))

#summary and stats
hist(dat_resp_short$Oxygen_adjusted, breaks = 15)
dat_resp_short %>%
  group_by(Treatment_PSU.x, Station.x) %>%
  summarise(Mean = mean(Oxygen_adjusted), SD = sd(Oxygen_adjusted)) 

stat_resp_short <- aov(Oxygen_adjusted ~ Station.x*Treatment_PSU.x, data = dat_resp_short)
summary(stat_resp_short)# sig diff for salinity not station

#check residuals
par(mfrow = c(2,2))
plot(stat_resp_short)
par(mfrow = c(1,1))




