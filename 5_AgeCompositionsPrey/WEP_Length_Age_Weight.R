# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)
  # Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# Citation: 
# Dorn MW and CL Barnes. In revision. Time-varying predation as a modifier of constant natural mortality for Gulf of Alaska walleye pollock. Fish Res. VSI Natural Mortality. Forthcoming.

# This script includes code necessary to identify age compositions of pollock consumed by major pollock predators in the Gulf of Alaska; specifically,  Arrowtooth Flounder, Pacific Cod, Pacific Halibut, Sablefish, and Walleye Pollock. Methods are described in Barnes et al. (2020). Age compositions were used to estimate changes age-specific predation of pollock through time. 

# First, we used survey data to estimate von Bertalanffy (1938) growth parameters to predict pollock age from length. We then quantified bias-corrected length-weight relationships (Brodziak 2012) for all pollock measured from stomach contents of our focal predators. We used multinomial logistic regression to estimates year-specific age compositions of pollock consumed. Small sample sizes precluded spatially-explicit age compositions.

# Bottom trawl survey data were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for details about survey design and data collection. 

# Food habits data were obtained from the Resource Ecology and Ecosystem Modeling (REEM) program of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See Livingston et al. (2017) for methods. 

# References:
# Brodziak J 2012. Fitting length-weight relationships with linear regression using the log-transformed allometric model with bias-correction. NOAA Technical Memorandum PIFSC-H-12-03. 
# Livingston PA, K Aydin, TW Buckley, GM Lang, M-S Yang, and BS Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fish. 100(4):443–470.
# von Bertalanffy, L. 1938. A quantitative theory of organic growth. Human Biology 10:181–213.
# von Szalay PG and NW Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. National Oceanic and Atmospheric Administration. Techn Mem NMFS-AFSC-325. 

setwd("~/Documents/AFSC/GOApollock_ESP/")
require(dplyr)
require(tidyr)
require(VGAM) 
require(ggplot2)

#########################################################
### QUANTIFY LENGTH-AGE RELATIONSHIPS, Walleye Pollock

# Input pollock length, weight, and age data from the AFSC bottom trawl survey (random samples):
specimen = read.csv(unz("1_Data/race_specimen.csv.zip", "race_specimen.csv"), header=T, skip=7, na.strings=c("","NA"))
  specimen = subset(specimen, Common.Name =- "walleye pollock") # select just WEP
# Exclude data from 1984 and 1987 (survey methods were standardized in 1990) and remove NA ages:
specimen = subset(specimen, Year >= 1990 & is.finite(Age..years.)) 
  specimen = subset(specimen, Age..years. = 0 & Length..mm. > 300) # remove obvious outliers (i.e., age 0 fish > 300 mm):

# Create unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to lengths data below):
specimen$Haul.Number_Join = with(specimen, paste(Vessel.Number, Cruise.Number, Haul.Number, sep="_"))  
  
# Exclude SEAK (GOA main assessment area only):
specimen.assess = subset(specimen, Starting.Longitude..dd. < -140)

# Set starting values for VBGF (based on raw data):
theta = c(650, 0.2, 0) 

# Estimate von Bertalannfy growth parameters:
SSQ = function(theta, x) {
  Linf = theta[1]
  K = theta[2]
  t0 = theta[3]
  epsilon = rep(0, length(specimen.assess$Age..years.))
  lpred = rep(0, length(specimen.assess$Age..years.))
  for (i in 1:length(specimen.assess$Age..years.)) {
    lpred[i] = Linf * (1 - exp(-K * (specimen.assess$Age..years.[i] - t0)))
    epsilon[i] = (specimen.assess$Length..mm.[i] - lpred[i])^2
  }
  ssq = sum(epsilon)
  return(ssq) }

# Solve using least squares:
out = optim(theta, fn = SSQ, method = "BFGS", x = specimen.assess$Age..years., hessian = TRUE)
  out$V = solve(out$hessian)  # solve the hessian
  out$SE = sqrt(diag(out$V))  # estimate SE
  out$R = out$V/(out$SE %o% out$SE)  # correlation
  out$par; out$SE

# Estimate pollock age from length:
Linf = out$par[1]
K = out$par[2]
t0 = out$par[3]
specimen.assess$pred.age = ((-(log(-(specimen.assess$Length..mm./Linf)+1)))/K)+t0

# Group pollock by age classes according to the Gulf of Alaska stock assessment:
specimen.assess$age.class = with(specimen.assess, 
    ifelse(pred.age < 1, "0", 
    ifelse(pred.age < 2, "1", 
    ifelse(pred.age < 3, "2", "3+")))) 

#########################################################
### QUANTIFY LENGTH-WEIGHT RELATIONSHIPS, Walleye Pollock

# Allometric model with bias-correction (Brodziak 2012):
specimen.assess$logW = log(specimen.assess$Weight..gm.)
specimen.assess$logL = log(specimen.assess$Length..mm.)
  specimen.assess = subset(specimen.assess, logW !="NA")

L_W = lm(logW ~ logL, data=specimen.assess)
summary(L_W)
  a = L_W$coefficients[1]
  b = L_W$coefficients[2]

# Calculate correction factor for predicting weight on original scale:
syx = summary(L_W)$sigma
cf = exp((syx^2)/2) 
specimen.assess$pred.logW = predict(L_W, 
              data.frame(logL=log(specimen.assess$Length..mm.)), interval="c") 
specimen.assess$pred.unbiW = cf *(exp(specimen.assess$pred.logW))

# Estimate pollock weight from length:
specimen.assess$logL = log(specimen.assess$Length..mm.)
specimen.assess$pred.logW = predict(L_W, 
                data.frame(logL=log(specimen.assess$Length..mm.)), interval="c")
specimen.assess$pred.unbiW = cf *(exp(specimen.assess$pred.logW)) 

#########################################################
### PREDICT WEIGHT FROM LENGTH OF POLLOCK PREY 

# Input food habits data for Arrowtooth Flounder:
prey.lengths = read.csv("1_Data/GOA_RawPL_AllDiets.csv")

# Create unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to specimen data above):
prey.lengths$Haul.Number_Join = with(prey.lengths, paste(VESSEL, CRUISE, HAUL, sep="_"))  

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
prey.lengths = prey.lengths %>%
  rename(Year = YEAR)
prey.lengths = subset(prey.lengths, Year >= 1990)

# Select pollock as prey:
prey.lengths = subset(prey.lengths, PREY_NODC == 8791030701)

# Link to specimen data (at add location information):
prey.lengths_loc = prey.lengths %>% left_join(specimen[,c("Haul.Number_Join", "Starting.Longitude..dd.", "Starting.Latitude..dd.")])

# Exclude SEAK (GOA main assessment area only):
prey.lengths_assess = subset(prey.lengths_loc, Starting.Longitude..dd. < -140)

# Predict pollock age from length:
prey.lengths_assess$pred.age = (log(-(prey.lengths_assess$PREY_SZ1/Linf)+1))/-K

# Group pollock by age classes according to the Gulf of Alaska stock assessment:
prey.lengths_assess$age.class = with(prey.lengths_assess, 
  ifelse(pred.age < 1, "0", 
     ifelse(pred.age < 2, "1", 
         ifelse(pred.age < 3, "2", "3+"))))
prey.lengths_assess$age.class = as.factor(prey.lengths_assess$age.class)
prey.lengths_assess$age.class[is.na(prey.lengths_assess$age.class)] = "3+"

# Predict pollock weight from length:
prey.lengths_assess$pred.logW = predict(L_W, 
              data.frame(logL=log(prey.lengths_assess$PREY_SZ1)), interval="c") 
prey.lengths_assess$pred.logW = as.vector(prey.lengths_assess$pred.logW[,1])
prey.lengths_assess$pred.unbiW = cf *(exp(prey.lengths_assess$pred.logW))

# Remove predators not included in assessment-based estimates of total predator biomass:
WEP.prey_ATF = subset(prey.lengths_assess, 
                      PRED_NODC == 8857040102 & PRED_LEN >= 19)
WEP.prey_PC = subset(prey.lengths_assess, 
                      PRED_NODC == 8791030401 & PRED_LEN >= 0)
WEP.prey_PH = subset(prey.lengths_assess, 
                      PRED_NODC == 8857041901 & PRED_LEN >= 82)
WEP.prey_SBL = subset(prey.lengths_assess, 
                      PRED_NODC == 8827020101 & PRED_LEN >= 45)
WEP.prey_WEP = subset(prey.lengths_assess, 
                      PRED_NODC == 8791030701 & PRED_LEN >= 37)

#########################################################
### CALCULATE PROPORTIONS OF POLLOCK CONSUMED, BY AGE CLASS 

# Calculate raw proportions of pollock consumed by age class:
### ATF:
WEP.prop_ATF = WEP.prey_ATF %>%
  group_by(Year, age.class) %>%
  summarise(sum.WT = sum(pred.unbiW))

### PC:
WEP.prop_PC = WEP.prey_PC %>%
  group_by(Year, age.class) %>%
  summarise(sum.WT = sum(pred.unbiW))

### PH:
WEP.prop_PH = WEP.prey_PH %>%
  group_by(Year, age.class) %>%
  summarise(sum.WT = sum(pred.unbiW))

### SBL:
WEP.prop_SBL = WEP.prey_SBL %>%
  group_by(Year, age.class) %>%
  summarise(sum.WT = sum(pred.unbiW))

### WEP:
WEP.prop_WEP = WEP.prey_WEP %>%
  group_by(PRED_NODC, Year, age.class) %>%
  summarise(sum.WT = sum(pred.unbiW))

# Estimate age compositions of pollock consumed, by survey year:
### ATF:
WEP.age_ATF = WEP.prey_ATF %>%
  pivot_wider(id_cols = Year, 
              names_from = age.class, 
              values_from = pred.unbiW, 
              values_fn = sum, 
              values_fill = 0) %>%
  as.data.frame()

ATF.Year = WEP.age_ATF$Year 
ATF.matrix = as.matrix(round(WEP.age_ATF[ ,2:5], digits=0)) 
  rownames(ATF.matrix) = ATF.Year 

ATF.fit_gam = vgam(ATF.matrix ~ s(Year), multinomial, data=WEP.age_ATF)
  summary(ATF.fit_gam)

WEP.age_ATF.gam = as.data.frame(predict(ATF.fit_gam, type="response"))
WEP.age_ATF.gam$Year = ATF.Year

WEP.ages_ATF.gam = WEP.age_ATF.gam %>% 
  pivot_longer(cols = 1:4,
               names_to = "age.class", 
               values_to = "propWT.age") %>%
  mutate(Predator = "ATF") %>%
  as.data.frame()

WEP.ages_ATF.gam$age.class = ordered(WEP.ages_ATF.gam$age.class, levels = c("3+", "2", "1", "0"))
  saveRDS(WEP.ages_ATF.gam, 
    file = "5_AgeCompositionsPrey/propWEP_ages_ATF.rds")
# WEP.ages_ATF.gam = readRDS("5_AgeCompositionsPrey/propWEP_ages_ATF.rds")

  
### PC:
WEP.age_PC = WEP.prey_PC %>%
  pivot_wider(id_cols = Year, 
              names_from = age.class, 
              values_from = pred.unbiW, 
              values_fn = sum, 
              values_fill = 0) %>%
  as.data.frame()

PC.Year = WEP.age_PC$Year 
PC.matrix = as.matrix(round(WEP.age_PC[ ,2:5], digits=0)) 
  rownames(PC.matrix) = PC.Year 

PC.fit_gam = vgam(PC.matrix ~ s(Year), multinomial, data=WEP.age_PC)
  summary(PC.fit_gam)

WEP.age_PC.gam = as.data.frame(predict(PC.fit_gam, type="response"))
WEP.age_PC.gam$Year = PC.Year

WEP.ages_PC.gam = WEP.age_PC.gam %>% 
  pivot_longer(cols = 1:4,
               names_to = "age.class", 
               values_to = "propWT.age") %>%
  mutate(Predator = "PC") %>%
  as.data.frame()

WEP.ages_PC.gam$age.class = ordered(WEP.ages_PC.gam$age.class, levels = c("3+", "2", "1", "0"))
  saveRDS(WEP.ages_PC.gam, 
    file = "5_AgeCompositionsPrey/propWEP_ages_PC.rds")
# WEP.ages_PC.gam = readRDS("5_AgeCompositionsPrey/propWEP_ages_PC.rds")


### PH:
WEP.age_PH = WEP.prey_PH %>%
  pivot_wider(id_cols = Year, 
              names_from = age.class, 
              values_from = pred.unbiW, 
              values_fn = sum, 
              values_fill = 0) %>%
  as.data.frame()

PH.Year = WEP.age_PH$Year 
PH.matrix = as.matrix(round(WEP.age_PH[ ,2:5], digits=0)) 
  rownames(PH.matrix) = PH.Year 

PH.fit_gam = vgam(PH.matrix ~ s(Year), multinomial, data=WEP.age_PH)
  summary(PH.fit_gam)

WEP.age_PH.gam = as.data.frame(predict(PH.fit_gam, type="response"))
WEP.age_PH.gam$Year = PH.Year

# Assign overall mean age compositions when data were missing:
mean.WEP.age3_PH = mean(WEP.age_PH.gam$`3+`)
mean.WEP.age2_PH = mean(WEP.age_PH.gam$`2`)
mean.WEP.age1_PH = mean(WEP.age_PH.gam$`1`)
mean.WEP.age0_PH = mean(WEP.age_PH.gam$`0`)

WEP.age_PH.gam = add_row(WEP.age_PH.gam)
  WEP.age_PH.gam[nrow(WEP.age_PH.gam),] = c(mean.WEP.age3_PH, 
                                            mean.WEP.age1_PH, 
                                            mean.WEP.age0_PH,                                                     mean.WEP.age2_PH, 2005)

WEP.ages_PH.gam = WEP.age_PH.gam %>% 
  pivot_longer(cols = 1:4,
               names_to = "age.class", 
               values_to = "propWT.age") %>%
  mutate(Predator = "PH") %>%
  as.data.frame()

WEP.ages_PH.gam$age.class = ordered(WEP.ages_PH.gam$age.class, levels = c("3+", "2", "1", "0"))
  saveRDS(WEP.ages_PH.gam, 
    file = "5_AgeCompositionsPrey/propWEP_ages_PH.rds")
# WEP.ages_PH.gam = readRDS("5_AgeCompositionsPrey/propWEP_ages_PH.rds") 
  

### SBL:
# Remove survey years with fewer than three (two here) measurable pollock:
WEP.prey_SBL %>% 
  group_by(Year) %>% 
  summarise(length(PREY_SZ1))
  
WEP.age_SBL = WEP.prey_SBL %>%
  pivot_wider(id_cols = Year, 
              names_from = age.class, 
              values_from = pred.unbiW, 
              values_fn = sum, 
              values_fill = 0) %>%
  as.data.frame()

SBL.Year = WEP.age_SBL$Year 
SBL.matrix = as.matrix(round(WEP.age_SBL[ ,2:5], digits=0)) 
  rownames(SBL.matrix) = SBL.Year 

SBL.fit_gam = vgam(SBL.matrix ~ s(Year), multinomial, data=WEP.age_SBL)
  summary(SBL.fit_gam)

WEP.age_SBL.gam = as.data.frame(predict(SBL.fit_gam, type="response"))
WEP.age_SBL.gam$Year = SBL.Year

# Assign overall mean age compositions when data were missing:
mean.WEP.age3_SBL = mean(WEP.age_SBL.gam$`3+`)
mean.WEP.age2_SBL = mean(WEP.age_SBL.gam$`2`)
mean.WEP.age1_SBL = mean(WEP.age_SBL.gam$`1`)
mean.WEP.age0_SBL = mean(WEP.age_SBL.gam$`0`)

WEP.age_SBL.gam = add_row(WEP.age_SBL.gam)
  WEP.age_SBL.gam[nrow(WEP.age_SBL.gam),] = c(mean.WEP.age3_SBL, 
                                              mean.WEP.age2_SBL, 
                                              mean.WEP.age0_SBL, 
                                              mean.WEP.age1_SBL, 1996)   
WEP.age_SBL.gam = add_row(WEP.age_SBL.gam)
  WEP.age_SBL.gam[nrow(WEP.age_SBL.gam),] = c(mean.WEP.age3_SBL, 
                                              mean.WEP.age2_SBL, 
                                              mean.WEP.age0_SBL, 
                                              mean.WEP.age1_SBL, 2005)   
WEP.age_SBL.gam = add_row(WEP.age_SBL.gam)
  WEP.age_SBL.gam[nrow(WEP.age_SBL.gam),] = c(mean.WEP.age3_SBL, 
                                              mean.WEP.age2_SBL, 
                                              mean.WEP.age0_SBL, 
                                              mean.WEP.age1_SBL, 2013)
WEP.age_SBL.gam = add_row(WEP.age_SBL.gam)
  WEP.age_SBL.gam[nrow(WEP.age_SBL.gam),] = c(mean.WEP.age3_SBL, 
                                              mean.WEP.age2_SBL, 
                                              mean.WEP.age0_SBL, 
                                              mean.WEP.age1_SBL, 2015) 
WEP.age_SBL.gam = add_row(WEP.age_SBL.gam)
  WEP.age_SBL.gam[nrow(WEP.age_SBL.gam),] = c(mean.WEP.age3_SBL, 
                                              mean.WEP.age2_SBL, 
                                              mean.WEP.age0_SBL, 
                                              mean.WEP.age1_SBL, 2017) 
WEP.age_SBL.gam = add_row(WEP.age_SBL.gam)
  WEP.age_SBL.gam[nrow(WEP.age_SBL.gam),] = c(mean.WEP.age3_SBL, 
                                              mean.WEP.age2_SBL, 
                                              mean.WEP.age0_SBL, 
                                              mean.WEP.age1_SBL, 2019) 

WEP.ages_SBL.gam = WEP.age_SBL.gam %>% 
  pivot_longer(cols = 1:4,
               names_to = "age.class", 
               values_to = "propWT.age") %>%
  mutate(Predator = "SBL") %>%
  as.data.frame()

WEP.ages_SBL.gam$age.class = ordered(WEP.ages_SBL.gam$age.class, levels = c("3+", "2", "1", "0"))

  saveRDS(WEP.ages_SBL.gam, 
    file = "5_AgeCompositionsPrey/propWEP_ages_SBL.rds")
# WEP.ages_SBL.gam = readRDS("5_AgeCompositionsPrey/propWEP_ages_SBL.rds")


### WEP:
# Remove survey years with fewer than three measurable pollock:
WEP.prey_WEP %>% 
  group_by(Year) %>% 
  summarise(length(PREY_SZ1))
 
WEP.age_WEP = WEP.prey_WEP %>%
  pivot_wider(id_cols = Year, 
              names_from = age.class, 
              values_from = pred.unbiW, 
              values_fn = sum, 
              values_fill = 1) %>%
  as.data.frame()
WEP.age_WEP$`3+` = 0

WEP.Year = WEP.age_WEP$Year 
WEP.matrix = as.matrix(round(WEP.age_WEP[ ,2:4], digits=0)) 
  rownames(WEP.matrix) = WEP.Year 

  
WEP.fit_gam = vgam(WEP.matrix ~ s(Year), multinomial, data=WEP.age_WEP)
  summary(WEP.fit_gam)

WEP.age_WEP.gam = as.data.frame(predict(WEP.fit_gam, type="response"))
WEP.age_WEP.gam$Year = WEP.Year

# Assign overall mean age compositions when data were missing:
mean.WEP.age2_WEP = mean(WEP.age_WEP.gam$`2`)
mean.WEP.age1_WEP = mean(WEP.age_WEP.gam$`1`)
mean.WEP.age0_WEP = mean(WEP.age_WEP.gam$`0`)

WEP.age_WEP.gam = add_row(WEP.age_WEP.gam)
  WEP.age_WEP.gam[nrow(WEP.age_WEP.gam),] = c(mean.WEP.age0_WEP, 
                                              mean.WEP.age1_WEP, 
                                              mean.WEP.age2_WEP, 1999)   
WEP.age_WEP.gam = add_row(WEP.age_WEP.gam)
  WEP.age_WEP.gam[nrow(WEP.age_WEP.gam),] = c(mean.WEP.age0_WEP, 
                                              mean.WEP.age1_WEP, 
                                              mean.WEP.age2_WEP, 2005)   
WEP.age_WEP.gam = add_row(WEP.age_WEP.gam)
  WEP.age_WEP.gam[nrow(WEP.age_WEP.gam),] = c(mean.WEP.age0_WEP, 
                                              mean.WEP.age1_WEP, 
                                              mean.WEP.age2_WEP, 2011)

WEP.age_WEP.gam$`3+` = 0 # no age 3+ WEP consumed

WEP.ages_WEP.gam = WEP.age_WEP.gam %>% 
  pivot_longer(cols = c(1:3,5),
               names_to = "age.class", 
               values_to = "propWT.age") %>%
  mutate(Predator = "WEP") %>%
  as.data.frame()

WEP.ages_WEP.gam$age.class = ordered(WEP.ages_WEP.gam$age.class, levels = c("3+", "2", "1", "0"))

  saveRDS(WEP.ages_WEP.gam, 
    file = "5_AgeCompositionsPrey/propWEP_ages_WEP.rds")
# WEP.ages_WEP.gam = readRDS("5_AgeCompositionsPrey/propWEP_ages_WEP.rds")