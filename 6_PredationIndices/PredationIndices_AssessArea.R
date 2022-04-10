# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)
  # Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# Citation: 
# Dorn MW and CL Barnes. In revision. Time-varying predation as a modifier of constant natural mortality for Gulf of Alaska walleye pollock. Fish Res. VSI Natural Mortality. Forthcoming.

# This script includes code necessary to calculate an index of predation for Walleye Pollock within the main stock assessment assessment area of the Gulf of Alaska (MT per year; 1990 to 2019). 

# Total biomass estimates were obtained from the most recent stock assessment for Arrowtooth Flounder (Spies et al. 2019), Pacific Cod (Barbeaux et al. 2019), Pacific Halibut (Stewart and Hicks 2020), Sablefish (Hanselmann et al. 2019), and Walleye Pollock (Dorn et al. 2019). Bottom trawl survey data were collected by the Resource Assessment and Conservation Engineering (RACE) Division and are publicly accessible at https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm (for methods, see von Szalay et al. 2016) and used to predict relative predator density. Food habits data were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php (for methods, see Livingston et al. 2017). 

# References:
# Barbeaux S, K Aydin, B Fissel, K Holsman, B Laurel, W Palsson, L Rogers, K Shotwell, Q Yang, and S Zador. 2019. Assessment of the Pacific cod stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report. 140pp.
# Barnes CL, Beaudreau AH, Dorn MW, Holsman KK, and Mueter FJ. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecol Appl. 30(7):e02141.
# Clark WG and SR Hare. 2006. Assessment and management of Pacific halibut: data, methods, and policy. IPHC Scientific Report 83. 
# Dorn M, AL Deary, BE Fissell, DT Jones, NE Lauffenburger, W.A. Palsson, LA Rogers, SK Shotwell, KA Spalinger, and S Zador. 2019. Assessment of the Walleye Pollock stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report. 161pp.
# Hanselman DH, CJ Rodgveller, KH Fenske, SK Shotwell, KB Echave, PW Malecha, and CR Lunsford. 2019. Assessment of the Sablefish stock in Alaska. North Pacific Fishery Management Council Bering Sea, Aleutian Islands, and Gulf of Alaska SAFE Report. 263pp.
# Livingston PA, K Aydin, TW Buckley, GM Lang, M-S Yang, and BS Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes 100(4):443–470.
# Spies I, K Aydin, JN Ianelli, and W Palsson. 2019. Assessment of the Arrowtooth Flounder stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report. 92pp.
# Stewart I and A Hicks. 2020. Assessment of the Pacific halibut (Hippoglossus stenolepis) stock at the end of 2019. International Pacific Halibut Commission Report. 32pp.
# von Szalay PG and NW Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. National Oceanic and Atmospheric Administration. Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/AFSC/GOApollock_ESP/")
require(dplyr)
require(tidyr)
require(readxl)
require(ggplot2)

##################################################################
# Input total biomass estimates from most recent stock assessment (ATF = age 1+, Gulf-wide, MT) Model 19.0:
ATF.M = data.frame("Year" = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019), "TotBio_t" = c(
  1580890,
  1700620,
  1701290,
  1753150,
  1865510,
  1937550,
  1964790,
  1940850,
  1844400,
  1708910,
  1586650,
  1463320,
  1378080,
  1333540))


# PC = age 0+, Gulf-wide (MT)
  # Barbeaux et al. 2019 (Model 19.14.48c)
PC.M = data.frame("Year" = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019), "TotBio_t" = c(
  759213,
  654597,
  518686,
  375787,
  318708,
  325922,
  280816,
  268873,
  342596,
  425866,
  445224,
  381875,
  155394,
  141458))


# PH = age 8+, GOA  as a fraction of coast-wide biomass (mill lb -> MT)
  # Stewart and Hicks 2019 [CW long model] 
# Coast-wide estimates were adjusted by proportions of > 32" fish sampled in the GOA (IPHC Regulatory Areas 4A, 3B, 3A, and 2C) as part of IPHC's setline survey. No setline survey data were available prior to 1998; overall mean proportions were assigned to those years. Assessment-based estimates were unavailable for 1990. Thus, biomass in 1990 was calculated (from 1996) based on trends from the bottom trawl survey, which were highly correlated with setline survey biomass. 

PH.M = data.frame("Year" = c(1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019), "TotBio_mill.lb" = c(
  1510, 
  1978,
  1786, 
  1423, 
  1296, 
  1057, 
  995,
  853,
  766,  
  807, 
  704,  
  611,
  579)) 

# IPHC setline data (to calculate proportions of biomass in GOA):
setline_98_09 = read_excel("1_Data/SetPacificHalibutData_1998_2009.xlsx")
setline_10_20 = read_excel("1_Data/SetPacificHalibutData_2010_2020.xlsx")
setline = rbind(setline_98_09, setline_10_20)
setline = setline %>%
  rename(StatArea = `IPHC Stat Area`)
  
# Group statistical areas (i.e., as defined in 2003) by IPHC regulatory area:
setline$StatArea = as.numeric(setline$StatArea)
setline$RegArea = with(setline, 
  ifelse(StatArea >= 006 & StatArea <= 050, "2A",
  ifelse(StatArea >= 060 & StatArea <= 135, "2B",
  ifelse(StatArea >= 140 & StatArea <= 184, "2C", 
  ifelse(StatArea >= 185 & StatArea <= 281, "3A",
  ifelse(StatArea >= 290 & StatArea <= 340, "3B",
  ifelse(StatArea >= 350 & StatArea <= 395, "4A",
  ifelse(StatArea >= 400 & StatArea <= 510, "4B", "Other"))))))))

setline$Region = with(setline,
     ifelse(RegArea == "2C" | 
            RegArea == "3A" | 
            RegArea == "3B" | 
            RegArea == "4A", "GOA", "Not_GOA"))
setline = as.data.frame(setline)
  
GOA.prop = setline %>%
  group_by(Year, Region) %>%
  summarise(sum.Reg = `O32 Pacific halibut weight`) %>%
  pivot_wider(names_from = Region, values_from = sum.Reg, values_fn = sum) %>%
  mutate(prop.GOA = GOA / (GOA + Not_GOA))
GOA.prop$Year = as.numeric(GOA.prop$Year)

PH.M = PH.M %>% left_join(GOA.prop[,c("Year", "prop.GOA")])

# Assign overall mean proportions before the start of IPHC survey (GOA-wide):
PH.M$prop.GOA = with(PH.M, 
                     ifelse(is.na(prop.GOA), mean(prop.GOA, na.rm = T), 
                            prop.GOA))

# Convert mill lb to metric tons:
PH.M$TotBio_t = PH.M$TotBio_mill.lb *  PH.M$prop.GOA * 453.59237

# Estimate 1990 biomass based on trend in bottom trawl survey:
load("2_RelativePredatorDensities/trawl_survey_data_sub.rda")
  trawl.sum = trawl_wide %>%
    group_by(Year) %>%
    summarise(PH.bio = sum(PH))
PH.M = add_row(PH.M)
  PH.M[nrow(PH.M),] = c(1990, NA, NA, (PH.M[3,4] * (trawl.sum[1,2]/trawl.sum[4,2])))   

# SBL = age 2+, Gulf-wide as sum wGOA, cGOA, wYak, and eYak/SE (kilo MT)
  # Hanselman et al. 2019
SBL.M = data.frame("Year" = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019), "TotBio_t" = c(
  ((42+121+46+60)*1000),
  ((30+110+56+84)*1000),
  ((29+97+35+55)*1000),
  ((30+85+27+51)*1000),
  ((42+83+22+46)*1000),
  ((42+101+26+43)*1000),
  ((38+95+26+47)*1000),
  ((29+85+29+48)*1000),
  ((29+78+22+40)*1000),
  ((24+84+31+44)*1000),
  ((22+71+19+43)*1000),
  ((21+55+21+36)*1000),
  ((36+88+36+46)*1000),
  ((81+146+46+83)*1000)))


# WEP = age 3+, GOA (-140 westward; MT);
  # orn et al. 2019
WEP.M = data.frame("Year" = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019), "TotBio_t" = c(
  (1452 * 1000),
  (1733 * 1000),
  (1012 * 1000),
  (740 * 1000),
  (620 * 1000),
  (994 * 1000),
  (688 * 1000),
  (537 * 1000),
  (1046 * 1000),
  (1164 * 1000),
  (1087 * 1000),
  (2186 * 1000),
  (1563 * 1000),
  (941 * 1000)))

# Convert MT to g (comparable with units for annual ration and proportions of pollock consumed):
ATF.M$TotBio_g = ATF.M$TotBio_t * 1000000
PC.M$TotBio_g = PC.M$TotBio_t * 1000000
PH.M$TotBio_g = PH.M$TotBio_t * 1000000
SBL.M$TotBio_g = SBL.M$TotBio_t * 1000000
WEP.M$TotBio_g = WEP.M$TotBio_t * 1000000

#########################################################
# Input relative predator densities (unitless) predicted from GAM(M)s by survey year:
load("2_RelativePredatorDensities/normATFabun.rda")
  ATF.density = normAbun_ATF
ATF.density$Year = as.numeric(as.character(ATF.density$Year))

load("2_RelativePredatorDensities/normPCabun.rda")
  PC.density = normAbun_PC
PC.density$Year = as.numeric(as.character(PC.density$Year))

load("2_RelativePredatorDensities/normPHabun.rda")
  PH.density = normAbun_PH
PH.density$Year = as.numeric(as.character(PH.density$Year))

load("2_RelativePredatorDensities/normSBLabun.rda")
  SBL.density = normAbun_SBL
SBL.density$Year = as.numeric(as.character(SBL.density$Year))

load("2_RelativePredatorDensities/normWEPabun.rda")
  WEP.density = normAbun_WEP
WEP.density$Year = as.numeric(as.character(WEP.density$Year))

# Exclude locations outside the main (i.e., 'GOA') Walleye Pollock stock assessment (areas east of -140):
ATF.goa = subset(ATF.density, Starting.Longitude..dd. < -140)
PC.goa = subset(PC.density, Starting.Longitude..dd. < -140)
PH.goa = subset(PH.density, Starting.Longitude..dd. < -140)
SBL.goa = subset(SBL.density, Starting.Longitude..dd. < -140)
WEP.goa = subset(WEP.density, Starting.Longitude..dd. < -140)
  
# Sum year-specific relative predator densities: 
ATF.goa_yr = ATF.goa %>%
  group_by(Year) %>%
  summarise(ATF.Abun = sum(ATF.normAbun))

PC.goa_yr = PC.goa %>%
  group_by(Year) %>%
  summarise(PC.Abun = sum(PC.normAbun))

PH.goa_yr = PH.goa %>%
  group_by(Year) %>%
  summarise(PH.Abun = sum(PH.normAbun)) %>%
  as.data.frame()

SBL.goa_yr = SBL.goa %>%
  group_by(Year) %>%
  summarise(SBL.Abun = sum(SBL.normAbun))

WEP.goa_yr = WEP.goa %>%
  group_by(Year) %>%
  summarise(WEP.Abun = sum(WEP.normAbun))

# Join to predator biomass data:
ATF.M = ATF.M %>% left_join(ATF.goa_yr)
PC.M = PC.M %>% left_join(PC.goa_yr)
PH.M = PH.M %>% left_join(PH.goa_yr)
SBL.M = SBL.M %>% left_join(SBL.goa_yr)
WEP.M = WEP.M %>% left_join(WEP.goa_yr)

#########################################################
# Input mean annual rations (g/g/y) calculated using Wisconsin bioenergetics models; already restricted to GOA west of -140:
ATF.ration = readRDS("3_MeanAnnualRations/Rations_ATF.rds")
PC.ration = readRDS("3_MeanAnnualRations/Rations_PC.rds")
PH.ration = readRDS("3_MeanAnnualRations/Rations_PH.rds")
SBL.ration = readRDS("3_MeanAnnualRations/Rations_SBL.rds")
WEP.ration = readRDS("3_MeanAnnualRations/Rations_WEP.rds")

# Calculate theorized maximum consumption (Cmax) with relative foraging rate (RFR) = 1:
ATF.Cmax = ATF.ration %>%
  group_by(Year) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy, na.rm = T))

PC.Cmax = PC.ration %>%
  group_by(Year) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy, na.rm = T))

PH.Cmax = PH.ration %>%
  group_by(Year) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy, na.rm = T))

SBL.Cmax = SBL.ration %>%
  group_by(Year) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy, na.rm = T))

# Add mean estimates for years without data:
mean.Cmax_SBL = mean(SBL.Cmax$Annual_Cmax)
SBL.Cmax = add_row(SBL.Cmax)
  SBL.Cmax[nrow(SBL.Cmax),] = list(2013, mean.Cmax_SBL)
SBL.Cmax = add_row(SBL.Cmax)
  SBL.Cmax[nrow(SBL.Cmax),] = list(2015, mean.Cmax_SBL) 
SBL.Cmax = add_row(SBL.Cmax)
  SBL.Cmax[nrow(SBL.Cmax),] = list(2017, mean.Cmax_SBL)  
SBL.Cmax = add_row(SBL.Cmax)
  SBL.Cmax[nrow(SBL.Cmax),] = list(2019, mean.Cmax_SBL) 

WEP.Cmax = WEP.ration %>%
  group_by(Year) %>%
  summarise(Annual_Cmax = mean(Cmax_ggy, na.rm = T))

# Join to total predator biomass and relative predator density estimates:
ATF.M = ATF.M %>% left_join(ATF.Cmax)
PC.M = PC.M %>% left_join(PC.Cmax)
PH.M = PH.M %>% left_join(PH.Cmax)
SBL.M = SBL.M %>% left_join(SBL.Cmax)
WEP.M = WEP.M %>% left_join(WEP.Cmax)

#########################################################
# Input year-specific proportions of pollock consumed (length- and biomass-weighted; unitless); already restricted to GOA west of -140:
load("4_ProportionsPollockConsumed/propWEP_ATF_19.rda")
ATF.prey = subset(ATF.prey_weighted.assess, select=-c(geometry))

load("4_ProportionsPollockConsumed/propWEP_PC_0.rda")
PC.prey = subset(PC.prey_weighted.assess, select=-c(geometry))

load("4_ProportionsPollockConsumed/propWEP_PH_82.rda")
PH.prey = subset(PH.prey_weighted.assess, select=-c(geometry))

load("4_ProportionsPollockConsumed/propWEP_SBL_45.rda")
SBL.prey = subset(SBL.prey_weighted.assess, select=-c(geometry))

load("4_ProportionsPollockConsumed/propWEP_WEP_37.rda")
WEP.prey = subset(WEP.prey_weighted.assess, select=-c(geometry))

# Calculate year-specific proportions of prey consumed:
ATF.prey = ATF.prey %>%
  group_by(Year, poll) %>% 
  summarise(sumWT = sum(prey.WT, na.rm = T)) %>%
  mutate(propWT = sumWT / sum(sumWT))

PC.prey = PC.prey %>%
  group_by(Year, poll) %>% 
  summarise(sumWT = sum(prey.WT, na.rm = T)) %>%
  mutate(propWT = sumWT / sum(sumWT))

PH.prey = PH.prey %>%
  group_by(Year, poll) %>% 
  summarise(sumWT = sum(prey.WT, na.rm = T)) %>%
  mutate(propWT = sumWT / sum(sumWT))

SBL.prey = SBL.prey %>%
  group_by(Year, poll) %>% 
  summarise(sumWT = sum(prey.WT, na.rm = T)) %>%
  mutate(propWT = sumWT / sum(sumWT)) %>%
  as.data.frame()

# Assign overall mean proportions within assessment area (west of -140) for years that Sablefish were not sampled for food habits:
mean.WEP_SBL = mean(SBL.prey$propWT)
SBL.prey = add_row(SBL.prey)
  SBL.prey[nrow(SBL.prey),] = c(2013, "pollock", NA, mean.WEP_SBL)
SBL.prey = add_row(SBL.prey)
  SBL.prey[nrow(SBL.prey),] = c(2015, "pollock", NA, mean.WEP_SBL)
SBL.prey = add_row(SBL.prey)
  SBL.prey[nrow(SBL.prey),] = c(2017, "pollock", NA, mean.WEP_SBL)
SBL.prey = add_row(SBL.prey)
  SBL.prey[nrow(SBL.prey),] = c(2019, "pollock", NA, mean.WEP_SBL)
SBL.prey = SBL.prey %>%
  mutate_at(c("Year", "sumWT", "propWT"), as.numeric)
  
WEP.prey = WEP.prey %>%
  group_by(Year, poll) %>% 
  summarise(sumWT = sum(prey.WT, na.rm = T)) %>%
  mutate(propWT = sumWT / sum(sumWT)) %>%
  as.data.frame()

# Add zeros for years when pollock were not consumed:
WEP.prey = add_row(WEP.prey)
  WEP.prey[nrow(WEP.prey),] = c(2005, "pollock", 0, 0)
WEP.prey = add_row(WEP.prey)
  WEP.prey[nrow(WEP.prey),] = c(2011, "pollock", 0, 0)
WEP.prey = WEP.prey %>%
  mutate_at(c("Year", "sumWT", "propWT"), as.numeric)

# Select only proportions of Walleye Pollock:
ATF.prey.sub = subset(ATF.prey, poll == "pollock")
PC.prey.sub = subset(PC.prey, poll == "pollock")
PH.prey.sub = subset(PH.prey, poll == "pollock")
SBL.prey.sub = subset(SBL.prey, poll == "pollock")
WEP.prey.sub = subset(WEP.prey, poll == "pollock")

# Join to total predator biomass, relative predator density, and mean annual ration estimates:
ATF.M = ATF.M %>% left_join(ATF.prey.sub)
PC.M = PC.M %>% left_join(PC.prey.sub)
PH.M = PH.M %>% left_join(PH.prey.sub)
SBL.M = SBL.M %>% left_join(SBL.prey.sub)
WEP.M = WEP.M %>% left_join(WEP.prey.sub)

#########################################################
# Calculate year-specific predation by Arrowtooth Flounder (g; all pollock age classes combined) and convert g of pollock consumed to metric tons::
ATF.M$WEP.cons_g = with(ATF.M, 
                       TotBio_g * ATF.Abun * Annual_Cmax * propWT)
PC.M$WEP.cons_g = with(PC.M, 
                       TotBio_g * PC.Abun * Annual_Cmax * propWT)
PH.M$WEP.cons_g = with(PH.M, 
                       TotBio_g * PH.Abun * Annual_Cmax * propWT)
SBL.M$WEP.cons_g = with(SBL.M, 
                       TotBio_g * SBL.Abun * Annual_Cmax * propWT)
WEP.M$WEP.cons_g = with(WEP.M, 
                       TotBio_g * WEP.Abun * Annual_Cmax * propWT)

#########################################################
# Calculate age-specific predation on pollock by Arrowtooth Flounder (g):
WEP_age_ATF = readRDS("5_AgeCompositionsPrey/propWEP_ages_ATF.rds")
WEP_age_PC = readRDS("5_AgeCompositionsPrey/propWEP_ages_PC.rds")
WEP_age_PH = readRDS("5_AgeCompositionsPrey/propWEP_ages_PH.rds")
WEP_age_SBL = readRDS("5_AgeCompositionsPrey/propWEP_ages_SBL.rds")
WEP_age_WEP = readRDS("5_AgeCompositionsPrey/propWEP_ages_WEP.rds")

# Join to total predator biomass estimates, relative predator densities, mean annual rations, and proportions of pollock consumed:
ATF.M = ATF.M %>% left_join(WEP_age_ATF)
PC.M = PC.M %>% left_join(WEP_age_PC)
PH.M = PH.M %>% left_join(WEP_age_PH)
SBL.M = SBL.M %>% left_join(WEP_age_SBL)
WEP.M = WEP.M %>% left_join(WEP_age_WEP)

# Calculate age-specific consumption:
ATF.M$WEP.cons_g_age = round(with(ATF.M,  
                                  WEP.cons_g * propWT.age), digits=0)
PC.M$WEP.cons_g_age = round(with(PC.M,  
                                  WEP.cons_g * propWT.age), digits=0)
PH.M$WEP.cons_g_age = round(with(PH.M,  
                                  WEP.cons_g * propWT.age), digits=0)
SBL.M$WEP.cons_g_age = round(with(SBL.M,  
                                  WEP.cons_g * propWT.age), digits=0)
WEP.M$WEP.cons_g_age = round(with(WEP.M,  
                                  WEP.cons_g * propWT.age), digits=0)

# Save as csv:
GOA.poll_M = ATF.M %>% 
  full_join(PC.M) %>%
  full_join(., PH.M) %>%
  full_join(., SBL.M) %>%
  full_join(., WEP.M) %>%
  rename(WEP_age = age.class, M_g = WEP.cons_g_age) 
GOA.poll_M$M_t = round((GOA.poll_M$M_g / 1000000), digits=2)

# Estimate the number of pollock consumed based on mean weights:
GOA.poll_M$M_N = with(GOA.poll_M, 
              ifelse(WEP_age == "0",  M_g / 11.1,
              ifelse(WEP_age == "1",  M_g / 69.5,
              ifelse(WEP_age == "2",  M_g / 259,
              ifelse(WEP_age == "3+", M_g / 637, NA)))))

GOA.poll_M = GOA.poll_M[,c("Year", "Predator", "WEP_age", "M_t", "M_N")]
  write.csv(GOA.poll_M, "6_PredationIndices/GOApoll_M.csv")

# Figure 1.
  # data = read.csv("GOApoll_M.csv")

# top; age-specific predation (all predators; MT):
M_bio.age = GOA.poll_M %>%
  group_by(Year, WEP_age) %>%
  summarise(mill.MT = sum(M_t / 1000000))

M_bio.age$WEP_age = ordered(M_bio.age$WEP_age, levels = c("3+", "2", "1", "0"))
M_bio.age = M_bio.age %>%
  filter(WEP_age != "0")

M.bio_age = ggplot(M_bio.age, aes(x = Year, y = mill.MT, group = WEP_age, fill = WEP_age)) +
  geom_area(position="stack", show.legend=T, col="gray20", lwd=0.1) +
  scale_fill_manual(values = c("mediumblue", "dodgerblue1", "lightskyblue1", "azure")) +
  ggtitle("") +
    theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(1.0, "lines"),
        legend.position = c(0.89,0.85),
        legend.key.width = unit(4.0, "mm"),
        legend.key.height = unit(5.0, "mm"),
        axis.line = element_line(color="black", linetype="solid"),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(family="Arial", color="white", size=11),
        axis.title.x = element_text(family="Arial", color="white", size=11),
        axis.text.y = element_text(family="Arial", color="black", size=11),
        plot.margin=unit(c(0,0.5,0,0.5),"cm")) +
  scale_x_continuous(breaks=c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019), labels=c("1990","1993","1996","1999","2001","","2005","","2009","","2013","","2017",""), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,7), breaks=c(0,2,4,6)) +
  labs(x = "Survey Year", y = bquote('Predation'~(MT~10^6)), fill = "Age class (yr)")
ggsave(filename="6_PredationIndices/Fig1A.png", plot=M.bio_age, dpi=500, width=6, height=3.5, units="in")

# bottom; age-specific predation (all predators; N):
M_n.age = GOA.poll_M  %>%  
  filter(WEP_age != "0") %>%
  group_by(Year, WEP_age) %>%
  summarise(sum.N = sum(M_N / 1000000000))

M.n_age = ggplot(M_n.age, aes(x = Year, y = sum.N, group = WEP_age, fill = WEP_age)) +
  geom_area(position="stack", show.legend=T, col="gray20", lwd=0.1) +
  scale_fill_manual(values = c("mediumblue", "dodgerblue1", "lightskyblue1", "azure")) +
  ggtitle("") +
    theme(panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(1.0, "lines"),
        legend.position = "none",
        axis.line = element_line(color="black", linetype="solid"),
        axis.ticks = element_line(color="black"),
        axis.title.x = element_text(family="Arial", color="black", size=11, vjust = 0.8),
        axis.text = element_text(family="Arial", color="black", size=11),
        plot.margin=unit(c(0,0.5,0,0.5),"cm")) +
  scale_x_continuous(breaks=c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019), labels=c("1990","1993","1996","1999","2001","","2005","","2009","","2013","","2017",""), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0,0), limits=c(0,40), breaks=c(0,10,20,30,40)) +
  labs(x = "Survey Year", y = bquote('Predation'~(N~10^9)))
ggsave(filename="6_PredationIndices/Fig1B.png", plot=M.n_age, dpi=500, width=6, height=3.5, units="in")
