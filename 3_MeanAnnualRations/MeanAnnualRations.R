# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)
  # Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# Citation: 
# Dorn MW and CL Barnes. In revision. Time-varying predation as a modifier of constant natural mortality for Gulf of Alaska walleye pollock. Fish Res. VSI Natural Mortality. Forthcoming.

# This script includes code necessary to estimate mean annual rations (g per g per yr) for major pollock predators in the Gulf of Alaska; specifically,  Arrowtooth Flounder, Pacific Cod, Pacific Halibut, Sablefish, and Walleye Pollock. Methods are described in Barnes et al. (2020). Mean annual rations were multiplied by total predator biomass, relative predator density, and proportions of pollock consumed to estimate year- and age-specific predation on pollock. Time-varying predation mortality was then tested as a modifier of assumed constant natural mortality in the main stock assessment model.

# Published bioenergetics parameters were used to estimate maximum daily consumption (Cmax; g per g per d), the temperature scaling functions, and mean relative foraging rates (RFR) were obtained from ATF, PC, and WEP: Holsman and Aydin (2015), PH: Holsman et al. (2019), SBL: Harvey (2009) and Armstrong and Schindler (2011). Mean region- and size-specific number of foraging days per year were also specified in Holsman and Aydin (2015). 

# Bottom trawl survey data were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for details about survey design and data collection. 

# References:
# Armstrong JB and Schindler DE. 2011. Excess digestive capacity in predators reflects a life of feast and famine. Nature. 476:84–87.
# Beaudreau, A. H., and T. E. Essington. 2009. Development of a new field-based approach for estimating consumption rates of fishes and comparison with a bioenergetics model for lingcod (Ophiodon elongatus). Canadian Journal of Fisheries and Aquatic Sciences 66:565−578.
# Harvey, C. J. 2009. Effects of temperature change on demersal fisheries in the California Current: a bioenergetics approach. Canadian Journal of Fisheries and Aquatic Sciences 66:1449–1461.
# Holsman, K. K., and K. Aydin. 2015. Comparative methods for evaluating climate change impacts on the foraging ecology of Alaskan groundfish. Marine Ecology Progress Series 521:217–235.
# Holsman, K. K., K. Aydin, J. Sullivan, T. Hurst, and G. Kruse. 2019. Climate effects and bottom-up controls on growth and size-at-age of Pacific halibut (Hippoglossus stenolepis) in Alaska (USA). Fisheries Oceanography 28:345–358.
# von Szalay, P. G., and N. W. Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. NOAA Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/AFSC/GOApollock_ESP/")
require(tidyr)
require(raster)
require(rgeos)
require(rgbif)
require(viridis)
require(gridExtra)
require(rasterVis)
require(purrr)
require(mapproj)
require(devtools)
require(stringr)
require(maptools)
require(rgdal)
require(ggplot2)
require(ggmap)
require(dplyr)

#########################################################
### INITIAL DATA PREPARATION ###

# Read in haul data for the bottom trawl survey:
haul.data =read.csv("1_Data/race_cpue_by_haul.csv") 
haul.data = haul.data %>%
  rename(YEAR = Year, CRUISE = Cruise.Number, VESSEL = Vessel.Number, HAUL = Haul.Number, GEAR_DEPTH = Gear.Depth, GEAR_TEMPERATURE = Gear.Temperature...C.) 

# Read in and format food habits data (all predators):
food.habits = read.csv("1_Data/GOA_RawPP_AllDiets.csv") 

# Join haul information and food habits data:
all.data = haul.data %>% right_join(food.habits, by = c("CRUISE", "VESSEL", "HAUL", "YEAR"))

# Exclude SEAK (GOA main assessment area only):
all.data = subset(all.data, Starting.Longitude..dd. < -140)

# Exclude data prior to 1990, when survey methods were standardized:
all.data = all.data %>%
  rename(Year = YEAR)
all.data = subset(all.data, Year >= 1990)

# Remove hauls without temperature measurements:
data.sub = all.data %>%
  filter(is.finite(GEAR_TEMPERATURE))

##########################################################
### BIOENERGETICS CALCULATIONS ### 

### ATF:
# Select food habits for Arrowtooth Flounder >= 19 cm:
ATF.data = subset(data.sub, PRED_NODC == 8857040102 & PRED_LEN >= 19)

# Parameter estimates from Holsman and Aydin (2015):
Ca = 0.125
Cb = -0.199
Cq = 2.497
Tco = 20.512
Tcm = 26

Y = (log(Cq)) * (Tcm - Tco + 2)
Z = (log(Cq)) * (Tcm - Tco)
X = (Z^2 * (1 + (1 + 40/Y)^0.5)^2)/400
ATF.data$V = (Tcm - ATF.data$GEAR_TEMPERATURE)/(Tcm - Tco)

# Compute temperature scaling function:
ATF.data$fT = with(ATF.data, V^X * (exp(X*(1-V))))

# Calculate maximum daily consumption (Cmax; g/g/d):
ATF.data$Cmax_ggd = Ca * (ATF.data$PRED_WT^Cb) * ATF.data$fT

# Assume fish are feeding at their theoretical maximum (RFR = 1): 
  # RFR, juv(< 40 cm): 0.79 # Holsman and Aydin 2015
  # RFR, adult(>= 40 cm): 1.07 # Holsman and Aydin 2015

# Multiply Cmax by effective number of foraging days (EFD; Holsman and Aydin 2015) to scale to year:
  # EFD ATF <  40 cm = 346 d
  # EFD ATF >= 40 cm = 306 d
ATF.data$EFD = ifelse(ATF.data$PRED_LEN < 40, 346, 306)
ATF.data$Cmax_ggy = ATF.data$Cmax_ggd * ATF.data$EFD
ATF.data = subset(ATF.data, Cmax_ggy !="Inf")

  saveRDS(ATF.data, "3_MeanAnnualRations/Rations_ATF.rds")
# readRDS("3_MeanAnnualRations/Rations_ATF.rds")


### PC:
# Select food habits for Pacific Cod >= 0 cm:
PC.data = subset(data.sub, PRED_NODC == 8791030401 & PRED_LEN >= 0)

# Parameter estimates from Holsman and Aydin (2015) and Holsman et al. (in prep):
Ca = 0.035
Cb = -0.122
Cq = 3.079
Tco = 10.957
Tcm = 25.901

Y = (log(Cq)) * (Tcm - Tco + 2)
Z = (log(Cq)) * (Tcm - Tco)
X = (Z^2 * (1 + (1 + 40/Y)^0.5)^2)/400
PC.data$V = (Tcm - PC.data$GEAR_TEMPERATURE)/(Tcm - Tco)  

# Compute temperature scaling function:
PC.data$fT = with(PC.data, V^X * (exp(X*(1-V))))

# Calculate maximum daily consumption (Cmax; g/g/d):
PC.data$Cmax_ggd = Ca * (PC.data$PRED_WT^Cb) * PC.data$fT

# Assume fish are feeding at their theoretical maximum (RFR = 1): 
  # RFR, juv(< 55 cm): 0.41 # Holsman and Aydin 2015
  # RFR, adult(>= 55 cm): 0.47 # Holsman and Aydin 2015

# Multiply Cmax by effective number of foraging days (EFD; Holsman and Aydin 2015) to scale to year:
  # EFD PC <  55 cm = 365 d
  # EFD PC >= 55 cm = 329 d
PC.data$EFD = ifelse(PC.data$PRED_LEN < 55, 365, 329)
PC.data$Cmax_ggy = PC.data$Cmax_ggd * PC.data$EFD
PC.data = subset(PC.data, Cmax_ggy !="Inf")

  saveRDS(PC.data, "3_MeanAnnualRations/Rations_PC.rds")
# readRDS("3_MeanAnnualRations/Rations_PC.rds")

### PH:  
# Select food habits for Pacific Halibut >= 82 cm:
PH.data = subset(data.sub, PRED_NODC == 8857041901 & PRED_LEN >= 82)

# Parameter estimates from Holsman et al. (2019):
Ca = 0.0625
Cb = -0.1076
Cq = 3.084 # listed as Qc in Holsman et al. (2018)
Tco = 12.97
Tcm = 18

Y = (log(Cq)) * (Tcm - Tco + 2)
Z = (log(Cq)) * (Tcm - Tco)
X = (Z^2 * (1 + (1 + 40/Y)^0.5)^2)/400
PH.data$V = (Tcm - PH.data$GEAR_TEMPERATURE)/(Tcm - Tco)

# Compute temperature scaling function:
PH.data$fT = with(PH.data, V^X * (exp(X*(1-V))))

# Calculate maximum daily consumption (Cmax; g/g/d):
PH.data$Cmax_ggd = Ca * (PH.data$PRED_WT^Cb) * PH.data$fT

# Assume fish are feeding at their theoretical maximum (RFR = 1): 
  # RFR, juv (< 40 cm): 2C = 0.02; 3A = 0.46; 3B = 0.37; 4AS = 0.19 # Holsman et al. 2018
  # RFR, adults (40 to 120 cm): 2C = 0.33; 3A = 0.41; 3B = 0.53; 4AS = 0.33 # Holsman et al. 2018

# Multiply Cmax by effective number of foraging days (EFD; Holsman et al. 2019) to scale to year:
  # EFD PH assumed to be 365 d
PH.data$EFD = 365 
PH.data$Cmax_ggy = PH.data$Cmax_ggd * PH.data$EFD
PH.data = subset(PH.data, Cmax_ggy !="Inf")

  saveRDS(PH.data, "3_MeanAnnualRations/Rations_PH.rds")
# readRDS("3_MeanAnnualRations/Rations_PH.rds")

### SBL:
# Select food habits for Sablefish >= 45 cm:
SBL.data = subset(data.sub, PRED_NODC == 8827020101 & PRED_LEN >= 45)

# Parameter estimates from Harvey (2009):
Ca = 0.4200
Cb = -0.3300
Cq = 2.20
Tco = 18 # listed as CTO in Harvey (2009)
Tcm = 23 # listed as CTM in Harvey (2009)

Y = (log(Cq)) * (Tcm - Tco + 2)
Z = (log(Cq)) * (Tcm - Tco)
X = (Z^2 * (1 + (1 + 40/Y)^0.5)^2)/400
SBL.data$V = (Tcm - SBL.data$GEAR_TEMPERATURE)/(Tcm - Tco)

# Compute temperature scaling function:
SBL.data$fT = with(SBL.data, V^X * (exp(X*(1-V))))

# Calculate maximum daily consumption (Cmax; g/g/d):
SBL.data$Cmax_ggd = Ca * (SBL.data$PRED_WT^Cb) * SBL.data$fT

# Assume fish are feeding at their theoretical maximum (RFR = 1): 
  # RFR, juveniles (40-50 cm): 0.273 # Armstrong and Schindler 2011
  # RFR, adults (>= 50 cm): 0.259 # Armstrong and Schindler 2011

# Multiply Cmax by effective number of foraging days to scale to year:
SBL.data$EFD = 365 # assumed (no reference)
SBL.data$Cmax_ggy = SBL.data$Cmax_ggd * SBL.data$EFD
SBL.data = subset(SBL.data, Cmax_ggy !="Inf")

  saveRDS(SBL.data, "3_MeanAnnualRations/Rations_SBL.rds")
# readRDS("3_MeanAnnualRations/Rations_SBL.rds")

### WEP:
# Select food habits for Walleye Pollock >= 37 cm:
WEP.data = subset(data.sub, PRED_NODC == 8791030701 & PRED_LEN >= 37)

# Parameter estimates from Holsman and Aydin (2015):
Ca = 0.119
Cb = -0.46
Cq = 2.6
Tco = 10
Tcm = 15

Y = (log(Cq)) * (Tcm - Tco + 2)
Z = (log(Cq)) * (Tcm - Tco)
X = (Z^2 * (1 + (1 + 40/Y)^0.5)^2)/400
WEP.data$V = (Tcm - WEP.data$GEAR_TEMPERATURE)/(Tcm - Tco) 

# Compute temperature scaling function:
WEP.data$fT = with(WEP.data, V^X * (exp(X*(1-V))))

# Calculate maximum daily consumption (Cmax; g/g/d):
WEP.data$Cmax_ggd = Ca * (WEP.data$PRED_WT^Cb) * WEP.data$fT

# Assume fish are feeding at their theoretical maximum (RFR = 1): 
  # RFR, juveniles (< 40 cm): 0.49 # Holsman and Aydin 2015
  # RFR, adults (>= 40 cm): 0.56 # Holsman and Aydin 2015

# Multiply Cmax by effective number of foraging days (EFD; Holsman and Aydin 2015) to scale to year:
WEP.data$EFD = 365 
WEP.data$Cmax_ggy = WEP.data$Cmax_ggd * WEP.data$EFD
WEP.data = subset(WEP.data, Cmax_ggy !="Inf")

  saveRDS(WEP.data, "3_MeanAnnualRations/Rations_WEP.rds")
# readRDS("3_MeanAnnualRations/Rations_WEP.rds")
  
# check:
ATF.plot = ATF.data %>%
  group_by(Year) %>%
  summarise(mean.Cmax = mean(Cmax_ggy))
ggplot(ATF.plot, aes(x = Year, y = mean.Cmax)) +
  geom_line(lwd=1) +
  geom_point(size = 3.5, stroke=2, shape=3, alpha = 1) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "none") +
  labs(x="", y="Predation on Walleye Pollock (mill MT)", fill="Age Class (yr)") +
  scale_x_continuous(limits=c(1990,2019.5), breaks = c(1990,1996,2001,2005,2009,2013,2017), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0.01,0.01), limits=c(0,8), breaks=c(0,2,4,6,8))

PC.plot = PC.data %>%
  group_by(Year) %>%
  summarise(mean.Cmax = mean(Cmax_ggy))
ggplot(PC.plot, aes(x = Year, y = mean.Cmax)) +
  geom_line(lwd=1) +
  geom_point(size = 3.5, stroke=2, shape=3, alpha = 1) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "none") +
  labs(x="", y="Predation on Walleye Pollock (mill MT)", fill="Age Class (yr)") +
  scale_x_continuous(limits=c(1990,2019.5), breaks = c(1990,1996,2001,2005,2009,2013,2017), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0.01,0.01), limits=c(0,8), breaks=c(0,2,4,6,8))

PH.plot = PH.data %>%
  group_by(Year) %>%
  summarise(mean.Cmax = mean(Cmax_ggy))
ggplot(PH.plot, aes(x = Year, y = mean.Cmax)) +
  geom_line(lwd=1) +
  geom_point(size = 3.5, stroke=2, shape=3, alpha = 1) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "none") +
  labs(x="", y="Predation on Walleye Pollock (mill MT)", fill="Age Class (yr)") +
  scale_x_continuous(limits=c(1990,2019.5), breaks = c(1990,1996,2001,2005,2009,2013,2017), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0.01,0.01), limits=c(0,8), breaks=c(0,2,4,6,8))

SBL.plot = SBL.data %>%
  group_by(Year) %>%
  summarise(mean.Cmax = mean(Cmax_ggy)) %>%
  as.data.frame()

# Add mean estimates for years without data:
SBL.plot$mean.Cmax = as.numeric(SBL.plot$mean.Cmax)
  SBL.plot[11,] = c(2013, mean(SBL.plot$mean.Cmax)) 
  SBL.plot[12,] = c(2015, mean(SBL.plot$mean.Cmax)) 
  SBL.plot[13,] = c(2017, mean(SBL.plot$mean.Cmax))  
  SBL.plot[14,] = c(2019, mean(SBL.plot$mean.Cmax))  

ggplot(SBL.plot, aes(x = Year, y = mean.Cmax)) +
  geom_line(lwd=1) +
  geom_point(size = 3.5, stroke=2, shape=3, alpha = 1) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "none") +
  labs(x="", y="Predation on Walleye Pollock (mill MT)", fill="Age Class (yr)") +
  scale_x_continuous(limits=c(1990,2019.5), breaks = c(1990,1996,2001,2005,2009,2013,2017), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0.01,0.01), limits=c(0,8), breaks=c(0,2,4,6,8))

WEP.plot = WEP.data %>%
  group_by(Year) %>%
  summarise(mean.Cmax = mean(Cmax_ggy))
ggplot(WEP.plot, aes(x = Year, y = mean.Cmax)) +
  geom_line(lwd=1) +
  geom_point(size = 3.5, stroke=2, shape=3, alpha = 1) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "none") +
  labs(x="", y="Predation on Walleye Pollock (mill MT)", fill="Age Class (yr)") +
  scale_x_continuous(limits=c(1990,2019.5), breaks = c(1990,1996,2001,2005,2009,2013,2017), expand = c(0.01, 0.01)) +
  scale_y_continuous(expand=c(0.01,0.01), limits=c(0,8), breaks=c(0,2,4,6,8))