# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)
  # Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# Citation: 
# Dorn MW and CL Barnes. In revision. Time-varying predation as a modifier of constant natural mortality for Gulf of Alaska walleye pollock. Fish Res. VSI Natural Mortality. Forthcoming.

# This script includes code necessary to calculate fork length- and biomass-weighted proportions of pollock (Chipps and Garvey 2007) consumed by major predators in the Gulf of Alaska; specifically,  Arrowtooth Flounder, Pacific Cod, Pacific Halibut, Sablefish, and Walleye Pollock. Methods are described in Barnes et al. (2020). Proportions of pollock consumed were multiplied by total predator biomass, relative predator density, and mean annual rations to estimate year- and age-specific predation on pollock. Time-varying predation mortality was then tested as a modifier of assumed constant natural mortality in the main stock assessment model.

# Bottom trawl survey data were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for details about survey design and data collection. 

# Food habits data were obtained from the Resource Ecology and Ecosystem Modeling (REEM) program of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See Livingston et al. (2017) for methods. 

# References:
# Chipps SR and JE Garvey. 2007. Assessment of diets and feeding patterns. In: Analysis and interpretation of freshwater fisheries data. CS Guy and ML Brown, eds. Bethesda, MD. Amer Fish Soc: 473–514.
# Livingston PA, K Aydin, TW Buckley, GM Lang, M-S Yang, and BS Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fish. 100(4):443–470.
# von Szalay, PG and NW Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. NOAA Techn Mem NMFS-AFSC-325. 

setwd("~/Documents/AFSC/GOApollock_ESP/")
require(tidyr)
require(dplyr)
require(sf)
require(ggplot2)

#########################################################
### INITIAL DATA PREPARATION ###

# Read in and format predator length data from AFSC bottom trawls:
lengths = read.csv(unz("1_Data/race_length_by_haul.csv.zip", "race_length_by_haul.csv"), header=T, skip=7, na.strings=c("","NA"))

# Create a unique haul identifier by concatenating vessel, cruise, and haul:
lengths$Haul.Number_Join = with(lengths, paste(Vessel.Number, Cruise.Number, Haul.Number, sep=""))

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
lengths = subset(lengths, Year >= 1990)

# Convert predator lengths to cm and bin:
lengths$PredLength = as.numeric(as.character(lengths$Length..mm.)) / 10
lengths$Frequency = as.numeric(as.character(lengths$Frequency))

### ATF:
ATF.lengths = subset(lengths, Common.Name == "arrowtooth flounder")
ATF.lengths$FL.bin = cut(ATF.lengths$PredLength, breaks = c(0,18,28,38,48,58,68,78,88,100))
  levels(ATF.lengths$FL.bin) = c("<=18", "19-28", "29-38", "39-48", "49-58", "59-68", "69-78", "79-88", ">88")

### PC:
PC.lengths = subset(lengths, Common.Name == "Pacific cod") 
PC.lengths$FL.bin = cut(PC.lengths$PredLength, breaks = c(0,10,20,30,40,50,60,70,80,90,100,110))
  levels(PC.lengths$FL.bin) = c("<=10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100", "101-110")

### PH:
PH.lengths = subset(lengths, Common.Name == "Pacific halibut")
PH.lengths$FL.bin = cut(PH.lengths$PredLength, breaks = c(0,21,41,61,81,101,121,141,161,181,201,211))
  levels(PH.lengths$FL.bin) = c("<=21", "22-41", "42-61", "62-81", "82-101", "102-121", "122-141", "142-161", "162-181", "182-201", ">=202")

### SBL:
SBL.lengths = subset(lengths, Common.Name == "sablefish")
SBL.lengths$FL.bin = cut(SBL.lengths$PredLength, breaks = c(0,14,24,34,44,54,64,74,84,94,110))
  levels(SBL.lengths$FL.bin) = c("<=14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85-94", ">95")

### WEP:
WEP.lengths = subset(lengths, Common.Name == "walleye pollock")
WEP.lengths$FL.bin = cut(WEP.lengths$PredLength, breaks = c(0,16,26,36,46,56,66,76,86,99))
  levels(WEP.lengths$FL.bin) = c("<=16", "17-26", "27-36", "37-46", "47-56", "57-66", "67-76", "77-86", ">=97")
  
##########################################################
### SUMMARIZE FOOD HABITS DATA ###
# Read in and format food habits data from subsampled stomachs:
food.habs = read.csv("1_Data/GOA_RawPP_AllDiets.csv") 
  
# Create a unique haul identifier by concatenating vessel, cruise, and haul:
food.habs$Haul.Number_Join = with(food.habs, paste(VESSEL, CRUISE, HAUL, sep=""))
food.habs = food.habs %>%
  rename(PredLength = PRED_LEN)

# Link to survey data for lon, lat information:
food.habits = food.habs %>% 
  left_join(lengths[,c("Haul.Number_Join", "Starting.Longitude..dd.", "Starting.Latitude..dd.")])
  
# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
food.habits = food.habits %>%
  rename(Year = YEAR)
food.habits = subset(food.habits, Year >= 1990)

# Remove empty stomachs:
pos.prey = subset(food.habits, PREY_NODC != "0")

# Create 10-cm predator length bins:
### ATF:
ATF.prey = subset(pos.prey, PRED_NODC == 8857040102)
ATF.prey$FL.bin = cut(ATF.prey$PredLength, breaks = c(0,18,28,38,48,58,68,78,88,100))
  levels(ATF.prey$FL.bin) = c("<=18", "19-28", "29-38", "39-48", "49-58", "59-68", "69-78", "79-88", ">88")

### PC:
PC.prey = subset(pos.prey, PRED_NODC == 8791030401)
PC.prey$FL.bin = cut(PC.prey$PredLength, breaks = c(0,10,20,30,40,50,60,70,80,90,100,110))
  levels(PC.prey$FL.bin) = c("<=10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100", "101-110")

### PH:
PH.prey = subset(pos.prey, PRED_NODC == 8857041901)
PH.prey$FL.bin = cut(PH.prey$PredLength, breaks = c(0,21,41,61,81,101,121,141,161,181,201,211))
  levels(PH.prey$FL.bin) = c("<=21", "22-41", "42-61", "62-81", "82-101", "102-121", "122-141", "142-161", "162-181", "182-201", ">=202")

### SBL:
SBL.prey = subset(pos.prey, PRED_NODC == 8827020101)
SBL.prey$FL.bin = cut(SBL.prey$PredLength, breaks = c(0,14,24,34,44,54,64,74,84,94,110))
levels(SBL.prey$FL.bin) = c("<=14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85-94", ">95")

### WEP:
WEP.prey = subset(pos.prey, PRED_NODC == 8791030701)
WEP.prey$FL.bin = cut(WEP.prey$PredLength, breaks = c(0,16,26,36,46,56,66,76,86,99))
  levels(WEP.prey$FL.bin) = c("<=16", "17-26", "27-36", "37-46", "47-56", "57-66", "67-76", "77-86", ">=97")

##########################################################
### CALCULATE WEIGHTED PROPORTIONS OF POLLOCK CONSUMED ###
  
### ATF:
# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (food habits):
ATF.sizes_diets = ATF.prey %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_diets = length(PredLength))

ATF.props_diets = ATF.sizes_diets %>%
  group_by(Haul.Number_Join) %>%
  mutate(length.prop_diets = length.freq_diets / sum(length.freq_diets)) 

# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (bottom trawl survey):
ATF.sizes_trawl = ATF.lengths %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_trawl = sum(Frequency))

ATF.props_trawl = ATF.sizes_trawl %>%
  group_by(Haul.Number_Join) %>% 
  mutate(length.prop_trawl = length.freq_trawl / sum(length.freq_trawl))

# Calculate sample weights from dividing length-based proportions of fish caught by those subsampled for gut content analysis in each survey year and grid cell:
ATF.lengths.prop = ATF.props_trawl %>% right_join(ATF.props_diets)
  ATF.lengths.prop[is.na(ATF.lengths.prop)] = 0
ATF.lengths.prop$length.wt = with(ATF.lengths.prop, (length.prop_trawl / length.prop_diets))

# Account for length-structured stomach sampling:
  # Multiply prey mass (g) by sample weights according to predator length frequencies:
ATF.prey.FL = ATF.prey %>% left_join(ATF.lengths.prop)
  ATF.prey.FL$prey_FL.wt = ATF.prey.FL$PREY_TWT * ATF.prey.FL$length.wt
ATF.prey.FL = subset(ATF.prey.FL, !is.na(Starting.Longitude..dd.))

# Account for unequal sampling across the Gulf of Alaska
  # Multiply prey mass (g) by sample weights according to grid-specific biomass:
load("2_RelativePredatorDensities/normATFabun.rda")
  normAbun_ATF = as.data.frame(subset(normAbun_ATF, select = -geometry))

ATF.prey.Bio = normAbun_ATF %>%
  group_by(Year) %>% 
  mutate(prey_Bio.wt = ATF.predAbun / mean(ATF.predAbun))
ATF.prey.Bio = ATF.prey.Bio %>%
  as.data.frame() %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
         crs = 4326)

# ATF.prey.Bio = unique(ATF.prey.Bio[,c("id2", "Year", "ATF.prey.Bio")])
ATF.prey.FL$Year = as.factor(ATF.prey.FL$Year)
ATF.prey_FL = ATF.prey.FL %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
           crs = 4326)
 
ATF.prey_weighted = list()
for(i in unique(ATF.prey_FL$Year)) {
  diet.yr = subset(ATF.prey_FL, Year == i)
  bio.yr = subset(ATF.prey.Bio, Year == i)
  ATF.prey_weighted[[i]] = st_join(diet.yr, bio.yr,
                                   join = st_nearest_feature, left = T) }
ATF.prey_weighted.df = as.data.frame(dplyr::bind_rows(ATF.prey_weighted))

ATF.prey_weighted.df$prey.WT = with(ATF.prey_weighted.df,  
                                    (prey_FL.wt * prey_Bio.wt))

# Select only sizes of fish encompassed in total predator biomass estimates:
ATF.prey_weighted.df = subset(ATF.prey_weighted.df, PredLength >= 19)

# Categorize prey (pollock vs other):  
ATF.prey_weighted.df$poll = with(ATF.prey_weighted.df, 
              ifelse(PREY_NODC == 8791030701, "pollock", "other"))
ATF.prey_weighted.df$poll = ordered(ATF.prey_weighted.df$poll, levels = c("pollock", "other"))
ATF.prey_weighted.df$Year = as.numeric(as.character(ATF.prey_weighted.df$Year.x))

ATF.coords = st_coordinates(ATF.prey_weighted.df$geometry)
  ATF.prey_weighted.assess = cbind(ATF.prey_weighted.df, ATF.coords) %>%
    rename(Starting.Longitude..dd. = X, Starting.Latitude..dd. = Y) %>%
    filter(Starting.Longitude..dd. < -140)

ATF.prey_weighted.assess = unique(ATF.prey_weighted.assess)  
  save(ATF.prey_weighted.assess, 
    file = "4_ProportionsPollockConsumed/propWEP_ATF_19.rda")
# load("4_ProportionsPollockConsumed/propWEP_ATF_19.rda")
  
  
### PC:
# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (food habits):
PC.sizes_diets = PC.prey %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_diets = length(PredLength))

PC.props_diets = PC.sizes_diets %>%
  group_by(Haul.Number_Join) %>%
  mutate(length.prop_diets = length.freq_diets / sum(length.freq_diets)) 

# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (bottom trawl survey):
PC.sizes_trawl = PC.lengths %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_trawl = sum(Frequency))

PC.props_trawl = PC.sizes_trawl %>%
  group_by(Haul.Number_Join) %>% 
  mutate(length.prop_trawl = length.freq_trawl / sum(length.freq_trawl))

# Calculate sample weights from dividing length-based proportions of fish caught by those subsampled for gut content analysis in each survey year and grid cell:
PC.lengths.prop = PC.props_trawl %>% right_join(PC.props_diets)
  PC.lengths.prop[is.na(PC.lengths.prop)] = 0
PC.lengths.prop$length.wt = with(PC.lengths.prop, 
                                 (length.prop_trawl / length.prop_diets))

# Account for length-structured stomach sampling:
  # Multiply prey mass (g) by sample weights according to predator length frequencies:
PC.prey.FL = PC.prey %>% left_join(PC.lengths.prop)
  PC.prey.FL$prey_FL.wt = PC.prey.FL$PREY_TWT * PC.prey.FL$length.wt
PC.prey.FL = subset(PC.prey.FL, !is.na(Starting.Longitude..dd.))

# Account for unequal sampling across the Gulf of Alaska
  # Multiply prey mass (g) by sample weights according to grid-specific biomass:
load("2_RelativePredatorDensities/normPCabun.rda")
  normAbun_PC = as.data.frame(subset(normAbun_PC, select = -geometry))

PC.prey.Bio = normAbun_PC %>%
  group_by(Year) %>% 
  mutate(prey_Bio.wt = PC.predAbun / mean(PC.predAbun))
PC.prey.Bio = PC.prey.Bio %>%
  as.data.frame() %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
         crs = 4326)

# PC.prey.Bio = unique(PC.prey.Bio[,c("id2", "Year", "PC.prey.Bio")])
PC.prey.FL$Year = as.factor(PC.prey.FL$Year)
PC.prey_FL = PC.prey.FL %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
           crs = 4326)
 
PC.prey_weighted = list()
for(i in unique(PC.prey_FL$Year)) {
  diet.yr = subset(PC.prey_FL, Year == i)
  bio.yr = subset(PC.prey.Bio, Year == i)
  PC.prey_weighted[[i]] = st_join(diet.yr, bio.yr,
                                   join = st_nearest_feature, left = T) }
PC.prey_weighted.df = as.data.frame(dplyr::bind_rows(PC.prey_weighted))
  PC.prey_weighted.df$prey.WT = with(PC.prey_weighted.df, (prey_FL.wt * prey_Bio.wt))

# Select only sizes of fish encompassed in total predator biomass estimates:
PC.prey_weighted.df = subset(PC.prey_weighted.df, PredLength >= 0)
  
# Categorize prey (pollock vs other):  
PC.prey_weighted.df$poll = with(PC.prey_weighted.df, 
              ifelse(PREY_NODC == 8791030701, "pollock", "other"))
PC.prey_weighted.df$poll = ordered(PC.prey_weighted.df$poll, levels = c("pollock", "other"))
PC.prey_weighted.df$Year = as.numeric(as.character(PC.prey_weighted.df$Year.x))

PC.coords = st_coordinates(PC.prey_weighted.df$geometry)
  PC.prey_weighted.assess = cbind(PC.prey_weighted.df, PC.coords) %>%
    rename(Starting.Longitude..dd. = X, Starting.Latitude..dd. = Y) %>%
    filter(Starting.Longitude..dd. < -140)

PC.prey_weighted.assess = unique(PC.prey_weighted.assess)  
  save(PC.prey_weighted.assess, 
    file = "4_ProportionsPollockConsumed/propWEP_PC_0.rda")
# load("4_ProportionsPollockConsumed/propWEP_PC_0.rda")
  
  
### PH:
# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (food habits):
PH.sizes_diets = PH.prey %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_diets = length(PredLength))

PH.props_diets = PH.sizes_diets %>%
  group_by(Haul.Number_Join) %>%
  mutate(length.prop_diets = length.freq_diets / sum(length.freq_diets)) 

# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (bottom trawl survey):
PH.sizes_trawl = PH.lengths %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_trawl = sum(Frequency))

PH.props_trawl = PH.sizes_trawl %>%
  group_by(Haul.Number_Join) %>% 
  mutate(length.prop_trawl = length.freq_trawl / sum(length.freq_trawl))

# Calculate sample weights from dividing length-based proportions of fish caught by those subsampled for gut content analysis in each survey year and grid cell:
PH.lengths.prop = PH.props_trawl %>% right_join(PH.props_diets)
  PH.lengths.prop[is.na(PH.lengths.prop)] = 0
PH.lengths.prop$length.wt = with(PH.lengths.prop, 
                                 (length.prop_trawl / length.prop_diets))

# Account for length-structured stomach sampling:
  # Multiply prey mass (g) by sample weights according to predator length frequencies:
PH.prey.FL = PH.prey %>% left_join(PH.lengths.prop)
  PH.prey.FL$prey_FL.wt = PH.prey.FL$PREY_TWT * PH.prey.FL$length.wt
PH.prey.FL = subset(PH.prey.FL, !is.na(Starting.Longitude..dd.))

# Account for unequal sampling across the Gulf of Alaska
  # Multiply prey mass (g) by sample weights according to grid-specific biomass:
load("2_RelativePredatorDensities/normPHabun.rda")

PH.prey.Bio = normAbun_PH %>%
  group_by(Year) %>% 
  mutate(prey_Bio.wt = PH.predAbun / mean(PH.predAbun))

PH.prey.Bio = PH.prey.Bio %>%
  as.data.frame() %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
         crs = 4326)

PH.prey.FL$Year = as.factor(PH.prey.FL$Year)
PH.prey_FL = PH.prey.FL %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
           crs = 4326)
 
PH.prey_weighted = list()
for(i in unique(PH.prey_FL$Year)) {
  diet.yr = subset(PH.prey_FL, Year == i)
  bio.yr = subset(PH.prey.Bio, Year == i)
  PH.prey_weighted[[i]] = st_join(diet.yr, bio.yr,
                                   join = st_nearest_feature, left = T) }
PH.prey_weighted.df = dplyr::bind_rows(PH.prey_weighted)
PH.prey_weighted.df$prey.WT = with(PH.prey_weighted.df, (prey_FL.wt * prey_Bio.wt))

# Select only sizes of fish encompassed in total predator biomass estimates:
PH.prey_weighted.df = subset(PH.prey_weighted.df, PredLength >= 82)

# Categorize prey (pollock vs other):  
PH.prey_weighted.df$poll = with(PH.prey_weighted.df, 
              ifelse(PREY_NODC == 8791030701, "pollock", "other"))
PH.prey_weighted.df$poll = ordered(PH.prey_weighted.df$poll, levels = c("pollock", "other"))
PH.prey_weighted.df$Year = as.numeric(as.character(PH.prey_weighted.df$Year.x))

PH.coords = st_coordinates(PH.prey_weighted.df$geometry)
  PH.prey_weighted.assess = cbind(PH.prey_weighted.df, PH.coords) %>%
    rename(Starting.Longitude..dd. = X, Starting.Latitude..dd. = Y) %>%
    filter(Starting.Longitude..dd. < -140) %>%
    as.data.frame()

PH.prey_weighted.assess = unique(PH.prey_weighted.assess)  
  save(PH.prey_weighted.assess, 
    file = "4_ProportionsPollockConsumed/propWEP_PH_82.rda")
# load("4_ProportionsPollockConsumed/propWEP_PH_82.rda")

  
### SBL:
# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (food habits):
SBL.sizes_diets = SBL.prey %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_diets = length(PredLength))

SBL.props_diets = SBL.sizes_diets %>%
  group_by(Haul.Number_Join) %>%
  mutate(length.prop_diets = length.freq_diets / sum(length.freq_diets)) 

# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (bottom trawl survey):
SBL.sizes_trawl = SBL.lengths %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_trawl = sum(Frequency))

SBL.props_trawl = SBL.sizes_trawl %>%
  group_by(Haul.Number_Join) %>% 
  mutate(length.prop_trawl = length.freq_trawl / sum(length.freq_trawl))

# Calculate sample weights from dividing length-based proportions of fish caught by those subsampled for gut content analysis in each survey year and grid cell:
SBL.lengths.prop = SBL.props_trawl %>% right_join(SBL.props_diets)
  SBL.lengths.prop[is.na(SBL.lengths.prop)] = 0
SBL.lengths.prop$length.wt = with(SBL.lengths.prop, 
                                  (length.prop_trawl / length.prop_diets))

# Account for length-structured stomach sampling:
  # Multiply prey mass (g) by sample weights according to predator length frequencies:
SBL.prey.FL = SBL.prey %>% left_join(SBL.lengths.prop)
  SBL.prey.FL$prey_FL.wt = SBL.prey.FL$PREY_TWT * SBL.prey.FL$length.wt
SBL.prey.FL = subset(SBL.prey.FL, !is.na(Starting.Longitude..dd.))

# Account for unequal sampling across the Gulf of Alaska
  # Multiply prey mass (g) by sample weights according to grid-specific biomass:
load("2_RelativePredatorDensities/normSBLabun.rda")

SBL.prey.Bio = normAbun_SBL %>%
  group_by(Year) %>% 
  mutate(prey_Bio.wt = SBL.predAbun / mean(SBL.predAbun))
SBL.prey.Bio = SBL.prey.Bio %>%
  as.data.frame() %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
         crs = 4326)

# SBL.prey.Bio = unique(SBL.prey.Bio[,c("id2", "Year", "SBL.prey.Bio")])
SBL.prey.FL$Year = as.factor(SBL.prey.FL$Year)
SBL.prey_FL = SBL.prey.FL %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
           crs = 4326)
 
SBL.prey_weighted = list()
for(i in unique(SBL.prey_FL$Year)) {
  diet.yr = subset(SBL.prey_FL, Year == i)
  bio.yr = subset(SBL.prey.Bio, Year == i)
  SBL.prey_weighted[[i]] = st_join(diet.yr, bio.yr,
                                   join = st_nearest_feature, left = T) }
SBL.prey_weighted.df = as.data.frame(dplyr::bind_rows(SBL.prey_weighted))
  SBL.prey_weighted.df$prey.WT = with(SBL.prey_weighted.df, (prey_FL.wt * prey_Bio.wt))

# Select only sizes of fish encompassed in total predator biomass estimates:
SBL.prey_weighted.df = subset(SBL.prey_weighted.df, PredLength >= 45)

# Categorize prey (pollock vs other):  
SBL.prey_weighted.df$poll = with(SBL.prey_weighted.df, 
              ifelse(PREY_NODC == 8791030701, "pollock", "other"))
SBL.prey_weighted.df$poll = ordered(SBL.prey_weighted.df$poll, levels = c("pollock", "other"))
SBL.prey_weighted.df$Year = as.numeric(as.character(SBL.prey_weighted.df$Year.x))

SBL.coords = st_coordinates(SBL.prey_weighted.df$geometry)
  SBL.prey_weighted.assess = cbind(SBL.prey_weighted.df, SBL.coords) %>%
    rename(Starting.Longitude..dd. = X, Starting.Latitude..dd. = Y) %>%
    filter(Starting.Longitude..dd. < -140)
  
# Missing 2013 onward - to fill with mean value as part of index:
SBL.prey_weighted.assess = unique(SBL.prey_weighted.assess) 
  save(SBL.prey_weighted.assess, 
    file = "4_ProportionsPollockConsumed/propWEP_SBL_45.rda")
# load("4_ProportionsPollockConsumed/propWEP_SBL_45.rda")
  
### WEP:
# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (food habits):
WEP.sizes_diets = WEP.prey %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_diets = length(PredLength))

WEP.props_diets = WEP.sizes_diets %>%
  group_by(Haul.Number_Join) %>%
  mutate(length.prop_diets = length.freq_diets / sum(length.freq_diets)) 

# Calculate proportions of fish subsampled within each size bin, survey year, and grid cell (bottom trawl survey):
WEP.sizes_trawl = WEP.lengths %>%
  group_by(Haul.Number_Join, FL.bin) %>% 
  summarise(length.freq_trawl = sum(Frequency))

WEP.props_trawl = WEP.sizes_trawl %>%
  group_by(Haul.Number_Join) %>% 
  mutate(length.prop_trawl = length.freq_trawl / sum(length.freq_trawl))

# Calculate sample weights from dividing length-based proportions of fish caught by those subsampled for gut content analysis in each survey year and grid cell:
WEP.lengths.prop = WEP.props_trawl %>% right_join(WEP.props_diets)
  WEP.lengths.prop[is.na(WEP.lengths.prop)] = 0
WEP.lengths.prop$length.wt = with(WEP.lengths.prop, (length.prop_trawl / length.prop_diets))

# Account for length-structured stomach sampling:
  # Multiply prey mass (g) by sample weights according to predator length frequencies:
WEP.prey.FL = WEP.prey %>% left_join(WEP.lengths.prop)
  WEP.prey.FL$prey_FL.wt = WEP.prey.FL$PREY_TWT * WEP.prey.FL$length.wt
WEP.prey.FL = subset(WEP.prey.FL, !is.na(Starting.Longitude..dd.))
  WEP.prey.FL$length.wt = with(WEP.prey.FL, ifelse(!is.finite(length.wt), 0, length.wt))
  
# Account for unequal sampling across the Gulf of Alaska
  # Multiply prey mass (g) by sample weights according to grid-specific biomass:
load("2_RelativePredatorDensities/normWEPabun.rda")
  normAbun_WEP = as.data.frame(subset(normAbun_WEP, select = -geometry))

WEP.prey.Bio = normAbun_WEP %>%
  group_by(Year) %>% 
  mutate(prey_Bio.wt = WEP.predAbun / mean(WEP.predAbun))
WEP.prey.Bio = WEP.prey.Bio %>%
  as.data.frame() %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
         crs = 4326)

# WEP.prey.Bio = unique(WEP.prey.Bio[,c("id2", "Year", "WEP.prey.Bio")])
WEP.prey.FL$Year = as.factor(WEP.prey.FL$Year)
WEP.prey_FL = WEP.prey.FL %>%
  st_as_sf(coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
           crs = 4326)
 
WEP.prey_weighted = list()
for(i in unique(WEP.prey_FL$Year)) {
  diet.yr = subset(WEP.prey_FL, Year == i)
  bio.yr = subset(WEP.prey.Bio, Year == i)
  WEP.prey_weighted[[i]] = st_join(diet.yr, bio.yr,
                                   join = st_nearest_feature, left = T) }
WEP.prey_weighted.df = as.data.frame(dplyr::bind_rows(WEP.prey_weighted))
  WEP.prey_weighted.df$prey.WT = with(WEP.prey_weighted.df, (prey_FL.wt * prey_Bio.wt))

# Select only sizes of fish encompassed in total predator biomass estimates:
WEP.prey_weighted.df = subset(WEP.prey_weighted.df, PredLength >= 37)

# Categorize prey (pollock vs other):  
WEP.prey_weighted.df$poll = with(WEP.prey_weighted.df, 
              ifelse(PREY_NODC == 8791030701, "pollock", "other"))
WEP.prey_weighted.df$poll = ordered(WEP.prey_weighted.df$poll, levels = c("pollock", "other"))
WEP.prey_weighted.df$Year = as.numeric(as.character(WEP.prey_weighted.df$Year.x))

WEP.coords = st_coordinates(WEP.prey_weighted.df$geometry)
  WEP.prey_weighted.assess = cbind(WEP.prey_weighted.df, WEP.coords) %>%
    rename(Starting.Longitude..dd. = X, Starting.Latitude..dd. = Y) %>%
    filter(Starting.Longitude..dd. < -140)
  
WEP.prey_weighted.assess = unique(WEP.prey_weighted.assess)  
  save(WEP.prey_weighted.assess, 
    file = "4_ProportionsPollockConsumed/propWEP_WEP_37.rda")
# load("4_ProportionsPollockConsumed/propWEP_WEP_37.rda")