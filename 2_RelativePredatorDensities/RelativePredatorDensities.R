# Code by: Cheryl Barnes (cheryl.barnes@noaa.gov)
  # Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# Citation: 
# Dorn MW and CL Barnes. In revision. Time-varying predation as a modifier of constant natural mortality for Gulf of Alaska walleye pollock. Fish Res. VSI Natural Mortality. Forthcoming.

# This script includes code necessary to estimate relative densities of major pollock predators in the Gulf of Alaska; specifically, Arrowtooth Flounder, Pacific Cod, Pacific Halibut, Sablefish, and Walleye Pollock. Methods are described in Barnes et al. (2020). Relative predator densities (kg per ha) were multiplied by total predator biomass, mean annual rations, and proportions of pollock consumed to estimate year- and age-specific predation on pollock. Time-varying predation mortality was then tested as a modifier of assumed constant natural mortality in the main stock assessment model.

# Standardized survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]) and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for details about survey design and data collection. We used delta generalized additive mixed models (GAMMs) to estimate probabilities of occurrence and biomass across a uniform grid spanning the study area. Predictions were then normalized to estimate relative predator density in each unique combination of survey year and grid cell.

# References:
# Barbeaux S, K Aydin, B Fissel, K Holsman, B Laurel, W Palsson, L Rogers, K Shotwell, Q Yang, and S Zador. 2019. Assessment of the Pacific cod stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report. 140pp.
# Barnes CL, Beaudreau AH, Dorn MW, Holsman KK, and Mueter FJ. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecol Appl. 30(7):e02141.
# Clark WG and SR Hare. 2006. Assessment and management of Pacific halibut: data, methods, and policy. IPHC Scientific Report 83. 
# Dorn M, AL Deary, BE Fissell, DT Jones, NE Lauffenburger, W.A. Palsson, LA Rogers, SK Shotwell, KA Spalinger, and S Zador. 2019. Assessment of the Walleye Pollock stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report. 161pp.
# Hanselman DH, CJ Rodgveller, KH Fenske, SK Shotwell, KB Echave, PW Malecha, and CR Lunsford. 2019. Assessment of the Sablefish stock in Alaska. North Pacific Fishery Management Council Bering Sea, Aleutian Islands, and Gulf of Alaska SAFE Report. 263pp.
# Spies I, K Aydin, JN Ianelli, and W Palsson. 2019. Assessment of the Arrowtooth Flounder stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report. 92pp.
# Stewart I and A Hicks. 2020. Assessment of the Pacific halibut (Hippoglossus stenolepis) stock at the end of 2019. International Pacific Halibut Commission Report. 32pp.
# von Szalay PG and NW Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. National Oceanic and Atmospheric Administration. Technical Memorandum NMFS-AFSC-325. 

setwd("~/Documents/AFSC/GOApollock_ESP/")
require(tidyverse)
require(tidyr)
require(dplyr)
require(ggplot2)
require(lme4)
require(MuMIn)
   options(na.action = "na.fail") 
require(sp)
require(maps)
require(sf)
require(mapdata)
require(PBSmapping)
require(mgcv)
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
require(ggmap)
require(FSA)
require(measurements)
require(readxl)
require(rnaturalearth)
   
#########################################################
### INITIAL DATA PREPARATION ###

##### ATF, PC, and WEP:
# Prepare and format AFSC bottom trawl survey data (Arrowtooth Flounder, Pacific Cod, and Walleye Pollock only); provided by Martin Dorn via AKFIN:
# These data include all survey tows conducted between 1984 and 2017.
trawl = read.csv("1_Data/race_cpue_by_haul.csv") 

# Manually assign statistical areas defined by the International North Pacific Fisheries Commission based on survey strata. Note: The second number from the right corresponds with individual statistical areas (e.g., Stratum entry '120' = INPFC StatArea '620'; Stratum entry '251' = StatArea '650').
strata = c(unique(trawl$Stratum))
stat.area = c("650", "650", "650", "610", "610", "610", "610", "610", "620", "620", "620", "630", "630", "630", "630", "640", "630", "630", "630", "630", "630", "610", "610", "610", "610", "620", "630", "620", 
              "630", "640", "640", "640", "640", "640", "650", "650", "650", "650", "630", "640", "620", "630", "620", "650", "620", "620", "610", "630", "640", "630", "620", "640", "640", "640", "630", "610", 
              "620", "650", "640") 
trawl$StatArea = stat.area[match(trawl$Stratum, strata)]

# Exclude data from 1984 and 1987 (survey methods were standardized in 1990):
trawl = subset(trawl, Year >= 1990)

# Treat survey Year as a factor:
trawl$Year = as.factor(trawl$Year)

# Create a unique haul identifier by concatenating Vessel.Number, Cruise.Number, and Haul.Number (allows for joining to proportional length data below):
trawl$Haul.Number_Join = with(trawl, paste(Vessel.Number, Cruise.Number, Haul.Number, sep=""))

# Covert long to wide data format:
trawl_wide = trawl %>%
  pivot_wider(id_cols = c(Survey, 
                          Year, 
                          Starting.Latitude..dd., 
                          Starting.Longitude..dd., 
                          Gear.Depth, 
                          Gear.Temperature...C., 
                          StatArea, 
                          Haul.Number_Join), 
              names_from = Species.Code, 
              values_from = Weight.CPUE..kg.km2., 
              values_fn = sum, 
              values_fill = 0)

# Relabel species codes: 10110 = Arrowtooth Flounder (ATF), 21720 = Pacific Cod (PC), 10120 = Pacific Halibut (PH), 20510 = Sablefish (SBL), 21740 = Walleye Pollock (WEP):
trawl_wide = trawl_wide %>%
  rename(ATF = `10110`, PC = `21720`, PH = `10120`, SBL = `20510`, WEP = `21740`)
  save(trawl_wide, 
    file = "2_RelativePredatorDensities/trawl_survey_data_sub.rda")
# load("2_RelativePredatorDensities/trawl_survey_data_sub.rda")

##### SBL:
# Prepare and format AFSC longline survey data (Sablefish only); provided by Dana Hanselman:
LL = read.csv(unz("1_Data/SBL_catch_summary_view_with_nulls.csv.zip", 
                  "catch_summary_view_with_nulls.csv"), skip=6, header=T, na.strings=c("","NA")) 

# Eliminate extraneous survey years, areas, and surveys (for comparability with AFSC bottom trawl survey):
LL.sub = subset(LL, Year >= 1990 & Year < 2020)
LL.sub = subset(LL.sub, subset = NMFS.Mgmt.Area %in% c("610", "620", "630", "640", "650"))
LL.sub = subset(LL.sub, Survey.Country != "Japan")

# Create unique identifier by concatenating year and station: 
LL.sub$Set.Number_Join = with(LL.sub, paste(Year, Vessel.Number, Station.Number, sep="_"))

# Calculate mean depth per set (m):
LL.sub$Starting.Depth..m. = as.numeric(as.character(LL.sub$Starting.Depth..m.))
LL.sub$Ending.Depth..m. = as.numeric(as.character(LL.sub$Ending.Depth..m.))
LL.sub$Depth_m = (LL.sub$Starting.Depth..m. + LL.sub$Ending.Depth..m.) / 2

# Set NAs to 0:
LL.sub$SumFreq[is.na(LL.sub$SumFreq)] = 0

# Sum numbers of fish and calculate mean depth (m) by Unique ID:
LL_enviro = LL.sub %>%
  group_by(Set.Number_Join) %>%
  mutate(sumFish = sum(SumFreq)) %>%
  mutate(Starting.Latitude..dd. = mean(Start.Latitude..DD.)) %>%
  mutate(Starting.Longitude..dd. = mean(Start.Longitude..DD.)) %>%
  mutate(Gear.Depth = mean(Depth_m))

# Remove unnecessary columns and rename:
LL_enviro = LL_enviro[,c("Year", "Station.Number", "NMFS.Mgmt.Area", "Starting.Latitude..dd.", "Starting.Longitude..dd.", "sumFish", "Set.Number_Join", "Gear.Depth")]
LL_enviro = LL_enviro %>%
  rename(Station = Station.Number, StatArea = NMFS.Mgmt.Area)

##### PH:
# Prepare and format IPHC setline survey data (Pacific Halibut only); available at https://www.iphc.int/datatest/fiss-pacific-halibut-data (need to separately download Excel files due to connection limits and formatting issues):
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

# Select only the regulatory areas encompassed by the Gulf of Alaska:
setline$Region = with(setline,
   ifelse(RegArea == "2C" | RegArea == "3A" | RegArea == "3B" | RegArea == "4A", "GOA", "Not_GOA"))
setline = as.data.frame(setline)
  setline = subset(setline, Region == "GOA")
setline$Year = as.factor(setline$Year)

# Create unique identifier by concatenating year and station: 
setline$Set.Number_Join = with(setline, paste(Year, `Vessel code`, Setno, sep="_"))

# Convert from lb to kg and remove commas from data frame:
setline$O32_kg = setline$`O32 Pacific halibut weight` * 0.453592

# Calculate mean depths (fm) for each set and convert to km:
setline$Gear.Depth = setline$`AvgDepth (fm)` * 1.8288

# Ensure that the same tows are included in each model by removing those with incomplete environmental data (i.e., hauls with missing depths):
setline = subset(setline, !is.na(Gear.Depth))
  setline = setline %>%
    rename(Starting.Longitude..dd. = `MidLon fished`, 
           Starting.Latitude..dd. = `MidLat fished`,
           IPHC.Reg.Area = `IPHC Reg Area`)

SL_enviro = setline[,c("Year", "Station", "IPHC.Reg.Area", "Starting.Latitude..dd.", "Starting.Longitude..dd.", "O32_kg", "Set.Number_Join", "Gear.Depth")]

#########################################################
### INCREASE COMPARABILITY WITH TOTAL BIOMASS ESTIMATES 

##### ATF, PC, and WEP:
# Adjust haul-specific CPUE (kg per ha) to exclude fish not encompassed within assessment-based estimates of total biomass (ATF < 1 yr or < 19 cm; PC < 0 cm or 0 yr; WEP < 37 cm or 3 yr). Note: 100 to 200 fish were subsampled for length measurements per species and haul. 

# Read in and format length data from AFSC bottom trawl survey:
lengths = read.csv(unz("1_Data/race_length_by_haul.csv.zip", "race_length_by_haul.csv"), header=T, skip=7, na.strings=c("","NA"))

# Select and rename columns:
lengths = lengths[,c("Survey","Year","Cruise.Join.ID","Haul.Join.ID","Catch.Join.ID","Cruise.Number", "Vessel.Number","Haul.Number","Starting.Latitude..dd.","Starting.Longitude..dd.","Stratum.INPFC.Area","Gear.Depth","Species.Code","Common.Name","Scientific.Name","Length..mm.","Frequency")]
colnames(lengths) = c("Survey", "Year", "Cruise.Join", "Haul.Join", "Catch.Join", "Cruise.Number", "Vessel.Number", "Haul.Number", "Starting.Latitude..dd.", "Starting.Longitude..dd.", "INPFC.Area", "Gear.Depth", "Species.Code", "Common.Name", "Sci.Name", "Length..mm.", "Frequency")

# Create unique haul identifier by concatenating Vessel.Number, Cruise.Number, and Haul.Number (allows for joining to survey data, prepared above):
lengths$Haul.Number_Join = with(lengths, paste(Vessel.Number, Cruise.Number, Haul.Number, sep=""))

# Exclude data prior to 1990, when survey methods were standardized:
lengths = subset(lengths, Year >= 1990)
  lengths$Year = as.factor(lengths$Year)

# Convert fork length measurements (mm to cm):
lengths$Length..mm. = as.numeric(as.character(lengths$Length..mm.))
lengths$PredLength = lengths$Length..mm. / 10

# Create 10-cm fork length bins (ATF, PC, and WEP only):

# ATF:
ATF = subset(lengths, Common.Name == "arrowtooth flounder")
ATF$FLBin = cut(ATF$PredLength, breaks = c(0,18,28,38,48,58,68,78,88,102))
levels(ATF$FLBin) = c("<=18", "19-28", "29-38", "39-48", "49-58", "59-68", "69-78", "79-88", ">88")

# PC:
PC = subset(lengths, Common.Name == "Pacific cod")
PC$FLBin = cut(PC$PredLength, breaks = c(0,10,20,30,40,50,60,70,80,90,100,110))
levels(PC$FLBin) = c("<=10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100", "101-110")

# WEP:
WEP = subset(lengths, Common.Name == "walleye pollock")
WEP$FLBin = cut(WEP$PredLength, breaks = c(0,16,26,36,46,56,66,76))
levels(WEP$FLBin) = c("<=16", "17-26", "27-36", "37-46", "47-56", "57-66", "67-76")


##### SBL:
# Adjust set-specific CPUE estimates (kg per station) to exclude fish not encompassed within assessment-based estimates of total biomass (Sablefish < 2 yr or < 45 cm). 

# Read in and format length data from AFSC longline survey:
SBL.lengths = read.csv("1_Data/SBL_length_summary_view.csv", header=T, skip=6) 

# Select and rename columns:
SBL.lengths = SBL.lengths[,c("Country", "Year", "Vessel.Number", "Station.Number", "NPFMC.Sablefish.Mgmt.Area", "NMFS.Area.Code", "Start.Latitude..DD.", "Start.Longitude..DD.", "Species.Code", "Common.Name", "Sex", "Length", "Frequencey")]
SBL.lengths = SBL.lengths %>%
  rename(Starting.Latitude..dd. = Start.Latitude..DD., Starting.Longitude..dd = Start.Longitude..DD., Length.cm = Length, Frequency = Frequencey)

# Create unique identifier by concatenating year and station:
SBL.lengths$Set.Number_Join = paste(SBL.lengths$Year, SBL.lengths$Vessel.Number, SBL.lengths$Station, sep="_")

# Eliminate extraneous survey years, areas, and surveys (for comparability with AFSC bottom trawl survey):
SBL.lengths = subset(SBL.lengths, Year >= 1990 & Year < 2020)
  SBL.lengths$Year = as.factor(SBL.lengths$Year)
SBL.lengths = subset(SBL.lengths, NMFS.Area.Code %in% 
                       c("610", "620", "630", "640", "650")) 
SBL.lengths = subset(SBL.lengths, Country != "Japan")

# Create fork length bins (cm):
SBL.lengths$PredLength = as.numeric(as.character(SBL.lengths$Length.cm))
SBL.lengths$FLBin = cut(SBL.lengths$PredLength, 
                        breaks = c(0,14,24,34,44,54,64,74,84,94,120))
levels(SBL.lengths$FLBin) = c("<=14", "15-24", "25-34", "35-44", "45-54", "55-64", "65-74", "75-84", "85-94", ">=95")


##### PH:
# Weights already divvied up by size class.

#######################################################
# Calculate proportions of fish (by weight) sampled in each length bin and haul:

### ATF:
# Predict weight (g) from length (cm) (parameter estimates from Spies et al. 2017):
ATF$WT = 0.004312 * ATF$PredLength ^ 3.186

# Convert to kg and calculate total weight per haul:
ATF$WT_kg = ATF$WT / 1000
ATF$predWT = as.numeric(as.character(ATF$Frequency)) * ATF$WT_kg

# Calculate total weight (kg; all size classes) of ATF, by haul:
ATF_all = ATF %>%
  group_by(Haul.Number_Join) %>%
  summarise(Total_WT = sum(predWT))

# Limit data to sizes of interest (ATF >= 19 cm):
ATF_red = subset(ATF, PredLength >= 19)

# Calculate total weight (kg; >= 19 cm) of ATF, by haulr:
ATF_19 = ATF_red %>%
  group_by(Haul.Number_Join) %>%
  summarise(Recruit_WT = sum(predWT))

# Calculate proportion of haul measuring >= 19 cm:
ATF.prop = merge(ATF_all, ATF_19, all=T); ATF.prop[is.na(ATF.prop)] = 0
ATF.prop$ATF.prop = ATF.prop$Recruit_WT / ATF.prop$Total_WT

# Join bottom trawl survey and proportional length data:
trawl_ATF.prop = trawl_wide %>% left_join(ATF.prop)

# Assign mean proportion of >= 19 cm fish to haulswithout length data:
trawl_ATF.prop$ATF.prop[is.na(trawl_ATF.prop$ATF.prop)] = mean(trawl_ATF.prop$ATF.prop, na.rm=T)

# Adjust CPUE (kg per ha) based on proportion of catch >= 19 cm:
trawl_ATF.prop$adjCPUE = trawl_ATF.prop$ATF * trawl_ATF.prop$ATF.prop

# Ensure that the same tows are included in each model by removing those with incomplete environmental data (i.e., hauls with missing depths or bottom temperatures):
trawl_ATF.comp = subset(trawl_ATF.prop, !is.na(Gear.Depth))
trawl_ATF.comp = subset(trawl_ATF.comp, !is.na(Gear.Temperature...C.))


### PC:
# Predict weight (g) from length (cm) (parameter estimates from Barbeaux et al. 2017):
PC$WT = 0.000005631 * PC$PredLength ^ 3.1306

# Convert to kg and calculate total weight per haul:
PC$WT_kg = PC$WT / 1000
PC$predWT = as.numeric(as.character(PC$Frequency)) * PC$WT_kg

# Calculate total weight (kg; all size classes) of PC, by haul:
PC_all = PC %>%
  group_by(Haul.Number_Join) %>%
  summarise(Total_WT = sum(predWT))

# Limit data to sizes of interest (PC >= 0 cm):
PC_red = subset(PC, PredLength >= 0)

# Calculate total weight (kg; >= 0 cm) of PC, by haul:
PC_19 = PC_red %>%
  group_by(Haul.Number_Join) %>%
  summarise(Recruit_WT = sum(predWT))

# Calculate proportion of haul measuring >= 0 cm:
PC.prop = merge(PC_all, PC_19, all=T); PC.prop[is.na(PC.prop)] = 0
PC.prop$PC.prop = PC.prop$Recruit_WT / PC.prop$Total_WT

# Join bottom trawl survey and proportional length data:
trawl_PC.prop = trawl_wide %>% left_join(PC.prop)

# Assign mean proportion of >= 0 cm fish to hauls without length data:
trawl_PC.prop$PC.prop[is.na(trawl_PC.prop$PC.prop)] = mean(trawl_PC.prop$PC.prop, na.rm=T)

# Adjust CPUE (kg per ha) based on proportion of catch >= 0 cm:
trawl_PC.prop$adjCPUE = trawl_PC.prop$PC * trawl_PC.prop$PC.prop

# Ensure that the same tows are included in each model by removing those with incomplete environmental data (i.e., hauls with missing depths or bottom temperatures):
trawl_PC.comp = subset(trawl_PC.prop, !is.na(Gear.Depth))
trawl_PC.comp = subset(trawl_PC.comp, !is.na(Gear.Temperature...C.))


### WEP:
# Predict weight (g) from length (cm)  using an allometric model with bias-correction (Brodziak 2012):
LW = read.csv(unz("1_Data/race_specimen.csv.zip", 
                  "race_specimen.csv"), skip=7, header=T, na.strings=c("","NA"))  # 1990 to 2019
WEP.LW = subset(LW, Common.Name == "walleye pollock")
  WEP.LW$logW = log(WEP.LW$Weight..gm.)
  WEP.LW$logL = log(WEP.LW$Length..mm.)
WEP.LW = subset(WEP.LW, logW !="NA")

L_W = lm(logW ~ logL, data=WEP.LW)
  summary(L_W); exp(coef(L_W)[1])

# Calculate bias correction factor for predicting on original scale:
syx = summary(L_W)$sigma
cf = exp((syx^2)/2) 

WEP$WT_log = predict(L_W, data.frame(logL=log(WEP$Length..mm.)), interval="c")
WEP$WT_log = as.numeric(WEP$WT_log[1])
WEP$WT_unbi = cf *(exp(WEP$WT_log)) 

# Convert to kg and calculate total weight per haul:
WEP$WT_kg = WEP$WT_unbi / 1000
WEP$predWT = as.numeric(as.character(WEP$Frequency)) * WEP$WT_kg

# Calculate total weight (kg; all size classes) of WEP, by haul:
WEP_all = WEP %>%
  group_by(Haul.Number_Join) %>%
  summarise(Total_WT = sum(predWT))

# Limit data to sizes of interest (WEP >= 37 cm):
WEP_red = subset(WEP, PredLength >= 37)

# Calculate total weight (kg; >= 37 cm) of WEP, by haul:
WEP_19 = WEP_red %>%
  group_by(Haul.Number_Join) %>%
  summarise(Recruit_WT = sum(predWT))

# Calculate proportion of haul measuring >= 37 cm:
WEP.prop = merge(WEP_all, WEP_19, all=T); WEP.prop[is.na(WEP.prop)] = 0
WEP.prop$WEP.prop = WEP.prop$Recruit_WT / WEP.prop$Total_WT

# Join bottom trawl survey and proportional length data:
trawl_WEP.prop = trawl_wide %>% left_join(WEP.prop)

# Assign mean proportion of >= 37 cm fish to hauls without length data:
trawl_WEP.prop$WEP.prop[is.na(trawl_WEP.prop$WEP.prop)] = mean(trawl_WEP.prop$WEP.prop, na.rm=T)

# Adjust CPUE (kg per ha) based on proportion of catch >= 37 cm:
trawl_WEP.prop$adjCPUE = trawl_WEP.prop$WEP * trawl_WEP.prop$WEP.prop

# Ensure that the same tows are included in each model by removing those with incomplete environmental data (i.e., hauls with missing depths or bottom temperatures):
trawl_WEP.comp = subset(trawl_WEP.prop, !is.na(Gear.Depth))
trawl_WEP.comp = subset(trawl_WEP.comp, !is.na(Gear.Temperature...C.))


### SBL:
# Estimate weight (kg) from length (cm) (parameter estimates from Hanselman et al. 2007 - Table 5; 1 = male, 2 = female):
SBL.lengths$WT_est = with(SBL.lengths, 
  ifelse(Sex == "1", (0.0000124 * PredLength ^ 2.960),
  ifelse(Sex == "2", (0.0000101 * PredLength ^ 3.015), NA)))

# Calculate mean weight at length for fish with no sex identified:
SBL.lengths = SBL.lengths %>%
  group_by(PredLength) %>%
  mutate(WT_calc = mean(na.omit(WT_est)))

SBL.lengths$WT = SBL.lengths$WT_est
SBL.lengths$WT = with(SBL.lengths,
  ifelse(is.na(WT), WT_calc, WT_est))
SBL.lengths = na.omit(SBL.lengths)

# Calculate total weight (kg; all size classes) by year and station:
SBL.lengths$Frequency = as.numeric(as.character(SBL.lengths$Frequency))
SBL.lengths$predWT = with(SBL.lengths, Frequency * WT)

LL_all = SBL.lengths %>%
  group_by(Set.Number_Join) %>%
  summarise(Total_WT = sum(predWT))

# Limit data to sizes of interest (SBL >= 45 cm):
SBL.lengths_red = subset(SBL.lengths, PredLength >= 45)

# Calculate total weight (kg; >= 45 cm) of SBL, by year and station:
LL_45 = SBL.lengths_red %>%
  group_by(Set.Number_Join) %>%
  summarise(Recruit_WT = sum(WT))

# Calculate proportion of set measuring >= 45 cm:
SBL.prop = as.data.frame(LL_all %>% 
                           full_join(LL_45)); SBL.prop[is.na(SBL.prop)] = 0
SBL.prop$SBL.prop = SBL.prop$Recruit_WT / SBL.prop$Total_WT

# Join longline survey and proportional length data:
LL_prop = LL_enviro %>% left_join(SBL.prop)

# Assign weights with NAs to zero when no fish were caught:
LL_prop$Total_WT[is.na(LL_prop$Total_WT) & LL_prop$sumFish==0] = 0
LL_prop$Recruit_WT[is.na(LL_prop$Recruit_WT) & LL_prop$sumFish==0] = 0
LL_prop$SBL.prop[is.na(LL_prop$SBL.prop) & LL_prop$sumFish==0] = 0

# Remove sets with weights but no fish and stations without lengthed fish:
LL_prop = subset(LL_prop, subset =!(sumFish == 0 & Total_WT > 0))
  LL_prop = na.omit(LL_prop)

# Adjust CPUE (kg per station) based on proportions of catch >= 45 cm:
LL_prop$adjCPUE = LL_prop$Total_WT * LL_prop$SBL.prop

# Remove extreme outlier (n = 1, depth > 1000 m):
LL_prop = subset(LL_prop, Gear.Depth < 1000)
LL_prop$Year = as.factor(LL_prop$Year)

# Ensure that the same sets are included in each model by removing those with incomplete environmental data (i.e., sets with missing depths; temperature excluded from analyses):
ll_SBL.comp = subset(LL_prop, !is.na(Gear.Depth))


### PH:
sl_PH.comp = setline

#########################################################
### MODEL FITTING ###

### ATF, presence-absence 
# Label each haul as present or absent for Arrowtooth Flounder:
trawl_ATF.comp$ATFpa = as.numeric(trawl_ATF.comp$adjCPUE > 0)
  length(trawl_ATF.comp$ATFpa) # total number of hauls conducted
sum(na.omit(trawl_ATF.comp$ATFpa)) # total number of hauls that captured arrowtooth

# Run the full model, without spatial autocorrelation:
ATF.pa.gam_full = gam(ATFpa ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = trawl_ATF.comp, family = binomial(link=logit), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
ATF.pa.gam_select = dredge(ATF.pa.gam_full, beta=F, evaluate=T, rank="AIC", trace=F)
  print(ATF.pa.gam_select, abbrev.names=F, warnings=T) 
    summary(ATF.pa.gam_full) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey Year:
ATF.pa.gamm_full = gamm(ATFpa ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), correlation = corGaus(form = ~ Starting.Longitude..dd. + Starting.Latitude..dd. | Year), data = trawl_ATF.comp, family = binomial(link=logit))
  summary(ATF.pa.gamm_full)
  summary(ATF.pa.gamm_full$gam)
  summary(ATF.pa.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
ATF.pa.gamm_full_noSA = gamm(ATFpa ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = trawl_ATF.comp, family = binomial(link=logit))

AIC(ATF.pa.gamm_full, ATF.pa.gamm_full_noSA) 
    # GAMM with spatial autocorrelation = best-fit model
  save(ATF.pa.gamm_full, 
    file = "2_RelativePredatorDensities/ATF_pa_gamm_full.rda")
# load("2_RelativePredatorDensities/ATF_pa_gamm_full.rda")
    
### ATF, CPUE (where present) 
# Subset hauls to include only those that sampled Arrowtooth Flounder:
ATF = subset(trawl_ATF.comp, adjCPUE > 0)
ATF$logATF = log(ATF$adjCPUE) # log-transform CPUE data

# Run the full model, without spatial autocorrelation:
ATF.cpue.gam_full = gam(logATF ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = ATF, family = gaussian(link=identity), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
ATF.cpue.gam_select = dredge(ATF.cpue.gam_full, beta=F, evaluate=T, rank="AIC", trace=F)
  print(ATF.cpue.gam_select, abbrev.names=F, warnings=T) 
  summary(ATF.cpue.gam_full) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey Year:
ATF.cpue.gamm_full = gamm(logATF ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), correlation = corGaus(form = ~ Starting.Longitude..dd. + Starting.Latitude..dd. | Year), data = ATF, family = gaussian(link=identity))
  summary(ATF.cpue.gamm_full)
  summary(ATF.cpue.gamm_full$gam)
  summary(ATF.cpue.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
ATF.cpue.gamm_full_noSA = gamm(logATF ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = ATF, family = gaussian(link=identity))

AIC(ATF.cpue.gamm_full, ATF.cpue.gamm_full_noSA)
    # GAMM with spatial autocorrelation = best-fit model
  save(ATF.cpue.gamm_full, 
  file = "2_RelativePredatorDensities/ATF_cpue_gamm_full.rda")
# load("2_RelativePredatorDensities/ATF_cpue_gamm_full.rda")
      

### PC, presence-absence 
# Label each haul as present or absent for Pacific Cod:
trawl_PC.comp$PCpa = as.numeric(trawl_PC.comp$adjCPUE > 0)
  length(trawl_PC.comp$PCpa) # total number of hauls conducted
sum(na.omit(trawl_PC.comp$PCpa)) # total number of hauls that captured P. cod

# Run the full model, without spatial autocorrelation:
PC.pa.gam_full = gam(PCpa ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = trawl_PC.comp, family = binomial(link=logit), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
PC.pa.gam_select = dredge(PC.pa.gam_full, beta=F, evaluate=T, rank="AIC", trace=F)
  print(PC.pa.gam_select, abbrev.names=F, warnings=T) 
    summary(PC.pa.gam_full) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey Year:
PC.pa.gamm_full = gamm(PCpa ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), correlation = corGaus(form = ~ Starting.Longitude..dd. + Starting.Latitude..dd. | Year), data = trawl_PC.comp, family = binomial(link=logit))
  summary(PC.pa.gamm_full)
  summary(PC.pa.gamm_full$gam)
  summary(PC.pa.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
PC.pa.gamm_full_noSA = gamm(PCpa ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = trawl_PC.comp, family = binomial(link=logit))

AIC(PC.pa.gamm_full, PC.pa.gamm_full_noSA) 
    # GAMM with spatial autocorrelation = best-fit model
  save(PC.pa.gamm_full, 
    file = "2_RelativePredatorDensities/PC_pa_gamm_full.rda")
# load("2_RelativePredatorDensities/PC_pa_gamm_full.rda")

### PC, CPUE (where present) 
# Subset hauls to include only those that sampled Pacific Cod:
PC = subset(trawl_PC.comp, adjCPUE > 0)
PC$logPC = log(PC$adjCPUE) # log-transform CPUE data

# Run the full model, without spatial autocorrelation:
PC.cpue.gam_full = gam(logPC ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = PC, family = gaussian(link=identity), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
PC.cpue.gam_select = dredge(PC.cpue.gam_full, beta=F, evaluate=T, rank="AIC", trace=F)
  print(PC.cpue.gam_select, abbrev.names=F, warnings=T) 
  summary(PC.cpue.gam_full) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey Year:
PC.cpue.gamm_full = gamm(logPC ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), correlation = corGaus(form = ~ Starting.Longitude..dd. + Starting.Latitude..dd. | Year), data = PC, family = gaussian(link=identity))
  summary(PC.cpue.gamm_full)
  summary(PC.cpue.gamm_full$gam)
  summary(PC.cpue.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
PC.cpue.gamm_full_noSA = gamm(logPC ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = PC, family = gaussian(link=identity))

AIC(PC.cpue.gamm_full, PC.cpue.gamm_full_noSA)
    # GAMM with spatial autocorrelation = best-fit model
  save(PC.cpue.gamm_full, 
    file = "2_RelativePredatorDensities/PC_cpue_gamm_full.rda")  
# load("2_RelativePredatorDensities/PC_cpue_gamm_full.rda")
        
        
### WEP, presence-absence 
# Label each haul as present or absent for Walleye Pollock:
trawl_WEP.comp$WEPpa = as.numeric(trawl_WEP.comp$adjCPUE > 0)
  length(trawl_WEP.comp$WEPpa) # total number of hauls conducted
sum(na.omit(trawl_WEP.comp$WEPpa)) # total number of hauls that captured pollock

# Run the full model, without spatial autocorrelation:
WEP.pa.gam_full = gam(WEPpa ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = trawl_WEP.comp, family = binomial(link=logit), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
WEP.pa.gam_select = dredge(WEP.pa.gam_full, beta=F, evaluate=T, rank="AIC", trace=F)
  print(WEP.pa.gam_select, abbrev.names=F, warnings=T) 
    summary(WEP.pa.gam_full) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey Year:
WEP.pa.gamm_full = gamm(WEPpa ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), correlation = corGaus(form = ~ Starting.Longitude..dd. + Starting.Latitude..dd. | Year), data = trawl_WEP.comp, family = binomial(link=logit))
  summary(WEP.pa.gamm_full)
  summary(WEP.pa.gamm_full$gam)
  summary(WEP.pa.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
WEP.pa.gamm_full_noSA = gamm(WEPpa ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = trawl_WEP.comp, family = binomial(link=logit))

AIC(WEP.pa.gamm_full, WEP.pa.gamm_full_noSA)
    # GAMM without spatial autocorrelation = best-fit model
  save(WEP.pa.gam_full, 
    file = "2_RelativePredatorDensities/WEP_pa_gam_full.rda")
# load("2_RelativePredatorDensities/WEP_pa_gam_full.rda")

### WEP, CPUE (where present) 
# Subset hauls to include only those that sampled Walleye Pollock:
WEP = subset(trawl_WEP.comp, adjCPUE > 0)
WEP$logWEP = log(WEP$adjCPUE) # log-transform CPUE data

# Run the full model, without spatial autocorrelation:
WEP.cpue.gam_full = gam(logWEP ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = WEP, family = gaussian(link=identity), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
WEP.cpue.gam_select = dredge(WEP.cpue.gam_full, beta=F, evaluate=T, rank="AIC", trace=F)
  print(WEP.cpue.gam_select, abbrev.names=F, warnings=T) 
  summary(WEP.cpue.gam_full) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey Year:
WEP.cpue.gamm_full = gamm(logWEP ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), correlation = corGaus(form = ~ Starting.Longitude..dd. + Starting.Latitude..dd. | Year), data = WEP, family = gaussian(link=identity))
  summary(WEP.cpue.gamm_full)
  summary(WEP.cpue.gamm_full$gam)
  summary(WEP.cpue.gamm_full$mer)

# Rerun GAM as a GAMM without spatial autocorrelation. Use AIC to select best model (with penalty for additional term):
WEP.cpue.gamm_full_noSA = gamm(logWEP ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth) + s(Gear.Temperature...C., k=4), data = WEP, family = gaussian(link=identity))

AIC(WEP.cpue.gamm_full, WEP.cpue.gamm_full_noSA)
    # GAMM with spatial autocorrelation = best-fit model
  save(WEP.cpue.gamm_full, 
    file = "2_RelativePredatorDensities/WEP_cpue_gamm_full.rda")
# load("2_RelativePredatorDensities/WEP_cpue_gamm_full.rda")


### SBL, presence-absence 
# Label each longline set as present or absent for Sablefish:
ll_SBL.comp$SBLpa = as.numeric(ll_SBL.comp$adjCPUE > 0)
  length(ll_SBL.comp$SBLpa) # total number of hauls conducted
sum(ll_SBL.comp$SBLpa) # total number of hauls that captured SBL
  # Few stations (< 100) did not sample Sablefish.
      # Do not run separate presence-absence model.

### SBL CPUE (where present) 
SBL = ll_SBL.comp
SBL$logSBL = log(SBL$adjCPUE + 0.001) # log-transform CPUE data and add small constant for zeros

# Run the full model, without spatial autocorrelation:
SBL.cpue.gam_full = gam(logSBL ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth), data = SBL, family = gaussian(link=identity), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
SBL.cpue.gam_select = dredge(SBL.cpue.gam_full, beta=F, evaluate=T, rank="AIC", trace=F)
  print(SBL.cpue.gam_select, abbrev.names=F, warnings=T) 
  summary(SBL.cpue.gam_full); sum(SBL.cpue.gam_full$edf) # full model = best-fit model

# Rerun the best-fit model with a spatial autocorrelation term that varies by survey year:
SBL.cpue.gamm_full = gamm(logSBL ~ Year + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth), correlation = corGaus(form = ~ Starting.Longitude..dd. + Starting.Latitude..dd. | Year), data = SBL, family = gaussian(link=identity))
  # GAMM does not converge; GAM without spatial autocorrelation = best-fit model
  save(SBL.cpue.gam_full, 
    file = "2_RelativePredatorDensities/SBL_cpue_gam_full.rda")
# load("2_RelativePredatorDensities/SBL_cpue_gam_full.rda")
  

### PH, presence-absence:
# Label each setline set as present or absent for Pacific Halibut:
sl_PH.comp$PHpa = as.numeric(sl_PH.comp$O32_kg > 0)
  length(sl_PH.comp$PHpa) # total number of hauls conducted
sum(sl_PH.comp$PHpa) # total number of hauls that captured PH
  # Few stations (< 300) did not sample P. Halibut.
      # Do not run separate presence-absence model.

### PH CPUE (where present) 
PH = sl_PH.comp
PH$logPH = log(PH$O32_kg + 0.001) # log-transform CPUE data and add small constant for zeros

# Run the full model, without spatial autocorrelation (different representation of year to allow for predictions outside the range of years used in fitting):
PH$Year = as.integer(as.character(PH$Year))
PH.cpue.gam_full = gam(logPH ~ s(Year, bs="cr", m=1) + s(Starting.Longitude..dd., Starting.Latitude..dd.) + s(Gear.Depth), data = PH, family = gaussian(link=identity), method="GCV.Cp")

# Generate all possible alternative models and select the best-fit based on delta AIC:
PH.cpue.gam_select = dredge(PH.cpue.gam_full, beta=FALSE, evaluate=T, rank="AIC", trace=FALSE)
  print(PH.cpue.gam_select, abbrev.names=FALSE, warnings=T) # Table S2 
  summary(PH.cpue.gam_full) # full model = best-fit model
    
# Rerun the best-fit model with a spatial autocorrelation term that varies by survey year:
  # GAMM for Pacific Halibut does not converge. Use GAM for probabilities of occurrence.
  save(PH.cpue.gam_full, 
    file = "2_RelativePredatorDensities/PH.cpue.gam_full.rda")
# load("2_RelativePredatorDensities/PH.cpue.gam_full.rda")

# Bottom temperature was not included as a covariate in models due to limited data collection on the IPHC setline survey).
#########################################################
### CREATE PREDICTION GRIDS ###

# Convert bottom trawl survey data to a spatial data frame and join to grid:
coordinates(trawl_wide) = ~ Starting.Longitude..dd. + Starting.Latitude..dd.
  trawl_enviro = st_as_sf(trawl_wide)
    st_crs(trawl_enviro) = 4326
  trawl.utm = st_transform(trawl_enviro, 26935)

trawl.mask = st_buffer(trawl.utm, dist=25000)
trawl.grid_utm = st_make_grid(trawl.utm, square=T, cellsize = c(50000,50000), what="centers")

trawl.grid = st_intersection(trawl.mask, trawl.grid_utm) # limits grid extent to survey area
trawl.grid.dd = st_transform(trawl.grid, 4326)

# Isolate lon/lat, calculate mean depth at each location, and eliminate extraneous columns:
trawl.grid.df = as.data.frame(trawl.grid.dd)
  trawl.coords = st_coordinates(trawl.grid.df$geometry)
trawl_grid.dd = cbind(trawl.grid.df, trawl.coords)

trawl_grid.dd = trawl_grid.dd %>%
  group_by(StatArea, X, Y) %>%
  summarise(Gear.Depth_m = mean(Gear.Depth, na.rm=T)) %>%
  rename(Starting.Longitude..dd. = X, 
         Starting.Latitude..dd. = Y, 
         Gear.Depth = Gear.Depth_m) %>%
  na.omit()

# Add year-specific temperature data to grid at nearest location:
trawl_grid.sf = st_as_sf(trawl_grid.dd, 
              coords = c("Starting.Longitude..dd.", "Starting.Latitude..dd."),
              crss = 4326)
trawl_temp = trawl_enviro[,c("Year", 
                             "geometry",
                             "Gear.Temperature...C.")]
trawl.enviro.data = list()
for(i in unique(trawl_temp$Year)) {
  survey.yr = subset(trawl_temp, Year == i)
  trawl.enviro.data[[i]] = st_join(trawl_grid.sf, survey.yr, 
                                   join = st_nearest_feature, left = T) }
trawl.enviro.df = as.data.frame(dplyr::bind_rows(trawl.enviro.data))

# Reformat as data frame:
trawl.coords = st_coordinates(trawl.enviro.df$geometry)
  trawl.input.data = cbind(trawl.enviro.df, trawl.coords)
  trawl.input.data = trawl.input.data %>%
    rename(Starting.Longitude..dd. = X, Starting.Latitude..dd. = Y) 
  
  save(trawl.input.data, 
    file = "2_RelativePredatorDensities/trawl_grid.rda")
# load("2_RelativePredatorDensities/trawl_grid.rda")  
  
  
# Convert AFSC longline survey data to a spatial data frame and join to grid:
coordinates(LL_enviro) = ~ Starting.Longitude..dd. + Starting.Latitude..dd.
  LL.enviro = st_as_sf(LL_enviro)
    st_crs(LL.enviro) = 4326
  LL.utm = st_transform(LL.enviro, 26935)

LL.mask = st_buffer(LL.utm, dist=25000)
LL.grid = st_intersection(LL.mask, trawl.grid_utm) # limits grid extent to survey area
LL.grid.dd = st_transform(LL.grid, 4326)

# Reformat as data frame:
LL.grid.df = as.data.frame(LL.grid.dd)
LL.coords = st_coordinates(LL.grid.df$geometry)
  LL.grid.data = cbind(LL.grid.df, LL.coords)
LL.grid.data = LL.grid.data %>%
    group_by(StatArea, X, Y) %>%
    summarise(Gear.Depth_m = mean(Gear.Depth, na.rm=T)) %>%
    as.data.frame() %>%
    rename(Starting.Longitude..dd. = X, Starting.Latitude..dd. = Y, 
         Gear.Depth = Gear.Depth_m) %>%
  na.omit()
  
# Add survey years (no temperature data for this model):
LL.input.data = LL.grid.data[rep(seq_len(nrow(LL.grid.data)), each = 14),]
LL.input.data$Year = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019)

LL.input.data$Year = as.factor(LL.input.data$Year)
  LL.input.data$StatArea = as.character(LL.input.data$StatArea)
LL.input.data = unique(na.omit(LL.input.data))
  
  save(LL.input.data, 
    file = "2_RelativePredatorDensities/LL_grid.rda")
# load("2_RelativePredatorDensities/LL_grid.rda")

  
# Convert IPHC setline survey data to a spatial data frame and join to grid:
coordinates(SL_enviro) = ~ Starting.Longitude..dd. + Starting.Latitude..dd.
  SL.enviro = st_as_sf(SL_enviro)
    st_crs(SL.enviro) = 4326
  SL.utm = st_transform(SL.enviro, 26935)

SL.mask = st_buffer(SL.utm, dist=25000)
SL.grid = st_intersection(SL.mask, trawl.grid_utm) # limits grid extent to survey area
SL.grid.dd = st_transform(SL.grid, 4326)

# Reformat as data frame:
SL.grid.df = as.data.frame(SL.grid.dd)
SL.coords = st_coordinates(SL.grid.df$geometry)
  SL.grid.data = cbind(SL.grid.df, SL.coords)
SL.grid.data = SL.grid.data %>%
    group_by(IPHC.Reg.Area, X, Y) %>%
    summarise(Gear.Depth_m = mean(Gear.Depth, na.rm=T)) %>%
    as.data.frame() %>%
    rename(Starting.Longitude..dd. = X, Starting.Latitude..dd. = Y, 
         Gear.Depth = Gear.Depth_m) %>%
  na.omit()

# Add survey years (no temperature data for this model):
SL.input.data = SL.grid.data[rep(seq_len(nrow(SL.grid.data)), each = 14),]
SL.input.data$Year = c(1990,1993,1996,1999,2001,2003,2005,2007,2009,2011,2013,2015,2017,2019)

SL.input.data$Year = as.factor(SL.input.data$Year)
  SL.input.data$IPHC.Reg.Area = as.character(SL.input.data$IPHC.Reg.Area)
SL.input.data = unique(na.omit(SL.input.data))

  save(SL.input.data, 
    file = "2_RelativePredatorDensities/SL_grid.rda")
# load("2_RelativePredatorDensities/SL_grid.rda")

#########################################################
### MODEL PREDICTIONS ###
# Remove locations without temperature data (all spp. except SBL):
trawl.input.data_comp = na.omit(trawl.input.data)
    
# ATF, presence-absence:
ATFpa_predict = predict.gam(ATF.pa.gamm_full$gam, newdata=trawl.input.data_comp, type="response", se.fit=T)
ATFpa_pred = cbind(trawl.input.data_comp, ATFpa_predict)

ATFpa_pred.CI = within(ATFpa_pred, {
  lower = fit - 1.96 * se.fit
  upper = fit + 1.96 * se.fit })
 
ATFpa_pred.CI = ATFpa_pred.CI %>%
  rename(ATFpa_fit = fit, ATFpa_se.fit = se.fit, ATFpa_hiCI = upper, ATFpa_loCI = lower)

# ATF, CPUE:
ATF.cpue_predict = predict.gam(ATF.cpue.gamm_full$gam, newdata=trawl.input.data_comp, type="response", se.fit=T)
ATF.cpue_pred = cbind(trawl.input.data_comp, ATF.cpue_predict)

ATF.cpue_pred.CI = within(ATF.cpue_pred, {
  lower = fit - 1.96 * se.fit
  upper = fit + 1.96 * se.fit })

ATF.cpue_pred.CI = ATF.cpue_pred.CI %>%
  rename(ATF.cpue_fit = fit, ATF.cpue_se.fit = se.fit, ATF.cpue_hiCI = upper, ATF.cpue_loCI = lower)

# ATF (occurrence-informed CPUE):
ATF.predictions = merge(ATFpa_pred.CI, ATF.cpue_pred.CI)
  ATF.predictions$ATF.predAbun = with(ATF.predictions, ATFpa_fit * ATF.cpue_fit)

# Normalize so that all grid cells in a given survey year sum to one:
normAbun_ATF = ATF.predictions %>%
  group_by(Year) %>%
  mutate(ATF.normAbun = ATF.predAbun / sum(ATF.predAbun))
normAbun_ATF = as.data.frame(normAbun_ATF)
    save(normAbun_ATF, 
      file = "2_RelativePredatorDensities/normATFabun.rda")
  # load("2_RelativePredatorDensities/normATFabun.rda")


# PC, presence-absence:
PCpa_predict = predict.gam(PC.pa.gamm_full$gam, newdata=trawl.input.data_comp, type="response", se.fit=T)
PCpa_pred = cbind(trawl.input.data_comp, PCpa_predict)

PCpa_pred.CI = within(PCpa_pred, {
  lower = fit - 1.96 * se.fit
  upper = fit + 1.96 * se.fit })
 
PCpa_pred.CI = PCpa_pred.CI %>%
  rename(PCpa_fit = fit, PCpa_se.fit = se.fit, PCpa_hiCI = upper, PCpa_loCI = lower)

# PC, CPUE:
PC.cpue_predict = predict.gam(PC.cpue.gamm_full$gam, newdata=trawl.input.data, type="response", se.fit=T)
PC.cpue_pred = cbind(trawl.input.data, PC.cpue_predict)

PC.cpue_pred.CI = within(PC.cpue_pred, {
  lower = fit - 1.96 * se.fit
  upper = fit + 1.96 * se.fit })

PC.cpue_pred.CI = PC.cpue_pred.CI %>%
  rename(PC.cpue_fit = fit, PC.cpue_se.fit = se.fit, PC.cpue_hiCI = upper, PC.cpue_loCI = lower)

# PC (occurrence-informed CPUE):
PC.predictions = merge(PCpa_pred.CI, PC.cpue_pred.CI)
  PC.predictions$PC.predAbun = with(PC.predictions, PCpa_fit * PC.cpue_fit)

# Normalize so that all grid cells in a given survey year sum to one:
normAbun_PC = PC.predictions %>%
  group_by(Year) %>%
  mutate(PC.normAbun = PC.predAbun / sum(PC.predAbun))
normAbun_PC = as.data.frame(normAbun_PC)
    save(normAbun_PC, 
      file = "2_RelativePredatorDensities/normPCabun.rda")
  # load("2_RelativePredatorDensities/normPCabun.rda")


# WEP, presence-absence:
WEPpa_predict = predict.gam(WEP.pa.gam_full, newdata=trawl.input.data_comp, type="response", se.fit=T)
WEPpa_pred = cbind(trawl.input.data_comp, WEPpa_predict)

WEPpa_pred.CI = within(WEPpa_pred, {
  lower = fit - 1.96 * se.fit
  upper = fit + 1.96 * se.fit })
 
WEPpa_pred.CI = WEPpa_pred.CI %>%
  rename(WEPpa_fit = fit, WEPpa_se.fit = se.fit, WEPpa_hiCI = upper, WEPpa_loCI = lower)

# WEP, CPUE:
WEP.cpue_predict = predict.gam(WEP.cpue.gamm_full$gam, newdata=trawl.input.data, type="response", se.fit=T)
WEP.cpue_pred = cbind(trawl.input.data, WEP.cpue_predict)

WEP.cpue_pred.CI = within(WEP.cpue_pred, {
  lower = fit - 1.96 * se.fit
  upper = fit + 1.96 * se.fit })

WEP.cpue_pred.CI = WEP.cpue_pred.CI %>%
  rename(WEP.cpue_fit = fit, WEP.cpue_se.fit = se.fit, WEP.cpue_hiCI = upper, WEP.cpue_loCI = lower)

# WEP (occurrence-informed CPUE):
WEP.predictions = merge(WEPpa_pred.CI, WEP.cpue_pred.CI)
  WEP.predictions$WEP.predAbun = with(WEP.predictions, WEPpa_fit * WEP.cpue_fit)

# Normalize so that all grid cells in a given survey year sum to one:
normAbun_WEP = WEP.predictions %>%
  group_by(Year) %>%
  mutate(WEP.normAbun = WEP.predAbun / sum(WEP.predAbun))
normAbun_WEP = as.data.frame(normAbun_WEP)
    save(normAbun_WEP, 
      file = "2_RelativePredatorDensities/normWEPabun.rda")
  # load("2_RelativePredatorDensities/normWEPabun.rda")


### SBL, CPUE (presence-absence not modeled):
SBL.cpue_predict = predict.gam(SBL.cpue.gam_full, newdata=LL.input.data, type="response", se.fit=T)
SBL.cpue_pred = cbind(LL.input.data, SBL.cpue_predict)

# Calculate 95% confidence intervals:
SBL.cpue_pred.CI = within(SBL.cpue_pred, {
  lower = fit - 1.96 * se.fit
  upper = fit + 1.96 * se.fit
})

# Rename columns to match initial input data:
SBL.cpue_pred.CI = SBL.cpue_pred.CI %>%
  rename(SBL.cpue_fit = fit, SBL.cpue_se.fit = se.fit, SBL.cpue_hiCI = upper, SBL.cpue_loCI = lower)

# Rename to match all species (with hurdle models):
SBL.predictions = SBL.cpue_pred.CI
  SBL.predictions$SBL.predAbun = SBL.predictions$SBL.cpue_fit
  
# Normalize so that all grid cells in a given survey year sum to one:
normAbun_SBL = SBL.predictions %>%
  group_by(Year) %>%
  mutate(SBL.normAbun = SBL.predAbun / sum(SBL.predAbun))
normAbun_SBL = as.data.frame(normAbun_SBL)
    save(normAbun_SBL, 
      file = "2_RelativePredatorDensities/normSBLabun.rda")
  # load("2_RelativePredatorDensities/normSBLabun.rda")


### PH, CPUE (presence-absence not modeled):
PH.cpue_predict = predict.gam(PH.cpue.gam_full, newdata=SL.input.data, type="response", se.fit=T)
PH.cpue_pred = cbind(SL.input.data, PH.cpue_predict)

# Calculate 95% confidence intervals:
PH.cpue_pred.CI = within(PH.cpue_pred, {
  lower = fit - 1.96 * se.fit
  upper = fit + 1.96 * se.fit
})

# Rename columns to match initial input data:
PH.cpue_pred.CI = PH.cpue_pred.CI %>%
  rename(PH.cpue_fit = fit, PH.cpue_se.fit = se.fit, PH.cpue_hiCI = upper, PH.cpue_loCI = lower)

# Rename to match all species (with hurdle models):
PH.predictions = PH.cpue_pred.CI
  PH.predictions$PH.predAbun = PH.predictions$PH.cpue_fit
  
# Normalize so that all grid cells in a given survey year sum to one:
normAbun_PH = PH.predictions %>%
  group_by(Year) %>%
  mutate(PH.normAbun = PH.predAbun / sum(PH.predAbun))
normAbun_PH = as.data.frame(normAbun_PH)
    save(normAbun_PH, 
      file = "2_RelativePredatorDensities/normPHabun.rda")
  # load("2_RelativePredatorDensities/normPHabun.rda")