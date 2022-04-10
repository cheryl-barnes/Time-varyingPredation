# Time-varying predation as a modifier of constant natural mortality for Gulf of Alaska walleye pollock <br><br>

## Citation
Dorn MW and CL Barnes. In revision. Time-varying predation as a modifier of constant natural mortality for Gulf of Alaska walleye pollock. Fish Res. VSI Natural Mortality. Forthcoming. <br><br>

<b> Martin W. Dorn </b><br>
Alaska Fisheries Science Center, National Marine Fisheries Service (NOAA) <br>
martin.dorn@noaa.gov <br>

<b> Cheryl L. Barnes </b><br>
School of Aquatic and Fishery Sciences, University of Washington <br>
cheryl.barnes@noaa.gov <br>

## Overview
This repository details the methods used to estimate age- and year-specific predation on Walleye Pollock (<i>Gadus chalcogrammus</i>) in the Gulf of Alaska (MT per year; 1990 to 2019). Pollock predators included: Arrowtooth Flounder (<i>Atheresthes stomias</i>), Pacific Cod (<i>Gadus macrocephalus</i>), Pacific Halibut (<i>Hippoglossus stenolepis</i>), Sablefish (<i>Anoplopoma fimbria</i>), and Walleye Pollock conspecifics. We used predation indices to test the effects of including time-varying natural mortality within the stock assessment model for Gulf of Alaska pollock.

## File Structure
Input data (survey and food habits) can be found in Folder 1 ('1_Data' folder). Folders 2 through 5 specify analyses for each component of the predation index. Folder 6 contains predation indices (all predators combined) for the area encompassed by the main stock assessment model for Gulf of Alaska pollock, culminating in Figure 1 (Dorn and Barnes, in revision). 

All analyses were conducted using R v4.1 (R Core Team 2021).

## Data Sources
<b>Total Predator Biomass</b>: Total biomass estimates were obtained from the most recent stock assessment for each groundfish predator (Barbeaux et al. 2019, Dorn et al. 2019, Hanselman et al. 2019, Spies et al. 201, Stewart and Hicks 2020). Coast-wide estimates for Pacific Halibut were adjusted to reflect biomass in the Gulf of Alaska.

<b>Relative Predator Densities</b>: Bottom trawl survey data (all groundfish predators) were collected by the Resource Assessment and Conservation Engineering (RACE) Division of the Alaska Fisheries Science Center (AFSC, NOAA) and are publicly accessible at https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. See von Szalay et al. (2016) for information about bottom trawl survey design and data collection methods. Setline survey data (Pacific Halibut; 1998 to 2019) were collected by the International Pacific Halibut Commission and are publicly available at: https://iphc.int/data/fiss-data-query. For setline survey methods, see Clark and Hare (2006). Longline survey data (Sablefish; 1990 to 2019) were collected by the AFSC's Auke Bay Laboratories and can be found at https://www.afsc.noaa.gov/maps/longline/Map.php. See Sigler and Zenger (1989) for methods descriptions of the longline survey. 

<b>Mean Annual Rations and Age-specific Proportions of Pollock Consumed</b>: Food habits data (all groundfish predators; 1990 to 2019) were provided by the AFSC's Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible at: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php. For food habits data collection and processing methods, see Livingston et al. (2017).  

## Financial and Logistical Support
This project was funded by the Pollock Conservation Cooperative Research Center (G00009488) and the Rasmuson Fisheries Research Center associated with the University of Alaska Fairbanks. An anonymous donor supplied additional funds via the Northern Gulf of Alaska Applied Research Award. The University of Alaska (Juneau Fisheries Division and Southeast Sitka Campus) provided facilities and additional support. 

## Acknowledgments
Discussions with Anne Beaudreau helped guide the direction of this research. We appreciate assistance with data acquisition and processing from Kerim Aydin, Steve Barbeaux, Troy Buckley, Dana Hanselman, Tom Kong, Ned Laman, Geoff Lang, Wayne Palsson, and Ian Stewart. Paul Spencer, Cody Szuwalski, Anne Hollowed, and two anonymous reviewers provided valuable comments to improve upon an earlier draft of this manuscript. <br>

## References 
### Methods
&#8212; Barnes CL, Beaudreau AH, Dorn MW, Holsman KK, and Mueter FJ. 2020. Development of a predation index to assess trophic stability in the Gulf of Alaska. Ecol Appl. 30(7):e02141.<br>

### Stock Assessments
&#8212; Barbeaux S, K Aydin, B Fissel, K Holsman, B Laurel, W Palsson, L Rogers, K Shotwell, Q Yang, and S Zador. 2019. Assessment of the Pacific cod stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report. 140pp.<br>
&#8212; Clark WG and SR Hare. 2006. Assessment and management of Pacific halibut: data, methods, and policy. IPHC Scientific Report 83. <br>
&#8212; Dorn M, AL Deary, BE Fissell, DT Jones, NE Lauffenburger, W.A. Palsson, LA Rogers, SK Shotwell, KA Spalinger, and S Zador. 2019. Assessment of the Walleye Pollock stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report. 161pp.<br>
&#8212; Hanselman DH, CJ Rodgveller, KH Fenske, SK Shotwell, KB Echave, PW Malecha, and CR Lunsford. 2019. Assessment of the Sablefish stock in Alaska. North Pacific Fishery Management Council Bering Sea, Aleutian Islands, and Gulf of Alaska SAFE Report. 263pp.<br>
&#8212; Spies I, K Aydin, JN Ianelli, and W Palsson. 2019. Assessment of the Arrowtooth Flounder stock in the Gulf of Alaska. North Pacific Fishery Management Council Gulf of Alaska SAFE Report. 92pp.<br>
&#8212; Stewart I and A Hicks. 2020. Assessment of the Pacific halibut (Hippoglossus stenolepis) stock at the end of 2019. International Pacific Halibut Commission Report. 32pp.<br>
### Survey and Food Habits Data
&#8212; Clark WG and SR Hare. 2006. Assessment and management of Pacific halibut: data, methods, and policy. IPHC Scientific Report 83. <br> 
&#8212; Livingston PA, K Aydin, TW Buckley, GM Lang, M-S Yang, and BS Miller. 2017. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fish. 100(4):443–470. <br>
&#8212; Sigler MF and HH Zenger, Jr. 1989. Assessment of Gulf of Alaska Sablefish and other groundfish based on the domestic longline survey, 1987. NOAA Techn Mem NMFS-AFSC Report 169. <br>
&#8212; von Szalay PG and NW Raring. 2016. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle, WA. NOAA Techn Mem NMFS-AFSC-325. <br>
### Species Distribution Modeling
&#8212; Barnes CL, AH Beaudreau, ME Hunsicker, and L Ciannelli (2018). Assessing the potential for competition between Pacific Halibut (<i>Hippoglossus stenolepis</i>) and Arrowtooth Flounder (<i>Atheresthes stomias</i>) in the Gulf of Alaska. PLoS ONE 13(12):e0209402. <br>
&#8212; Hunsicker ME, L Ciannelli, KM Bailey, S Zador, and L Stige. 2013. Climate and demography dictate the strength of predator-prey overlap in a subarctic marine ecosystem. PLoS ONE 8(6):e66025. <br>
&#8212; Shelton AO, ME Hunsicker, EJ Ward, BE Feist, R Blake, CL Ward, et al. 2017. Spatio-temporal models reveal subtle changes to demersal communities following the Exxon Valdez oil spill. ICES J Mar Sci. doi:10.1093/icesjms/fsx079. <br>
### Bioenergetics
&#8212; Armstrong JB and Schindler DE. 2011. Excess digestive capacity in predators reflects a life of feast and famine. Nature. 476:84–87. <br>
&#8212; Beaudreau AH and TE Essington. 2009. Development of a new field-based approach for estimating consumption rates of fishes and comparison with a bioenergetics model for lingcod (<i>Ophiodon elongatus</i>). Can J Fish Aquat Sci. 66:565−578. <br>
&#8212; Harvey CJ 2009. Effects of temperature change on demersal fisheries in the California Current: a bioenergetics approach. Can J Fish Aquat Sci. 66:1449–1461. <br>
&#8212; Holsman KK and K Aydin. 2015. Comparative methods for evaluating climate change impacts on the foraging ecology of Alaskan groundfish. Mar Ecol Progr Ser. 521:217–235. <br>
&#8212; Holsman KK, K Aydin, J Sullivan, T Hurst, and G Kruse. 2019. Climate effects and bottom-up controls on growth and size-at-age of Pacific halibut (<i>Hippoglossus stenolepis</i>) in Alaska (USA). Fish Oceanogr. 28:345–358. <br>
### Miscellaneous
&#8212; Brodziak J. 2012. Fitting length-weight relationships with linear regression using the log-transformed allometric model with bias-correction. NOAA Technical Memorandum PIFSC-H-12-03. <br>
&#8212; Chipps SR and JE Garvey. 2007. Assessment of diets and feeding patterns. In: Analysis and interpretation of freshwater fisheries data. CS Guy and ML Brown, eds. Bethesda, MD. Amer Fish Soc. 473–514. <br>
&#8212; R Core Team. 2021. R: a language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.
