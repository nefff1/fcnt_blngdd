################################################################################.
# SETUP ------------------------------------------ #############################
################################################################################.

library(tidyverse); theme_set(theme_classic())
library(cowplot)
library(FD)
library(parallel)
library(vegan)
library(glmmTMB)

select <- dplyr::select
rename <- dplyr::rename

# define global variables ######################################################
################################################################################.

# define tree names and order
tree.order <- c(BU = "Fagus", AH = "Acer", BI = "Betula", HBU = "Carpinus", 
                PA = "Populus", LI = "Tilia",
                ES = "Fraxinus", EI = "Quercus", LA = "Larix", FI = "Picea",
                KI = "Pinus", DGL = "Pseudotsuga")


# plotting order of filtering variables
var.order <- c("prop_forest2000", "DWV", "NMDS1", "NMDS2", "ssc_mean",
               "Canopy_closure", "Cover_above2m", "Period1")

# filtering variables at site level
vars.site <- c("prop_forest2000", "DWV", "NMDS1", "NMDS2", "ssc_mean",
               "Canopy_closure")

# filtering variables at patch level
vars.patch <- c("DWV", "ssc_mean", "Cover_above2m")

# read data ####################################################################
################################################################################.

# beetle data from deadwood eclectors ------------------------------------------.
d.abund.stem <- read.csv("Data/BEXIS/31123_16_Dataset/31123_16_data.csv")

d.abund.stem <- d.abund.stem %>% 
  mutate(PlotID = substr(TrapID, 1, 5)) %>% # to adjust new data set
  rename(ToplaAdults = Adults) %>% # to adjust new data set
  mutate(PlotID = as.character(PlotID),
         SpeciesID = as.character(SpeciesID)) %>% 
  filter(ToplaAdults > 0)

# exclude Prunus (only present at 20 plots) and plot not sampled in other samplings
d.abund.stem <- d.abund.stem %>% 
  filter(substr(TreeID, 1, 2) != "KB",
         PlotID != "HEW42")

# edit data
d.abund.stem <- d.abund.stem %>% 
  filter(TreeID != "DGL080") %>%
  # correct two mistakes and one taxonomic inconsistency
  mutate(PlotID = ifelse(TreeID == "FI061", "SEW08", PlotID),
         PlotID = ifelse(TreeID == "DGL012", "HEW08", PlotID),
         SpeciesID = ifelse(SpeciesID == "Phloeopora bernhaueri", 
                            "Phloeopora corticalis", SpeciesID))


# beetle data flight interception traps (ground level) -------------------------.
d.abund.ffb <- read.csv("Data/BEXIS/22007_5_Dataset/22007_5_data.csv", sep = ";")

d.abund.ffb <- d.abund.ffb %>% 
  mutate(CollectionYear = as.Date(CollectionYear, origin = "1899-12-30"),
         CollectionYear = substr(CollectionYear, 1, 4),
         CollectionYear = as.integer(CollectionYear)) %>% 
  mutate(PlotID = ifelse(nchar(PlotID) == 4,
                         paste0(substr(PlotID, 1, 3), "0", substr(PlotID, 4, 4)),
                         PlotID)) %>% 
  rename(SpeciesID = Species) %>% 
  filter(Order == "Coleoptera",
         NumberAdults > 0,
         !PlotID %in% c("HEW10", "HEW13"))


# beetle data flight interception traps (crown level) --------------------------.
d.abund.ffk <- read.csv("Data/BEXIS/22026_2_Dataset/22026_2_data.csv", sep = ";")

d.abund.ffk <- d.abund.ffk %>% 
  mutate(PlotID = ifelse(nchar(PlotID) == 4,
                         paste0(substr(PlotID, 1, 3), "0", substr(PlotID, 4, 4)),
                         PlotID),
         NumberAdults = as.integer(NumberAdults)) %>% 
  filter(Order == "Coleoptera",
         NumberAdults > 0,
         !PlotID %in% c("HEW10", "HEW13")) 

# beetle data pitfall traps ----------------------------------------------------,
d.abund.bf <- read.csv("Data/BEXIS/21968_2_Dataset/21968_2_data.csv", sep = ";")

d.abund.bf <- d.abund.bf %>% 
  mutate(PlotID = ifelse(nchar(PlotID) == 4,
                         paste0(substr(PlotID, 1, 3), "0", substr(PlotID, 4, 4)),
                         PlotID),
         NumberAdults = as.integer(NumberAdults)) %>% 
  filter(Order == "Coleoptera",
         NumberAdults > 0,
         !PlotID %in% c("HEW10", "HEW13"))

# beetle morphometric trait data -----------------------------------------------.

d.traits.m <- read.csv("Data/BEXIS/31155_3_Dataset/31155_3_data.csv", sep = ";")

d.traits.m <- d.traits.m %>% 
  rename_at(vars(-c(ID, Order, Suborder, Family, Species, Author, Source,
                    Voucher_ID, BE_Plot_ID, BE_Traptype, BE_CollectionDate,
                    Sampling_coordinates)),
            ~ tolower(.)) %>% 
  rename(wing_length = hind_wing_length,
         wing_width = hind_wing_width,
         wing_area = hind_wing_area) %>% 
  mutate(sex = ifelse(sex == "f", "female",
                      ifelse(sex == "m", "male", sex)))

# similar species (to estimate missing traits) ---------------------------------.
d.similar.species <- read.csv("Data/Similar_Species.csv")

# body length data (from FHL) --------------------------------------------------.
d.bodylength <- read.csv("Data/Body_length.csv")

# lightness data ---------------------------------------------------------------.

# lightness derived from Dries (2016)
d.lightness.dries <- read.csv("Data/Hagge_et_al/Hagge_etal_saproxylic_beetle_morphological_trait_database_20210426.csv")

# lightness derived from measured individuals and from coleonet.de (Lompe 2002)
d.lightness.plus <- read.csv("Data/Lightness_completion.csv")


# list of saproxylic species (subset based on study organisms) -----------------.
d.xylo <- read.csv("Data/Saproxylic_species_sub.csv")


# tree species composition data ------------------------------------------------.
d.tree1 <- read.csv("Data/BEXIS/18270_2_Dataset/18270_2_data.csv",
                    sep = ";") %>% 
  rename(BA = Basal_Area,
         PlotID = EP) %>% 
  select(-c(Ep_siteid, Sampling_from, Sampling_till)) %>% 
  mutate(period = 1)

d.tree2 <- read.csv("Data/BEXIS/22907_2_Dataset/22907_2_data.csv", sep = ";") %>%
  rename(PlotID = EP) %>% 
  select(-Ep_siteid) %>% 
  mutate(period = 2)

# Canopy cover subplot-level (patch) -------------------------------------------.
d.canopy.sp <- read.csv("Data/BEXIS/30925_2_Dataset/30925_2_data.csv", sep = ";")
d.normalplots <- read.csv("Data/BELongDead_Subplots_Normal.csv")

d.canopy.sp <- d.canopy.sp %>% 
  select(-Plotid) %>% 
  rename(PlotID = EP,
         Subplot = SUBPLOT) %>% 
  left_join(mutate(d.normalplots, x = 1), by = c("PlotID", "Subplot")) %>% 
  filter(x == 1) %>%
  select(-c(x, Subplot))

# Deadwood amount --------------------------------------------------------------.

# two inventories
d.deadwood.1 <- read.csv("Data/BEXIS/24546_2_Dataset/24546_2_data.csv", sep = ";")
d.deadwood.2 <- read.csv("Data/BEXIS/24526_2_Dataset/24526_2_data.csv", sep = ";")

# summarise total deadwood amount
d.deadwood.1 <- d.deadwood.1 %>% 
  group_by(plotID) %>% 
  summarise(DWV = sum(DWV, na.rm  = T),
            .groups = "drop")

d.deadwood.2 <- d.deadwood.2 %>% 
  filter(stratum != "T7") %>% # stratum T7 was only included in the second inventory
  group_by(plotID) %>% 
  summarise(DWV = sum(DWV, na.rm  = T),
            .groups = "drop")

# combine the two data sets

d.deadwood <- d.deadwood.1 %>% 
  bind_rows(d.deadwood.2) %>% 
  mutate(PlotID = ifelse(nchar(plotID) == 4,
                         paste0(substr(plotID, 1, 3), "0", substr(plotID, 4, 4)),
                         plotID)) %>% 
  group_by(PlotID) %>% 
  summarise(DWV = mean(DWV), 
            .groups = "drop")

# stand structural complexity --------------------------------------------------.
d.ssc <- read.csv("Data/BEXIS/20607_2_Dataset/20607_2_data.csv", sep = ";")

d.ssc <- d.ssc %>% 
  rename(PlotID = EP) %>% 
  mutate(PlotID = ifelse(nchar(PlotID) == 4,
                         paste0(substr(PlotID, 1, 3), "0", substr(PlotID, 4, 4)),
                         PlotID)) %>% 
  arrange(PlotID)

# canopy closure ---------------------------------------------------------------.
d.canopy <- read.csv("Data/BEXIS/22927_2_Dataset/22927_2_data.csv", sep = ";")

d.canopy <- d.canopy %>% 
  rename(PlotID = EP) %>% 
  select(-c(Exploratory, Ep_siteid))



# landscape --------------------------------------------------------------------.
d.lands <- read.csv("Data/BEXIS/26007_2_Dataset/26007_2_data.csv", sep = ";")

d.lands <- d.lands %>%
  mutate(prop_forest2000 = F_DLM_2000 / 100) %>% 
  group_by(PlotID) %>% 
  summarise(prop_forest2000 = mean(prop_forest2000),
            .groups = "drop")

# EP data ----------------------------------------------------------------------.
d.EP <- read.csv("Data/BEXIS/20826_6_Dataset/20826_6_data.csv", sep = ";")
d.EP <- d.EP %>%
  rename(PlotID_long = PlotID) %>%
  rename(PlotID = Ep_siteID) %>%
  mutate(PlotID = ifelse(nchar(PlotID) == 4,
                         paste0(substr(PlotID, 1, 3), "0", substr(PlotID, 4, 4)),
                         PlotID)) %>%
  arrange(PlotID)


# Taxonomic name standardisation -----------------------------------------------.
d.std <- read.csv("Data/Name_standardisation.csv", row.names = "X")

# name standardisation #########################################################
################################################################################.

d.abund.stem <- d.abund.stem %>%
  select(-c(Order, EDV_CODE, Suborder, Family)) %>% 
  left_join(d.std, by = c("SpeciesID" = "Name_original"))

d.abund.ffb <- d.abund.ffb %>% 
  select(-c(Order, Suborder, Family)) %>% 
  left_join(d.std, by = c("SpeciesID" = "Name_original"))

d.abund.ffk <- d.abund.ffk %>% 
  select(-c(Order, Suborder, Family)) %>% 
  left_join(d.std, by = c("Species" = "Name_original"))

d.abund.bf <- d.abund.bf %>% 
  select(-c(Order, Suborder, Family)) %>% 
  left_join(d.std, by = c("Species" = "Name_original"))

d.traits.m <- d.traits.m %>% # data set already has standardised names
  rename(Name_std = Species)

d.xylo <- d.xylo %>% 
  left_join(d.std, by = c("species" = "Name_original"))

# exclude non-saproxylic species ###############################################
################################################################################.

d.abund.stem <- d.abund.stem %>% 
  filter(Name_std %in% d.xylo$Name_std)

d.abund.ffb <- d.abund.ffb %>% 
  filter(Name_std %in% d.xylo$Name_std)

d.abund.ffk <- d.abund.ffk %>% 
  filter(Name_std %in% d.xylo$Name_std) 

d.abund.bf <- d.abund.bf %>% 
  filter(Name_std %in% d.xylo$Name_std)

sp_list <- unique(c(d.abund.stem$Name_std, d.abund.ffb$Name_std, 
                    d.abund.ffk$Name_std, d.abund.bf$Name_std))

################################################################################.
# EDIT TRAITS --------------------------------- ################################
################################################################################.

# missing part estimation ######################################################
################################################################################.

# track number of estimates / what is estimated
n_estimates <- 0
d.estimates <- data.frame()

# 1st: estimate missing leg parts based on conspecifics and congenerics --------.
# if no conspecifics and congenerics available, switch to similar genera

missing_legparts <- d.traits.m %>%
  filter_at(vars(contains("femur_length"),
                 contains("tibia_length")),
            any_vars(is.na(.)))

# initiate data frame, in which corrections are made
# this is necessary to not take corrected values for further estimations, which would cause a bias
d.traits.m_corr <- d.traits.m
d.traits.m_corr$estimation <- F

for (i in 1:nrow(missing_legparts)){
  sp_i <- missing_legparts$Name_std[i] # species of observation
  sex_i <- missing_legparts$sex[i] # sex of observation
  voucher_i <- missing_legparts$ID[i] # voucher of observation
  genus_i <- missing_legparts$Genus[i] # species of observation
  
  
  missing_i <- missing_legparts[i, ] %>%  # missing leg traits of observation
    select(contains("femur_length"), contains("tibia_length")) %>% 
    select(which(is.na(.))) %>% 
    names()
  nonmissing_i <- missing_legparts[i, ] %>%  # NON missing leg traits of observation
    select(contains("femur_length"), contains("tibia_length")) %>% 
    select(which(!is.na(.))) %>% 
    names()
  
  # only do the replacements when at least some leg parts were measured
  if (length(missing_i) < 5){ 
    
    # loop through all missing traits
    for (trait_i in missing_i){
      
      subset <- d.traits.m[NULL, ]
      
      # subset observations of same species and sex
      if (!is.na(sex_i)){
        subset <-  d.traits.m %>%
          filter(ID != voucher_i,
                 Name_std == sp_i,
                 sex == sex_i | is.na(sex),
                 !is.na(!! sym(trait_i)))
      }
      
      # if sex is not given or if there are no observations of the same sex
      if (is.na(sex_i) | nrow(subset) == 0){
        subset <- d.traits.m %>% 
          filter(ID != voucher_i,
                 Name_std == sp_i,
                 !is.na(!! sym(trait_i)))
      }
      
      # if there are still no other specimens, take other species from genus
      if (nrow(subset) == 0){
        subset <- d.traits.m %>% 
          filter(ID != voucher_i,
                 Genus == genus_i,
                 !is.na(!! sym(trait_i)))
      }
      
      # if there are still no other specimens, take similar species
      if (nrow(subset) == 0){
        subset <- d.traits.m %>% 
          filter(ID != voucher_i,
                 Name_std %in% d.similar.species$Similar_Species[d.similar.species$Missing_SpeciesID == sp_i],
                 !is.na(!! sym(trait_i)))
      }
      
      if (nrow(subset) > 0){
        ratios <- data.frame(nonmissing = nonmissing_i) # initiate data frame for ratios
        
        # for each missing trait, calculate for other observations the ratios between measured leg traits and the missing trait
        ratios <- subset %>% 
          rename(missingtrait = trait_i) %>% 
          select(ID, !! enquo(nonmissing_i), missingtrait) %>% 
          gather(ratiotrait, traitvalue, -c(ID, missingtrait)) %>% 
          mutate(ratio =  traitvalue / missingtrait) %>% 
          group_by(ratiotrait) %>% 
          summarise(ratio = mean(ratio, na.rm = T),
                    .groups = "drop") # take mean across all specimens
        
        # now estimate trait value from all the different other leg traits that were measured
        estimate <- missing_legparts[i, ] %>% 
          select(one_of(ratios$ratiotrait)) %>% 
          gather(ratiotrait, traitvalue) %>% 
          left_join(ratios, by = "ratiotrait") %>% 
          mutate(estimate = traitvalue / ratio) %>% 
          summarise(estimate = mean(estimate, na.rm = T)) %>% 
          deframe()
        
        # make changes in the corrected data frame
        d.traits.m_corr[d.traits.m_corr$ID == voucher_i,
                        trait_i] <- estimate
        rm(estimate)
        
        if (sp_i %in% sp_list) {
          n_estimates <- n_estimates + 1 # track changes
          
          d.estimates <- data.frame(Estimated_for = sp_i,
                                    Estimated_traits = "legs",
                                    Estimated_from = unique(subset$Name_std)) %>% 
            bind_rows(d.estimates, .)
        }
        
        d.traits.m_corr$estimation[d.traits.m_corr$ID == voucher_i] <- T
        
      } else { # if no similar species, take other legs of same specimen
        
        # ratios of all leg pairs
        legratios <- missing_legparts[i, ] %>% 
          summarise(ratio1 = front_femur_length / front_tibia_length,
                    ratio2 = mid_femur_length / mid_tibia_length,
                    ratio3 = hind_femur_length / hind_tibia_length) %>% 
          unlist()
        
        if (grepl("tibia", trait_i)){
          part2 <- paste0(substr(trait_i, 1, regexpr("_", trait_i) - 1),
                          "_femur_length")
          
          estimate <- missing_legparts[i, part2] / mean(legratios, na.rm = T)
        } else {
          part2 <- paste0(substr(trait_i, 1, regexpr("_", trait_i) - 1),
                          "_tibia_length")
          
          estimate <- missing_legparts[i, part2] * mean(legratios, na.rm = T)
        }
        d.traits.m_corr[d.traits.m_corr$ID == voucher_i,
                        trait_i] <- estimate
        rm(estimate)
        
        if (sp_i %in% sp_list) n_estimates <- n_estimates + 1 # track changes
        
        d.traits.m_corr$estimation[d.traits.m_corr$ID == voucher_i] <- T
      }
    }
  }
}


d.traits.m_corr <- d.traits.m_corr %>% 
  select(-estimation)
d.traits.m <- d.traits.m_corr; rm(d.traits.m_corr)

# 2nd: estimate missing body dimensions based on conspecifics and congenerics --.

# only for cases in which ONE body dimension is missing
missing_bodydims <- d.traits.m %>% 
  filter(Name_std %in% c(d.abund.stem$Name_std, d.abund.ffb$Name_std,
                         d.abund.ffk$Name_std, d.abund.bf$Name_std,
                         d.similar.species$Similar_Species)) %>%
  filter(is.na(body_length) + is.na(body_width) + is.na(body_height) == 1)


# initiate data frame, in which corrections are made
# this is necessary to not take corrected values for further estimations, which would cause a bias
d.traits.m_corr <- d.traits.m
d.traits.m_corr$estimation <- F

for (i in 1:nrow(missing_bodydims)){
  sp_i <- missing_bodydims$Name_std[i] # species of observation
  sex_i <- missing_bodydims$sex[i] # sex of observation
  voucher_i <- missing_bodydims$ID[i] # voucher of observation
  genus_i <- missing_bodydims$Genus[i] # species of observation
  
  
  missing_i <- missing_bodydims[i, ] %>%  # missing body dimension of observation
    select(contains("body_")) %>% 
    select(which(is.na(.))) %>% 
    names()
  nonmissing_i <- missing_bodydims[i, ] %>%  # NON missing body dimension of observation
    select(contains("body_")) %>% 
    select(which(!is.na(.))) %>% 
    names()
  
  
  subset <- d.traits.m[NULL, ]
  
  # subset observations of same species and sex
  if (!is.na(sex_i)){
    subset <-  d.traits.m %>%
      filter(ID != voucher_i,
             Name_std == sp_i,
             sex == sex_i | is.na(sex),
             !is.na(!! sym(missing_i)))
  }
  
  # if sex is not given or if there are no observations of the same sex
  if (is.na(sex_i) | nrow(subset) == 0){
    subset <- d.traits.m %>% 
      filter(ID != voucher_i,
             Name_std == sp_i,
             !is.na(!! sym(missing_i)))
  }
  
  # if there are still no other specimens, take other species from genus
  if (nrow(subset) == 0){
    subset <- d.traits.m %>% 
      filter(ID != voucher_i,
             Genus == genus_i,
             !is.na(!! sym(missing_i)))
  }
  
  # if there are still no other specimens, take similar species
  if (nrow(subset) == 0){
    subset <- d.traits.m %>% 
      filter(ID != voucher_i,
             Name_std %in% d.similar.species$Similar_Species[d.similar.species$Missing_SpeciesID == sp_i],
             !is.na(!! sym(missing_i)))
  }
  
  if (nrow(subset) > 0){
    ratios <- data.frame(nonmissing = nonmissing_i) # initiate data frame for ratios
    
    # for each missing trait, calculate for other observations the ratios between measured traits and the missing trait
    ratios <- subset %>% 
      rename(missingtrait = missing_i) %>% 
      select(ID, !! enquo(nonmissing_i), missingtrait) %>% 
      gather(ratiotrait, traitvalue, -c(ID, missingtrait)) %>% 
      mutate(ratio =  traitvalue / missingtrait) %>% 
      group_by(ratiotrait) %>% 
      summarise(ratio = mean(ratio, na.rm = T),
                .groups = "drop") # take mean across all specimens
    
    # now estimate trait value from all the different other leg traits that were measured
    estimate <- missing_bodydims[i, ] %>% 
      select(one_of(ratios$ratiotrait)) %>% 
      gather(ratiotrait, traitvalue) %>% 
      left_join(ratios, by = "ratiotrait") %>% 
      mutate(estimate = traitvalue / ratio) %>% 
      summarise(estimate = mean(estimate, na.rm = T)) %>% 
      deframe()
    
    # make changes in the corrected data frame
    d.traits.m_corr[d.traits.m_corr$ID == voucher_i,
                    missing_i] <- estimate
    rm(estimate)
    
    if (sp_i %in% sp_list) {
      n_estimates <- n_estimates + 1 # track changes
      d.estimates <- data.frame(Estimated_for = sp_i,
                                Estimated_traits = "body dimensions",
                                Estimated_from = unique(subset$Name_std)) %>% 
        bind_rows(d.estimates, .)
    }
    
    d.traits.m_corr$estimation[d.traits.m_corr$ID == voucher_i] <- T
    
  } else {
    warning(paste(sp_i, "has no similar species."))
  }
}

d.traits.m_corr <- d.traits.m_corr %>% 
  select(-estimation)
d.traits.m <- d.traits.m_corr; rm(d.traits.m_corr)

# missing mass estimation ######################################################
################################################################################.

# for cases with missing mass, estimate mass from body volume
d.traits.m <- d.traits.m %>% 
  mutate(body.volume = body_length * body_width * body_height)

mod_mass <- lm(log(mass) ~ log(body.volume), d.traits.m)

missing_mass <- d.traits.m %>% 
  filter(Name_std %in% c(d.abund.stem$Name_std, d.abund.ffb$Name_std,
                         d.abund.ffk$Name_std, d.abund.bf$Name_std)) %>% 
  filter(is.na(mass))

# predict mass values
missing_mass$estimation <- exp(predict(mod_mass, missing_mass))

# make changes
for (voucher_i in missing_mass$voucher){
  d.traits.m$mass[d.traits.m$ID == voucher_i] <- missing_mass$estimation[missing_mass$ID == voucher_i]
}

n_estimates <- n_estimates + sum(missing_mass$Name_std %in% sp_list) # track

# calculate processed traits ###################################################
################################################################################.

d.traits.m <- d.traits.m %>% 
  rowwise() %>% 
  mutate(
    leg_length = mean(c(front_femur_length + front_tibia_length,
                        mid_femur_length + mid_tibia_length,
                        hind_femur_length + hind_tibia_length),
                      na.rm = T),
    wing.aspect = wing_length / wing_width, 
    wing.loading.rev = wing_area / mass,
    body.roundness = body_height / body_width,
    mandibular.aspect = inlever_length / mandible_length,
    hairiness = hairiness_pronotum_count / hairiness_pronotum_line)  %>% 
  ungroup()

# length regressions to standardise several absolute traits --------------------.

v.std.traits <- c("body_width", "head_length", "wing_length", "leg_length", 
                  "antenna_length", "eye_length")

for (tr_i in v.std.traits){
  dat_target <- d.traits.m
  dat_target$response <-  dat_target[, tr_i] %>% deframe
  
  # if there are 0s in the response, those are set to NA (see below why)
  if (any(dat_target$response == 0, na.rm = T)){
    warning(paste("Zeros in", tr_i))
    dat_target$response[dat_target$response == 0] <- NA
  }
  
  dat_target <- dat_target %>% 
    filter(!is.na(response), !is.na(body_length))
  
  mod_target <- lm(log(response) ~ log(body_length), dat_target)
  print(paste(tr_i, paste(coef(mod_target), collapse = "/")))
  
  res_trait <- paste0(gsub("_", ".", tr_i), ".res")
  
  dat_target[, res_trait] <- residuals(mod_target)
  
  d.traits.m <- dat_target %>% 
    select(ID, !! sym(res_trait)) %>% 
    left_join(d.traits.m, ., by = "ID")
  
  # if there are 0s in the response (wing length!), set relative values to the observed minimum
  d.traits.m <- d.traits.m %>% 
    mutate(!! sym(res_trait) := ifelse(!! sym(tr_i) == 0, 
                                       min(!! sym(res_trait), na.rm = T),
                                       !! sym(res_trait)))
}


v.sel.traits <- c("body.width.res", "body.roundness", "head.length.res", 
                  "eye.length.res", "antenna.length.res", "hairiness",
                  "leg.length.res", "wing.length.res", "wing.loading.rev", 
                  "wing.aspect", "mandibular.aspect")

# trait transformations ########################################################
################################################################################.

min.nonzero <- min(d.traits.m$hairiness[d.traits.m$hairiness != 0], na.rm = T)
d.traits.m$hairiness <- log(d.traits.m$hairiness + min.nonzero)

min.nonzero <- min(d.traits.m$wing.loading.rev[d.traits.m$wing.loading.rev != 0], na.rm = T)
d.traits.m$wing.loading.rev <- log(d.traits.m$wing.loading.rev + min.nonzero)

# aggregate traits #############################################################
################################################################################.

# use all species, even non-saproxylics, because might be important for missing-part estimations
v.species <- d.traits.m %>%
  select(Name_std) %>% 
  distinct() %>% 
  deframe()

f.agg.sex <- function(sp_i){

    target <- d.traits.m %>% 
    filter(Name_std == sp_i)
  
  if (nrow(target) == 1){ # if only one measured specimen
    out <- target %>% 
      select(Name_std, !! enquo(v.sel.traits))
  } else {
    if (all(is.na(target$sex))){
      out <- target %>% 
        select(Name_std, !! enquo(v.sel.traits)) %>% 
        group_by(Name_std) %>% 
        summarise_all(~mean(., na.rm = T)) %>% 
        ungroup()
    } else {
      if (!all(target$sex %in% c("male", "female") | is.na(target$sex))) warning(paste("Wrong sex in", sp_i))
      
      if (!any(is.na(target$sex))){
        out <- target %>% 
          select(Name_std, sex, !! enquo(v.sel.traits)) %>% 
          group_by(Name_std, sex) %>% 
          summarise_all(~mean(., na.rm = T)) %>% 
          ungroup() %>% 
          select(-sex) %>% 
          group_by(Name_std) %>% 
          summarise_all(~mean(., na.rm = T)) %>% 
          ungroup()
      } else {
        if (all(c("male", "female") %in% target$sex)){
          # if there is males, females, and NAs in the data
          out <- target %>% 
            mutate(sex = ifelse(is.na(sex), "female", sex)) %>% 
            # troick to distribute NA cases to both sexes but weight them only half
            bind_rows(target %>% 
                        mutate(sex = ifelse(is.na(sex), "male", sex))) %>% 
            select(Name_std, sex, !! enquo(v.sel.traits)) %>% 
            group_by(Name_std, sex) %>% 
            summarise_all(~mean(., na.rm = T)) %>% 
            ungroup() %>% 
            select(-sex) %>% 
            group_by(Name_std) %>% 
            summarise_all(~mean(., na.rm = T)) %>% 
            ungroup()
        } else {
          # if there is only one sex and NAs, attribute all NAs to the missing sex
          out <- target %>% 
            mutate(sex_adj = ifelse(is.na(sex),
                                    ifelse("female" %in% sex, "male", "female"), 
                                    sex)) %>% 
            select(-sex) %>% 
            rename(sex = sex_adj) %>% 
            select(Name_std, sex, !! enquo(v.sel.traits)) %>% 
            group_by(Name_std, sex) %>% 
            summarise_all(~mean(., na.rm = T)) %>% 
            ungroup() %>% 
            select(-sex) %>% 
            group_by(Name_std) %>% 
            summarise_all(~mean(., na.rm = T)) %>% 
            ungroup()
          
        }
      }
      
    }
  }
  out
}

d.traits.m.agg <- lapply(v.species, f.agg.sex) %>% 
  do.call(rbind, .)

# estimate missing values in processed traits ##################################
################################################################################.

# track changes
n_estimates_proc <- 0

missing_parts <- d.traits.m.agg %>% 
  filter(Name_std %in% c(d.abund.stem$Name_std, d.abund.ffb$Name_std,
                         d.abund.ffk$Name_std, d.abund.bf$Name_std,
                         d.similar.species$Similar_Species)) %>% 
  mutate(n.NA = rowSums(is.na(.)),
         n.NA = ifelse((wing.length.res == min(wing.length.res, na.rm = T)) & !is.na(wing.length.res), 
                       n.NA - 1, n.NA)) %>% #wing.aspect is NA per default for wing.length == 0
  filter(n.NA > 0) %>% 
  left_join(d.std %>% select(Name_std, Genus) %>% distinct(), by = "Name_std")

# initiate data frame, in which corrections are made
# this is necessary to not take corrected values for further estimations, which would cause a bias
d.traits.m.agg_corr <- d.traits.m.agg
d.traits.m.agg_corr$estimation <- F

for (i in 1:nrow(missing_parts)){
  sp_i <- missing_parts$Name_std[i] # species of observation
  genus_i <- missing_parts$Genus[i] # species of observation
  
  
  missing_i <- missing_parts[i, ] %>%  # missing measurements
    select(which(is.na(.))) %>% 
    names()
  
  # exclude wing shape from missing measurements if wing length is zero
  if ("wing.shape" %in% missing_i & (missing_parts$wing.length.res[i] == min(missing_parts$wing.length.res, na.rm = T) & 
                                     !is.na(missing_parts$wing.length.res[i]))){
    missing_i <- missing_i[missing_i != "wing.shape"]
  }
  
  
  subset <- d.traits.m.agg %>% 
    left_join(d.std %>% select(Name_std, Genus) %>% distinct(), 
              by = "Name_std") %>% 
    filter(Genus == genus_i,
           Name_std != sp_i)
  
  
  # if there are no other species of the genus, take similar species
  if (nrow(subset) == 0){
    subset <- d.traits.m.agg %>%
      filter(Name_std != sp_i,
             Name_std %in% d.similar.species$Similar_Species[d.similar.species$Missing_SpeciesID == sp_i])
  }
  
  if (nrow(subset) > 0){
    
    subset_agg <- subset %>% 
      select(!! enquo(v.sel.traits)) %>% 
      summarise_all(~mean(., na.rm = T))
    
    for (ms_i in missing_i){
      d.traits.m.agg_corr[d.traits.m.agg_corr$Name_std == sp_i, ms_i] <- subset_agg[, ms_i]
      
      if (sp_i %in% sp_list) {
        n_estimates_proc <- n_estimates_proc + 1 # track changes
        d.estimates <- data.frame(Estimated_for = sp_i,
                                  Estimated_traits = ms_i,
                                  Estimated_from = unique(subset$Name_std)) %>% 
          bind_rows(d.estimates, .)
      }
    }
    
    d.traits.m.agg_corr$estimation[d.traits.m.agg_corr$Name_std == sp_i] <- T
  } else {
    warning(paste(sp_i,"has no similar species."))
  }
}


d.traits.m.agg <- d.traits.m.agg_corr %>% 
  select(-estimation)

# estimate antenna length for Calyptomerus alpestris:

# measured from FHL: Antenna length relative to body length = 0.22

# now convert this into the model residuals
est_antenna_missing <- d.traits.m %>% 
  filter(Name_std == "Calyptomerus alpestris") %>% 
  select(body_length) %>% 
  mutate(antenna_length = 0.22 * body_length,
         antenna_length_exp = exp(-1.38724635243472 + 
                                    1.15532824494014 * log(body_length))) %>% # from length regression before
  mutate(antenna.length.res = antenna_length - antenna_length_exp) %>% 
  summarise(antenna.length.res = mean(antenna.length.res)) %>% 
  deframe()

# make change in data frame:
d.traits.m.agg <- d.traits.m.agg %>% 
  mutate(antenna.length.res = ifelse(Name_std == "Calyptomerus alpestris",
                                     est_antenna_missing, antenna.length.res))

# add trait data for completely missing species ################################
################################################################################.

missing_measurements <- d.abund.stem %>% 
  filter(!Name_std %in% d.traits.m$Name_std) %>% 
  group_by(Name_std, FHLCode) %>% 
  summarise(n_stem = sum(ToplaAdults, na.rm = T)) %>%
  ungroup() %>% 
  full_join(d.abund.ffb %>% 
              filter(!Name_std %in% d.traits.m$Name_std) %>% 
              group_by(Name_std, FHLCode) %>% 
              summarise(n_ffb = sum(NumberAdults, na.rm = T)) %>% 
              ungroup,
            by = c("Name_std", "FHLCode")) %>%
  full_join(d.abund.ffk %>% 
              filter(!Name_std %in% d.traits.m$Name_std) %>% 
              group_by(Name_std, FHLCode) %>% 
              summarise(n_ffk = sum(NumberAdults, na.rm = T)) %>% 
              ungroup,
            by = c("Name_std", "FHLCode")) %>%
  full_join(d.abund.bf %>% 
              filter(!Name_std %in% d.traits.m$Name_std) %>% 
              group_by(Name_std, FHLCode) %>% 
              summarise(n_bf = sum(NumberAdults, na.rm = T)) %>% 
              ungroup,
            by = c("Name_std", "FHLCode")) %>%
  arrange(Name_std)


for (sp_i in missing_measurements$Name_std){
  estimate <- d.traits.m.agg %>% 
    filter(Name_std %in% d.similar.species$Similar_Species[d.similar.species$Missing_SpeciesID == sp_i]) %>% 
    mutate(Name_std = sp_i) %>% 
    group_by(Name_std) %>% 
    summarise_all(~mean(., na.rm = T)) %>% 
    ungroup() 
  
  d.estimates <- data.frame(Estimated_for = sp_i,
                            Estimated_traits = "all but length and lightness",
                            Estimated_from = unique(d.similar.species$Similar_Species[d.similar.species$Missing_SpeciesID == sp_i])) %>% 
    bind_rows(d.estimates, .)
  
  if (nrow(estimate) == 1) {
    d.traits.m.agg <- bind_rows(d.traits.m.agg, estimate)
  } else {
    warning(paste("Check", sp_i))
  }
}

d.traits.m.agg <- d.traits.m.agg %>% 
  filter(Name_std %in% c(d.abund.stem$Name_std, d.abund.ffb$Name_std,
                         d.abund.ffk$Name_std, d.abund.bf$Name_std))

# add literature length data ###################################################
################################################################################.

d.traits.m.agg <-
  d.bodylength %>%  
  left_join(d.std, by = c(SpeciesID = "Name_original")) %>% 
  group_by(Name_std) %>% 
  summarise(body.length = mean(Body_Size)) %>% # 1 duplicate because of name standardisation
  ungroup() %>% 
  mutate(body.length = log(body.length)) %>%  # log transformation
  left_join(d.traits.m.agg, ., by = "Name_std")

v.sel.traits <- c("body.length", v.sel.traits)

# add lightness data ###########################################################
################################################################################.

# name standardisation
d.lightness.dries <-
  d.lightness.dries %>% 
  filter(!is.na(colour_lightness)) %>% 
  mutate(species = gsub("_", " ", species)) %>%
  # two standardisation steps necessary because of critical synonyms
  left_join(d.std %>% select(FHLCode, Name_std) %>% distinct(), 
            by = c(edv_code = "FHLCode")) %>% 
  left_join(d.std %>% select(Name_original, Name_std), by = c(species = "Name_original"),
            suffix = c("", ".plus")) %>% 
  mutate(Name_std = ifelse(is.na(Name_std), Name_std.plus, Name_std)) %>% 
  select(-Name_std.plus)

d.lightness.plus <- d.lightness.plus %>% 
  left_join(d.std, by = c(SpeciesID = "Name_original"))

# determine SDs and means to bring the newly measured data set to the same scale
sd_dries <- sd(d.lightness.dries$colour_lightness, na.rm = T)
mean_dries <- mean(d.lightness.dries$colour_lightness, na.rm = T)

sd_plus <- sd(d.lightness.plus$colour_lightness, na.rm = T)
mean_plus <- mean(d.lightness.plus$colour_lightness, na.rm = T)

# change scale
d.lightness.plus <- d.lightness.plus %>% 
  mutate(colour_lightness = (colour_lightness - mean_plus) / sd_plus * sd_dries + mean_dries)

# restrict additional data to species not measured from Dries. Aggregate per species
d.lightness.plus <- d.lightness.plus %>% 
  filter(!Name_std %in% d.lightness.dries$Name_std) %>% 
  group_by(Name_std) %>% 
  summarise(colour_lightness = mean(colour_lightness, na.rm = T),
            .groups = "drop")

# combine the two data sets
d.lightness <- d.lightness.dries %>% 
  bind_rows(d.lightness.plus) %>% 
  select(Name_std, colour_lightness) %>% 
  rename(lightness = colour_lightness)

# add to trait data frame
d.traits.m.agg <- d.traits.m.agg %>% 
  left_join(d.lightness, by = "Name_std")

# there are missing values. Estimate these based on congenerics / similar species:

lightness_missing <- d.traits.m.agg %>%
  filter(is.na(lightness)) %>% 
  select(Name_std) %>% 
  left_join(d.std %>% 
              select(Name_std, Genus) %>% 
              distinct(),
            by = "Name_std")

d.traits.m.agg_corr <- d.traits.m.agg
d.traits.m.agg_corr$estimation <- NA
for (i in 1:nrow(lightness_missing)){
  
  sp_i <- lightness_missing$Name_std[i] # species of observation
  genus_i <- lightness_missing$Genus[i] # genus of observation
  
  subset <- d.lightness %>% 
    left_join(d.std %>% select(Name_std, Genus) %>% distinct(), 
              by = "Name_std") %>% 
    filter(Genus == genus_i,
           Name_std != sp_i,
           !is.na(lightness))
  
  
  # if there are no other species of the genus, take similar species
  if (nrow(subset) == 0){
    subset <- d.lightness %>%
      filter(Name_std != sp_i,
             Name_std %in% d.similar.species$Similar_Species[d.similar.species$Missing_SpeciesID == sp_i],
             !is.na(lightness))
  }
  
  if (nrow(subset) > 0){
    
    subset_agg <- subset %>% 
      summarise(lightness = mean(lightness, na.rm = T))
    
    d.traits.m.agg_corr$lightness[d.traits.m.agg_corr$Name_std == sp_i] <- subset_agg$lightness
    
    d.traits.m.agg_corr$estimation[d.traits.m.agg_corr$Name_std == sp_i] <- T
    
    d.estimates <- data.frame(Estimated_for = sp_i,
                              Estimated_traits = "lightness",
                              Estimated_from = unique(subset$Name_std)) %>% 
      bind_rows(d.estimates, .)
    
  } else {
    warning(paste(sp_i,"has no similar species."))
  }
  
}

d.traits.m.agg <- d.traits.m.agg_corr %>% 
  select(-estimation)

v.sel.traits <- c(v.sel.traits, "lightness")
v.sel.traits.ext <- paste0(v.sel.traits,  ".CWM")


################################################################################.
# EDIT OTHER DATA ------------------------- ####################################
################################################################################.

# aggregate abundance data #####################################################
################################################################################.

# Site scale -------------------------------------------------------------------.

d.abund.site.agg <- d.abund.ffb %>% 
  group_by(PlotID, Name_std) %>% 
  summarise(NumberAdults = sum(NumberAdults)) %>% 
  ungroup()  %>% 
  mutate(abund_rel = NumberAdults / sum(NumberAdults)) %>% 
  bind_rows(d.abund.ffk %>% 
              group_by(PlotID, Name_std) %>% 
              summarise(NumberAdults = sum(NumberAdults)) %>% 
              ungroup()  %>% 
              mutate(abund_rel = NumberAdults / sum(NumberAdults))) %>% 
  bind_rows(d.abund.bf %>% 
              group_by(PlotID, Name_std) %>% 
              summarise(NumberAdults = sum(NumberAdults)) %>% 
              ungroup()  %>% 
              mutate(abund_rel = NumberAdults / sum(NumberAdults))) %>% 
  bind_rows(d.abund.stem %>% 
              group_by(PlotID, Name_std) %>% 
              summarise(NumberAdults = sum(ToplaAdults)) %>% 
              ungroup()  %>% 
              mutate(abund_rel = NumberAdults / sum(NumberAdults))) %>% 
  group_by(PlotID, Name_std) %>% 
  summarise(NumberAdults = sum(abund_rel)) %>%  # NumberAdults for convenience reasons (same term as in other dfs)
  ungroup() %>% 
  mutate(ID = PlotID)


# Patch scale ------------------------------------------------------------------.

d.abund.stem.agg <- d.abund.stem %>% 
  group_by(PlotID, TreeID, Name_std) %>% 
  summarise(NumberAdults = sum(ToplaAdults),
            .groups = "drop") %>% 
  mutate(ID = paste(PlotID, TreeID, sep = "_"))

d.abund.stem.agg.2p <- d.abund.stem %>% 
  mutate(Period = ifelse(CollectionYear %in% 2010:2013,
                         "1st", "2nd")) %>% 
  group_by(PlotID, TreeID, Period, Name_std) %>% 
  summarise(NumberAdults = sum(ToplaAdults),
            .groups = "drop") %>% 
  mutate(ID = paste(PlotID, TreeID, Period, sep = "_"))

d.abund.patch.agg <- d.abund.stem %>% 
  group_by(PlotID, Name_std) %>% 
  summarise(NumberAdults = sum(ToplaAdults),
            .groups = "drop") %>% 
  mutate(ID = PlotID)


# NMDS tree species composition ################################################
################################################################################.

set.seed(21)
nmds.trees <- d.tree1 %>%
  bind_rows(d.tree2) %>% 
  filter(PlotID %in% d.abund.patch.agg$PlotID) %>% 
  group_by(PlotID, Species) %>%
  summarise(BA = mean(BA, na.rm = T)) %>%
  spread(Species, BA, fill = 0) %>%
  column_to_rownames("PlotID") %>%
  metaMDS(distance = "bray", k = 2, trymax = 400)

d.nmds.trees <- data.frame(scores(nmds.trees, display = c("sites", "species"),
                                  shrink = FALSE)) %>%
  rownames_to_column("PlotID")


################################################################################.
# ANALYSES ------------------------------------ ################################
################################################################################.

# descriptive stats ############################################################
################################################################################.

d.abund.stem %>% 
  rename(NumberAdults = ToplaAdults) %>% 
  select(-CollectionYear) %>% 
  bind_rows(d.abund.ffk, d.abund.ffb, d.abund.bf) %>% 
  summarise(n = sum(NumberAdults),
            n_spec = length(unique(Name_std)))

d.abund.stem %>% 
  rename(NumberAdults = ToplaAdults) %>% 
  select(-CollectionYear) %>% 
  bind_rows(d.abund.ffk, d.abund.ffb, d.abund.bf) %>% 
  group_by(PlotID) %>% 
  summarise(n = sum(NumberAdults),
            n_spec = length(unique(Name_std)),
            .groups = "drop") %>% 
  summarise(n_mean = mean(n),
            n_sd = sd(n),
            n_spec_mean = mean(n_spec),
            n_spec_sd = sd(n_spec))


d.abund.stem %>% 
  rename(NumberAdults = ToplaAdults) %>% 
  summarise(n = sum(NumberAdults),
            n_spec = length(unique(Name_std)))

d.abund.stem %>% 
  rename(NumberAdults = ToplaAdults) %>% 
  group_by(PlotID) %>% 
  summarise(n = sum(NumberAdults),
            n_spec = length(unique(Name_std)),
            .groups = "drop") %>% 
  summarise(n_mean = mean(n),
            n_sd = sd(n),
            n_spec_mean = mean(n_spec),
            n_spec_sd = sd(n_spec))

d.abund.stem.agg.2p %>% 
  group_by(TreeID) %>% 
  summarise(n = sum(NumberAdults),
            n_spec = length(unique(Name_std)),
            .groups = "drop") %>% 
  summarise(n_mean = mean(n),
            n_sd = sd(n),
            n_spec_mean = mean(n_spec),
            n_spec_sd = sd(n_spec))


(tot <- d.abund.stem %>% 
    mutate(type = "TEK") %>% 
    rename(NumberAdults = ToplaAdults) %>% 
    select(-CollectionYear) %>% 
    bind_rows(d.abund.ffk %>% mutate(type = "FFK"),
              d.abund.ffb %>% mutate(type = "FFB"), 
              d.abund.bf %>% mutate(type = "BF")) %>% 
    group_by(type) %>% 
    summarise(n = sum(NumberAdults),
              n_spec = length(unique(Name_std))))

d.abund.stem %>% 
  mutate(type = "TEK") %>% 
  rename(NumberAdults = ToplaAdults) %>% 
  select(-CollectionYear) %>% 
  bind_rows(d.abund.ffk %>% mutate(type = "FFK"),
            d.abund.ffb %>% mutate(type = "FFB"), 
            d.abund.bf %>% mutate(type = "BF"),
            d.abund.stem.agg %>% mutate(PlotID = TreeID, type = "TEK_log")) %>% 
  group_by(type, PlotID) %>% 
  summarise(n = sum(NumberAdults),
            n_spec = length(unique(Name_std)),
            .groups = "drop") %>% 
  group_by(type) %>% 
  summarise(n_sum = sum(n),
            n_mean = mean(n),
            n_sd = sd(n),
            n_spec_mean = mean(n_spec),
            n_spec_sd = sd(n_spec)) %>% 
  left_join(tot %>% 
              rename(n_spec_sum = n_spec) %>% 
              select(type, n_spec_sum), by = "type") %>% 
  select(type, n_sum, n_mean, n_sd, n_spec_sum, n_spec_mean, n_spec_sd)

# ranges of environmental variables --------------------------------------------.

d.lands %>% 
  filter(PlotID %in% d.abund.stem$PlotID) %>% 
  summarise(min = min(prop_forest2000),
            mean = mean(prop_forest2000),
            sd = sd(prop_forest2000),
            max = max(prop_forest2000))

d.deadwood %>% 
  filter(PlotID %in% d.abund.stem$PlotID) %>%
  summarise(min = min(DWV),
            mean = mean(DWV),
            sd = sd(DWV),
            max = max(DWV))

d.ssc %>% 
  filter(PlotID %in% d.abund.stem$PlotID) %>%
  summarise(min = min(ssc_mean),
            mean = mean(ssc_mean),
            max = max(ssc_mean))

d.canopy %>% 
  filter(PlotID %in% d.abund.stem$PlotID) %>%
  summarise(min = min(Canopy_closure),
            mean = mean(Canopy_closure),
            sd = sd(Canopy_closure),
            max = max(Canopy_closure))

d.canopy.sp %>% 
  filter(PlotID %in% d.abund.stem$PlotID) %>%
  summarise(min = min(Cover_above2m),
            mean = mean(Cover_above2m),
            sd = sd(Cover_above2m),
            max = max(Cover_above2m))

# set basic functions & parameters #############################################
################################################################################.

# Define dbFD.edit function, which creates customized FD and CWM output.

dbFD.edit <- function(traits, d.traits, d.abund, calc.FDRao = F, ...){
  
  ids <- sort(unique(d.abund$ID))
  
  d.traits <- d.traits %>% 
    select(Name_std, one_of(traits))
  
  for (trait_i in traits){
    d.traits <- d.traits %>% 
      filter(!is.na(!! sym(trait_i)))
  }
  
  d.traits <- d.traits %>% 
    column_to_rownames("Name_std")
  
  d.abund <-
    d.abund %>% 
    filter(Name_std %in% rownames(d.traits)) %>% 
    select(ID, Name_std, NumberAdults) %>% 
    spread(Name_std, NumberAdults, fill = 0) %>%
    column_to_rownames("ID")
  
  d.traits <- d.traits[names(d.abund), , drop = F]
  
  d.abund <- d.abund[, rownames(d.traits), drop = F]
  
  
  if (calc.FDRao){
    temp <- dbFD(d.traits, d.abund, ...)
    out <- data.frame(ID = rownames(temp$CWM),
                      temp$CWM,
                      FDRao = temp$RaoQ,
                      stringsAsFactors = F) %>% 
      right_join(data.frame(ID = ids,
                            stringsAsFactors = F), by = "ID")
  } else {
    temp <- functcomp(d.traits, as.matrix(d.abund))
    out <- data.frame(ID = rownames(temp),
                      temp,
                      stringsAsFactors = F) %>% 
      right_join(data.frame(ID = ids,
                            stringsAsFactors = F), by = "ID")
  }
  
  
  
  out
  
}

# Define a vector containing all the traits that cannot contain NAs per default

v.sel.traits.noNA <- v.sel.traits[v.sel.traits != "wing.aspect"]

# compute CWMs and FD ##########################################################
################################################################################.

# Site level -------------------------------------------------------------------.

d.rao.site <- dbFD.edit(v.sel.traits.noNA, d.traits.m.agg, d.abund.site.agg, 
                        calc.FDRao = T,
                        calc.FRic = F, calc.FGR = F, calc.FDiv = F)

cl <- makeCluster(future::availableCores())
clusterExport(cl, c("d.abund.site.agg", "dbFD.edit", "v.sel.traits", "d.traits.m.agg"))
clusterEvalQ(cl, {library(tidyverse); library(FD)})
l.cwm.site <- parLapply(cl, v.sel.traits, dbFD.edit, d.traits.m.agg, d.abund.site.agg,
                        calc.FDRao = T,
                        calc.FRic = F, calc.FGR = F, calc.FDiv = F)
stopCluster(cl)


d.cwm.site <- do.call(cbind, lapply(l.cwm.site, 
                                    function(x)  {
                                      names(x)[3] <- paste(names(x)[2], 
                                                           "FD", sep = ".")
                                      names(x)[2] <- paste(names(x)[2], 
                                                           "CWM", sep = ".")
                                      x <- column_to_rownames(x, "ID")
                                      x
                                    })) %>%
  rownames_to_column("ID") %>% 
  left_join(select(d.rao.site, ID, FDRao), by = "ID")

d.cwm.site <- d.cwm.site %>% 
  mutate(PlotID = ID)

# Patch level ------------------------------------------------------------------.

d.rao.patch <- dbFD.edit(v.sel.traits.noNA, d.traits.m.agg, d.abund.patch.agg, 
                         calc.FDRao = T,
                         calc.FRic = F, calc.FGR = F, calc.FDiv = F)

cl <- makeCluster(future::availableCores())
clusterExport(cl, c("d.abund.patch.agg", "dbFD.edit", "v.sel.traits", "d.traits.m.agg"))
clusterEvalQ(cl, {library(tidyverse); library(FD)})
l.cwm.patch <- parLapply(cl, v.sel.traits, dbFD.edit, d.traits.m.agg, d.abund.patch.agg,
                         calc.FDRao = T,
                         calc.FRic = F, calc.FGR = F, calc.FDiv = F)
stopCluster(cl)

d.cwm.patch <- do.call(cbind, lapply(l.cwm.patch, 
                                     function(x)  {
                                       names(x)[3] <- paste(names(x)[2], 
                                                            "FD", sep = ".")
                                       names(x)[2] <- paste(names(x)[2], 
                                                            "CWM", sep = ".")
                                       x <- column_to_rownames(x, "ID")
                                       x
                                     })) %>%
  rownames_to_column("ID") %>% 
  left_join(select(d.rao.patch, ID, FDRao), by = "ID")


d.cwm.patch <- d.cwm.patch %>% 
  mutate(PlotID = ID) %>% 
  select(-ID)

# Object level -----------------------------------------------------------------.

d.rao.object <- dbFD.edit(v.sel.traits.noNA, d.traits.m.agg, d.abund.stem.agg.2p, 
                          calc.FDRao = T,
                          calc.FRic = F, calc.FGR = F, calc.FDiv = F)

cl <- makeCluster(future::availableCores())
clusterExport(cl, c("d.abund.stem.agg.2p", "dbFD.edit", "v.sel.traits", "d.traits.m.agg"))
clusterEvalQ(cl, {library(tidyverse); library(FD)})
l.cwm.object <- parLapply(cl, v.sel.traits, dbFD.edit, d.traits.m.agg, d.abund.stem.agg.2p,
                          calc.FDRao = T,
                          calc.FRic = F, calc.FGR = F, calc.FDiv = F)
stopCluster(cl)

d.cwm.object <- do.call(cbind, lapply(l.cwm.object, 
                                      function(x)  {
                                        names(x)[3] <- paste(names(x)[2], 
                                                             "FD", sep = ".")
                                        names(x)[2] <- paste(names(x)[2], 
                                                             "CWM", sep = ".")
                                        x <- column_to_rownames(x, "ID")
                                        x
                                      })) %>%
  rownames_to_column("ID") %>% 
  left_join(select(d.rao.object, ID, FDRao), by = "ID")


d.cwm.object <- d.cwm.object %>% 
  mutate(PlotID = substr(ID, 1, 5),
         Period = substr(ID, nchar(ID) - 2, nchar(ID)),
         TreeID = substr(ID, 7, nchar(ID) - 4),
         Tree_sp = gsub("[!0-9]", "", TreeID))


# nullmodels: Region to site ###################################################
################################################################################.

# Define function f.null.site, which will create one null model run for each site 
# (i.e., sample species from the regional species pool, stratified by regional species abundances) and calculate FD and CWM for that null model draw.

f.null.site <- function(run_i){
  
  out <- list()
  for (pl_i in unique(d.abund.site.agg$PlotID)){
    
    d.abund.regionpool <- d.abund.site.agg %>% 
      filter(substr(PlotID, 1, 3) == substr(pl_i, 1, 3)) %>% 
      group_by(Name_std) %>% 
      summarise(NumberAdults = sum(NumberAdults)) %>% 
      ungroup()
    
    d.abund.target <- d.abund.site.agg %>% 
      filter(PlotID == pl_i) %>% 
      arrange(-NumberAdults)
    
    # shuffle species names in abundance data frame. Start with most abundant species
    # probability of drawing a species corresponds to species abundance in the region
    d.abund.null <- d.abund.target %>% 
      mutate(Name_std = sample(d.abund.regionpool$Name_std, nrow(.), replace = F, 
                               prob = d.abund.regionpool$NumberAdults))
    
    
    d.rao.site.null <- dbFD.edit(v.sel.traits.noNA, d.traits.m.agg, d.abund.null, 
                                 calc.FDRao = T,
                                 calc.FRic = F, calc.FGR = F, calc.FDiv = F)
    
    l.cwm.site.null <- lapply(v.sel.traits, dbFD.edit, d.traits.m.agg, d.abund.null,
                              calc.FDRao = T,
                              calc.FRic = F, calc.FGR = F, calc.FDiv = F)
    
    out[[pl_i]] <- do.call(cbind, lapply(l.cwm.site.null, 
                                         function(x)  {
                                           names(x)[3] <- paste(names(x)[2], 
                                                                "FD", sep = ".")
                                           names(x)[2] <- paste(names(x)[2], 
                                                                "CWM", sep = ".")
                                           x <- column_to_rownames(x, "ID")
                                           x
                                         })) %>%
      rownames_to_column("ID") %>% 
      left_join(select(d.rao.site.null, ID, FDRao), by = "ID")
    
    out[[pl_i]] <- out[[pl_i]] %>% 
      mutate(PlotID = pl_i,
             run = run_i)
  }
  
  do.call(rbind, out)
  
}

# Apply the function (paralleling):

set.seed(35)
cl <- makeCluster(future::availableCores())
clusterExport(cl, c("d.abund.site.agg", "dbFD.edit", "v.sel.traits.noNA", "v.sel.traits", "d.traits.m.agg"))
clusterEvalQ(cl, {library(tidyverse); library(FD)})
d.cwm.site.null <- parLapply(cl, 1:1999, f.null.site)
d.cwm.site.null <- do.call(rbind, d.cwm.site.null)
stopCluster(cl)

# nullmodels Site to patch #####################################################
################################################################################.

# Define function f.null.patch, which will create one null model run for each patch
# (i.e., sample species from site pool, stratified by species abundances) and calculate FD and CWM for that null model draw.

f.null.patch <- function(run_i){
  
  out <- list()
  
  for (pl_i in unique(d.abund.site.agg$PlotID)){
    
    d.abund.sitepool <- d.abund.site.agg %>% 
      filter(PlotID == pl_i) 
    
    d.abund.sitepool <- d.abund.sitepool %>% 
      bind_rows(d.abund.patch.agg %>% 
                  filter(PlotID == pl_i,
                         !Name_std %in% d.abund.sitepool$Name_std) %>% 
                  mutate(NumberAdults = min(d.abund.sitepool$NumberAdults)))
    
    d.abund.target <- d.abund.patch.agg %>% 
      filter(PlotID == pl_i) %>% 
      arrange(-NumberAdults)
    
    # shuffle species names in abundance data frame. Start with most abundant species
    # probability of drawing a species corresponds to species abundance at the site
    d.abund.null <- d.abund.target %>% 
      mutate(Name_std = sample(d.abund.sitepool$Name_std, nrow(.), replace = F, 
                               prob = d.abund.sitepool$NumberAdults))
    
    
    d.rao.patch.null <- dbFD.edit(v.sel.traits.noNA, d.traits.m.agg, d.abund.null, 
                                  calc.FDRao = T,
                                  calc.FRic = F, calc.FGR = F, calc.FDiv = F)
    
    l.cwm.patch.null <- lapply(v.sel.traits, dbFD.edit, d.traits.m.agg, d.abund.null,
                               calc.FDRao = T,
                               calc.FRic = F, calc.FGR = F, calc.FDiv = F)
    
    out[[pl_i]] <- do.call(cbind, lapply(l.cwm.patch.null, 
                                         function(x)  {
                                           names(x)[3] <- paste(names(x)[2], 
                                                                "FD", sep = ".")
                                           names(x)[2] <- paste(names(x)[2], 
                                                                "CWM", sep = ".")
                                           x <- column_to_rownames(x, "ID")
                                           x
                                         })) %>%
      rownames_to_column("ID") %>% 
      left_join(select(d.rao.patch.null, ID, FDRao), by = "ID")
    
    out[[pl_i]] <- out[[pl_i]] %>% 
      mutate(PlotID = ID,
             run = run_i)
  }
  
  do.call(rbind, out) 
}

# Run this function (parallelized):

set.seed(38)
cl <- makeCluster(future::availableCores())
clusterExport(cl, c("d.abund.site.agg", "d.abund.patch.agg", "dbFD.edit", 
                    "v.sel.traits.noNA", "v.sel.traits", "d.traits.m.agg"))
clusterEvalQ(cl, {library(tidyverse); library(FD)})
d.cwm.patch.null <- parLapply(cl, 1:1999, f.null.patch)
d.cwm.patch.null <- do.call(rbind, d.cwm.patch.null)
stopCluster(cl)

# nullmodels: Patch to object ##################################################
################################################################################.

# Define function f.null.object, which will create one null model run for each object 
# (i.e., sample species from patch pool, stratified by species abundances) and calculate FD and CWM for that null model draw.

f.null.object <- function(run_i){
  
  out <- list()
  
  
  for (tr_i in unique(d.abund.stem.agg.2p$TreeID)){
    for (p_i in c("1st", "2nd")){
      
      pl_i <- d.abund.stem.agg.2p %>% 
        filter(TreeID == tr_i) %>% 
        {unique(.$PlotID)}
      
      d.abund.patchpool <- d.abund.patch.agg %>% 
        filter(PlotID == pl_i) 
      
      d.abund.target <- d.abund.stem.agg.2p %>% 
        filter(TreeID == tr_i,
               Period == p_i) %>% 
        arrange(-NumberAdults)
      
      cont = T
      while (cont){ 
        # shuffle species names in abundance data frame. Start with most abundant species
        # probability of drawing a species corresponds to species abundance on the patch
        d.abund.null <- d.abund.target %>% 
          mutate(Name_std = sample(d.abund.patchpool$Name_std, nrow(.), replace = F, 
                                   prob = d.abund.patchpool$NumberAdults))
        
        
        d.rao.object.null <- tryCatch({dbFD.edit(v.sel.traits.noNA, d.traits.m.agg, d.abund.null, 
                                                 calc.FDRao = T,
                                                 calc.FRic = F, calc.FGR = F, calc.FDiv = F)},
                                      error = function(e) "error")
        
        l.cwm.object.null <- tryCatch({lapply(v.sel.traits, dbFD.edit, d.traits.m.agg, d.abund.null,
                                              calc.FDRao = T,
                                              calc.FRic = F, calc.FGR = F, calc.FDiv = F)},
                                      error = function(e) "error")
        
        if (any(c(grepl("error", d.rao.object.null), grepl("error", l.cwm.object.null)))){
          cont <- T
        } else {
          cont <- F
        }
      }
      
      out[[paste(tr_i, p_i, sep = ".")]] <- do.call(cbind, lapply(l.cwm.object.null, 
                                                                  function(x)  {
                                                                    names(x)[3] <- paste(names(x)[2], 
                                                                                         "FD", sep = ".")
                                                                    names(x)[2] <- paste(names(x)[2], 
                                                                                         "CWM", sep = ".")
                                                                    x <- column_to_rownames(x, "ID")
                                                                    x
                                                                  })) %>%
        rownames_to_column("ID") %>% 
        left_join(select(d.rao.object.null, ID, FDRao), by = "ID")
      
      out[[paste(tr_i, p_i, sep = ".")]]  <- out[[paste(tr_i, p_i, sep = ".")]]  %>% 
        mutate(PlotID = pl_i,
               TreeID = tr_i,
               Period = p_i,
               run = run_i)
    }
  }
  
  do.call(rbind, out)
  
}

# Run the function (parallelized):

set.seed(98)
cl <- makeCluster(future::availableCores())
clusterExport(cl, c("d.abund.stem.agg.2p", "d.abund.patch.agg", "dbFD.edit", "v.sel.traits", "v.sel.traits.noNA",
                    "d.traits.m.agg"))
clusterEvalQ(cl, {library(tidyverse); library(FD)})
d.cwm.object.null <- parLapply(cl, 1:1999, f.null.object)
d.cwm.object.null <- do.call(rbind, d.cwm.object.null)
stopCluster(cl)


# prepare data for final analyses ##############################################
################################################################################.

# Site level -------------------------------------------------------------------.

d.cwm.site.nulldiff <- d.cwm.site.null %>%
  select(-c(ID, run)) %>% 
  group_by(PlotID) %>% 
  summarise_all(.funs = ~mean(.)) %>% 
  ungroup() %>% 
  left_join(select(d.cwm.site, one_of(names(d.cwm.site.null)), -ID), 
            by = c("PlotID"), suffix = c(".null", "")) %>% 
  gather(what, value, -c(PlotID)) %>% 
  mutate(kind = ifelse(grepl(".null", what), 
                       "null",
                       "obs"),
         variable = gsub(".null", "", what)) %>% 
  select(-what) %>% 
  spread(kind, value) %>% 
  mutate(diff = obs - null) %>% 
  select(-c(null, obs)) %>% 
  spread(variable, diff) %>% 
  left_join(d.EP, by = "PlotID") %>%
  left_join(d.nmds.trees %>% 
              filter(PlotID %in% d.cwm.site.null$PlotID) %>% 
              mutate_if(~is.numeric(.), ~(. - mean(.)) / (2 * sd(.))), by = "PlotID") %>%
  left_join(d.ssc %>% 
              filter(PlotID %in% d.cwm.site.null$PlotID) %>% 
              mutate_if(~is.numeric(.), ~(. - mean(.)) / (2 * sd(.))), by = "PlotID") %>%
  left_join(d.canopy %>% 
              filter(PlotID %in% d.cwm.site.null$PlotID) %>% 
              mutate_if(~is.numeric(.), ~(. - mean(.)) / (2 * sd(.))), by = "PlotID") %>%
  left_join(d.deadwood %>% 
              filter(PlotID %in% d.cwm.site.null$PlotID) %>% 
              mutate_if(~is.numeric(.), ~(. - mean(.)) / (2 * sd(.))), by = "PlotID") %>% 
  left_join(d.lands %>% 
              filter(PlotID %in% d.cwm.site.null$PlotID) %>% 
              mutate_if(~is.numeric(.), ~(. - mean(.)) / (2 * sd(.))), by = "PlotID") %>% 
  select(PlotID, Exploratory, 
         !! enquo(v.sel.traits.ext), FDRao,
         prop_forest2000, NMDS1, NMDS2, ssc_mean, Canopy_closure, DWV)

# Patch level ------------------------------------------------------------------.

d.cwm.patch.nulldiff <- d.cwm.patch.null %>% 
  select(-c(ID, run)) %>% 
  group_by(PlotID) %>% 
  summarise_all(.funs = ~mean(.)) %>% 
  ungroup() %>% 
  left_join(select(d.cwm.patch, one_of(names(d.cwm.patch.null))), 
            by = c("PlotID"), suffix = c(".null", "")) %>% 
  gather(what, value, -c(PlotID)) %>% 
  mutate(kind = ifelse(grepl(".null", what), 
                       "null",
                       "obs"),
         variable = gsub(".null", "", what)) %>% 
  select(-what) %>% 
  spread(kind, value) %>% 
  mutate(diff = obs - null) %>% 
  select(-c(null, obs)) %>% 
  spread(variable, diff) %>% 
  left_join(d.EP, by = "PlotID") %>%
  left_join(d.nmds.trees %>% 
              filter(PlotID %in% d.cwm.patch.null$PlotID) %>% 
              mutate_if(~is.numeric(.), ~(. - mean(.)) / (2 * sd(.))), by = "PlotID") %>%
  left_join(d.ssc %>% 
              filter(PlotID %in% d.cwm.patch.null$PlotID) %>% 
              mutate_if(~is.numeric(.), ~(. - mean(.)) / (2 * sd(.))), by = "PlotID") %>%
  left_join(d.canopy.sp %>% 
              filter(PlotID %in% d.cwm.patch.null$PlotID) %>% 
              mutate_if(~is.numeric(.), ~(. - mean(.)) / (2 * sd(.))), by = "PlotID") %>% 
  left_join(d.deadwood %>% 
              filter(PlotID %in% d.cwm.patch.null$PlotID) %>% 
              mutate_if(~is.numeric(.), ~(. - mean(.)) / (2 * sd(.))), by = "PlotID") %>%
  select(PlotID, Exploratory,
         !! enquo(v.sel.traits.ext), FDRao,
         ssc_mean, Cover_above2m, DWV)

# Object level -----------------------------------------------------------------.

d.cwm.object.nulldiff <- d.cwm.object.null %>% 
  select(-c(ID, run)) %>% 
  group_by(PlotID, TreeID, Period) %>% 
  summarise_all(.funs = ~mean(.)) %>% 
  ungroup() %>% 
  left_join(select(d.cwm.object, one_of(names(d.cwm.object.null))), 
            by = c("PlotID", "TreeID", "Period"), suffix = c(".null", "")) %>% 
  select(-ID) %>% 
  gather(what, value, -c(PlotID, TreeID, Period)) %>% 
  mutate(kind = ifelse(grepl(".null", what), 
                       "null",
                       "obs"),
         variable = gsub(".null", "", what)) %>% 
  select(-what) %>% 
  spread(kind, value) %>% 
  mutate(diff = obs - null) %>% 
  select(-c(null, obs)) %>% 
  spread(variable, diff) %>% 
  mutate(Tree_sp = gsub("[!0-9]", "", TreeID),
         Tree_sp = factor(Tree_sp, levels = names(tree.order))) 

# main plots ###################################################################
################################################################################.

d.nullmeans.site <- d.cwm.site.null %>% 
  select(-c(ID, PlotID, run)) %>% 
  summarise_all(list(mn = ~mean(.),
                     sd = ~sd(.)))

d.nullmeans.patch <- d.cwm.patch.null %>% 
  select(-c(ID, PlotID, run)) %>% 
  summarise_all(list(mn = ~mean(.),
                     sd = ~sd(.)))

d.nullmeans.object <- d.cwm.object.null %>% 
  select(-c(ID, PlotID, TreeID, Period, run)) %>% 
  summarise_all(list(mn = ~mean(.),
                     sd = ~sd(.)))


# define function that prepares data for plotting for each trait
f_levelplot <- function(trait_i){
  
  # null models ################################################################.
  
  # compute grand means for the three levels
  mns_i <- data.frame(mean = c(d.nullmeans.site[, paste0(trait_i, "_mn")],
                               d.nullmeans.patch[, paste0(trait_i, "_mn")],
                               d.nullmeans.object[, paste0(trait_i, "_mn")]),
                      level = c("Site", "Patch", "Object"))
  
  # determine density site-level null models
  density.null.site <- data.frame(density(d.cwm.site.null[, trait_i])[c("x", "y")])
  
  # highest definition interval (only for subsetting)
  range_hdi_site <- bayestestR::hdi(d.cwm.site.null[, trait_i], ci = .95)
  
  # subset based on highest definition interval
  density.null.site <- density.null.site %>%
    filter(x >= range_hdi_site$CI_low,
           x <= range_hdi_site$CI_high) %>%
    mutate(y = y / max(y), # scale to max 1
           x_diff = median(diff(x))) # add differences for later plotting (width of rects)
  
  # determine density patch-level null models
  density.null.patch <- data.frame(density(d.cwm.patch.null[, trait_i])[c("x", "y")])
  
  # highest definition interval (only for subsetting)
  range_hdi_patch <- bayestestR::hdi(d.cwm.patch.null[, trait_i], ci = .95)
  
  # subset based on highest definition interval
  density.null.patch <- density.null.patch %>%
    filter(x >= range_hdi_patch$CI_low,
           x <= range_hdi_patch$CI_high) %>%
    mutate(y = y / max(y), # scale to max 1
           x_diff = median(diff(x))) # add differences for later plotting (width of rects)
  
  
  # determine density object-level null models
  density.null.object <- data.frame(density(d.cwm.object.null[, trait_i])[c("x", "y")])
  
  # highest definition interval (only for subsetting)
  range_hdi_object <- bayestestR::hdi(d.cwm.object.null[, trait_i], ci = .95)
  
  # subset based on highest definition interval
  density.null.object <- density.null.object %>% 
    filter(x >= range_hdi_object$CI_low,
           x <= range_hdi_object$CI_high) %>% 
    mutate(y = y / max(y), # scale to max 1
           x_diff = median(diff(x))) # add differences for later plotting (width of rects)
  
  
  # observations ###############################################################.
  
  # Model site-level -----------------------------------------------------------.
  
  # subset for the respective trait
  d.target.site <- d.cwm.site.nulldiff %>% 
    filter(!is.na(!! sym(trait_i)))
  
  formula_site <- 
    as.formula(paste(trait_i, " ~ NMDS1 + NMDS2 + prop_forest2000 + DWV +
          Canopy_closure + ssc_mean"))
  
  # linear model
  mod.site <- lm(formula_site, data = d.target.site, na.action = na.fail)
  
  # Model patch-level ----------------------------------------------------------.
  
  # subset for the respective trait
  d.target.patch <- d.cwm.patch.nulldiff %>% 
    filter(!is.na(!! sym(trait_i)))
  
  formula_patch <- 
    as.formula(paste(trait_i, " ~ DWV + Cover_above2m + ssc_mean"))
  
  # linear model
  mod.patch <- lm(formula_patch, data = d.target.patch, na.action = na.fail)
  
  # Model object-level ---------------------------------------------------------.
  
  # subset for the respective trait
  d.target.object <- d.cwm.object.nulldiff %>% 
    filter(!is.na(!! sym(trait_i)))
  
  formula_object <- 
    as.formula(paste(trait_i, " ~ Period * Tree_sp  + (1 | PlotID/TreeID)"))
  
  # change contrasts to sum (necessary for later plotting)
  options(contrasts = c("contr.sum", "contr.poly"))
  
  # linear mixed effects model
  mod.object <- glmmTMB(formula_object,
                        na.action = na.fail, 
                        data = d.target.object)
  
  pval.object <- car::Anova(mod.object, type = "III") %>%
    as.data.frame %>%
    rownames_to_column("var") %>%
    rename(p_val = `Pr(>Chisq)`) %>%
    mutate(sign = ifelse(p_val <= 0.05, "sign", "ns"))
  
  # extract coefficients
  d.plotting <- summary(mod.site)$coefficients %>% 
    as.data.frame() %>% 
    rownames_to_column("var") %>% 
    mutate(level = "Site") %>% 
    bind_rows(summary(mod.patch)$coefficients %>% 
                as.data.frame() %>% 
                rownames_to_column("var") %>% 
                mutate(level = "Patch")) %>% 
    rename(p_val = "Pr(>|t|)",
           se = "Std. Error") %>% 
    bind_rows(summary(mod.object)$coefficients$cond %>% 
                as.data.frame() %>% 
                rownames_to_column("var") %>% 
                mutate(level = "Object") %>% 
                rename(p_val = "Pr(>|z|)",
                       se = "Std. Error"))
  
  # change tree numbers (from sum contrasts) to names
  for (i in 1:nlevels(d.cwm.object.nulldiff$Tree_sp)){
    d.plotting <- d.plotting %>% 
      mutate(var = ifelse(var %in% c(paste0("Tree_sp", i), paste0("Period1:Tree_sp", i)),
                          gsub(paste0("Tree_sp", i), 
                               tree.order[levels(d.cwm.object.nulldiff$Tree_sp)[i]], 
                               var),
                          var))
  }
  
  # extract intercepts
  intercept_i <- d.plotting %>% 
    filter(var == "(Intercept)") %>% 
    select(Estimate, level) %>% 
    rename(intercept = Estimate)
  
  
  # add values for the one missing tree species (Pseudotsuga)
  d.plotting$lastlevel <- NA    
  
  # main effect
  d.plotting <- d.plotting %>% 
    filter(var %in% tree.order) %>% 
    summarise(Estimate = -sum(Estimate),
              se = mean(se), # approximation
              var = tree.order[length(tree.order)],
              level = "Object",
              lastlevel = T) %>% 
    mutate(p_val = ifelse(Estimate - 1.96 * se > 0 |
                            Estimate + 1.96 * se < 0,
                          0.05, 1)) %>% # p value approximation (confidence intervals crossing zero?)
    bind_rows(d.plotting, .) 
  
  # interaction
  d.plotting <- d.plotting %>% 
    filter(var %in% paste0("Period1:",tree.order)) %>% 
    summarise(Estimate = -sum(Estimate),
              se = mean(se), # approximation
              var = paste0("Period1:", tree.order[length(tree.order)]),
              level = "Object",
              lastlevel = T) %>% 
    mutate(p_val = ifelse(Estimate - 1.96 * se > 0 |
                            Estimate + 1.96 * se < 0,
                          0.05, 1)) %>% # p value approximation (confidence intervals crossing zero?)
    bind_rows(d.plotting, .)
  
  
  
  # last edits
  d.plotting <-
    d.plotting %>% 
    filter(var != "(Intercept)") %>% # not used for plotting
    left_join(mns_i, by = "level") %>% # add level means
    left_join(intercept_i, by = "level") %>% # add level intercepts
    mutate(Tree_sp = ifelse(grepl(paste(tree.order, collapse = "|"), var),
                            gsub("Period1:", "", var), NA), # add a tree species variable
           period1_plus = Estimate[var == "Period1"]) %>% # add estimate of Period1 as an extra variable
    filter(var != "Period1") %>% # exclude the Period1 estimate (will be included in tree-species effects)
    group_by(Tree_sp) %>%  # determine estimates (predictions) per tree species for first and second period
    mutate(treesp_plus = unique(Estimate[!grepl(":", var)]), # tree species main effect estimate
           int_plus = ifelse(is.na(Tree_sp), NA, unique(Estimate[grepl(":", var)])), # tree species x period interaction estimate
           Estimate = ifelse(!grepl(":", var), # calculate estimates for the different cases
                             ifelse(is.na(Tree_sp),
                                    mean + intercept + Estimate,                              
                                    mean + intercept + Estimate - period1_plus - int_plus),
                             mean + intercept + Estimate + period1_plus + treesp_plus)) %>% 
    ungroup() %>% 
    mutate(var = ifelse(!is.na(Tree_sp),
                        ifelse(grepl("Period1", var), "Period1", "Period2"), # change variable names for all tree species to Period1 and Period2
                        var),
           level = factor(level, levels = c("Site", "Patch", "Object")),
           sign = ifelse(p_val <= .05, "sign", "ns"), # significant effect?
           sign = ifelse(is.na(sign), "ns", sign),
           varlevel = interaction(var, level), # varlevel will be the x-axis of the plot
           varlevel = factor(varlevel, levels = c(paste0(vars.site, ".Site"),
                                                  paste0(vars.patch, ".Patch"),
                                                  paste0(c("Period1", "Period2"), ".Object"))))
  
  
  
  # prepare output
  out <- list(d.plotting = d.plotting, # main data frame for plotting
              density.null.site =  density.null.site, # null model density site level
              density.null.patch = density.null.patch, # null model density patch level
              density.null.object = density.null.object, # null model density object level
              intercept_mns = intercept_i %>% left_join(mns_i, by = "level"), # means and intercepts
              pval.object = pval.object) # p values for object-level model
  
  # add trait column to output
  lapply(out, function(x) bind_cols(x, trait = trait_i))
}

# define variable labels for plotting
var.labels <- c(prop_forest2000 = "Forest surroundings",
                DWV = "Deadwood volume",
                NMDS1 = "Broadleaf share",
                NMDS2 = "Conifer type",
                ssc_mean = "Structural complexity",
                Canopy_closure = "Canopy closure site",
                Cover_above2m = "Canopy closure patch",
                Period1 = "1st period",
                Period2 = "2nd period")



# FDRao plot -------------------------------------------------------------------.

# compile data using the function
l.raoplot <- f_levelplot("FDRao")

# define x-axis labels
xlabels <- var.labels[gsub(".Site|.Patch|.Object", "", 
                           levels(l.raoplot$d.plotting$varlevel))]
names(xlabels) <- levels(l.raoplot$d.plotting$varlevel)

# plot
p <- l.raoplot$d.plotting %>%
  ggplot() +
  geom_point(aes(x = varlevel, y = Estimate), col = NA)  + # dummy for numeric scaling of x
  geom_rect(data = l.raoplot$density.null.site, # null model shading site
            aes(ymin = x+x_diff/2, ymax = x-x_diff/2,
                xmin = -Inf, xmax = length(vars.site) + .5, alpha = y),
            fill = "grey60", col = NA) +
  geom_segment(data = l.raoplot$intercept_mns[l.raoplot$intercept_mns$level == "Site", ], # intercept site
               aes(x = -Inf,     
                   xend = length(vars.site)  + .5,
                   y = mean + intercept,
                   yend = mean + intercept,
                   col = "Site"), lty = 1) +
  geom_rect(data = l.raoplot$density.null.patch, # null model shading patch
            aes(ymin = x+x_diff/2, ymax = x-x_diff/2,
                xmin = length(vars.site) + .5,
                xmax = length(vars.site) + length(vars.patch) + .5,
                alpha = y),
            fill = "grey60", col = NA) +
  geom_segment(data = l.raoplot$intercept_mns[l.raoplot$intercept_mns$level == "Patch", ],  # intercept patch
               aes(x = length(vars.site) + .5,
                   xend = length(vars.site) + length(vars.patch) + .5,
                   y = mean + intercept,
                   yend = mean + intercept,
                   col = "Patch"), lty = 1) +
  geom_rect(data = l.raoplot$density.null.object, # null model shading object
            aes(ymin = x+x_diff/2, ymax = x-x_diff/2,
                xmin = length(vars.site) + length(vars.patch) + .5,
                xmax = Inf,
                alpha = y),
            fill = "grey60", col = NA) +
  geom_segment(data = l.raoplot$intercept_mns[l.raoplot$intercept_mns$level == "Object", ],  # intercept object
               aes(x = length(vars.site) + length(vars.patch) + .5,
                   xend = Inf,
                   y = mean + intercept,
                   yend = mean + intercept), col = "#7570B3", lty = 1) +
  geom_segment(data = . %>% filter(is.na(Tree_sp)), # confidence intervals for plot and patch effects
               aes(x = varlevel, y = Estimate, col = level,
                   yend = Estimate - 1.96 * se, xend = varlevel,
                   size = sign)) +
  geom_segment(data = . %>% filter(is.na(Tree_sp)), # confidence intervals for plot and patch effects
               aes(x = varlevel, y = Estimate, col = level,
                   yend = Estimate + 1.96 * se, xend = varlevel,
                   size = sign)) +
  geom_point(aes(x = varlevel, y = Estimate, shape = level), # white background for points
             col = "white", size = 2.75) +
  geom_point(aes(x = varlevel, y = Estimate, col = level, shape = level), # actual points
             size = 2) +
  geom_line(data = . %>% filter(!is.na(Tree_sp)), # connect points of same species
            aes(x = varlevel, y = Estimate, group = Tree_sp),
            alpha = .5, col = "#7570b3") +
  geom_vline(xintercept = seq(.5, 8.5, 1), lty = 3, size = .25) + # vertical separation lines
  geom_vline(xintercept = c(10.5, 11.5), lty = 3, size = .25) + # vertical separation lines
  geom_vline(xintercept = c(6.5, 9.5), lty = 1, size = .5) + # vertical separation lines
  scale_colour_manual(values = c(Site = "#1b9e77", Patch = "#d95f02", 
                                 Object = rgb(117/255,112/255,179/255, .5)),
                      breaks = c("Site", "Patch", "Object"), 
                      labels = c("Site", "Patch", "Object"),
                      name = "Filtering step") +
  scale_shape_manual(values = c(Site = 19, Patch = 15, Object = 17), 
                     labels = c("Site", "Patch", "Object"),
                     name = "Filtering step") +
  scale_size_manual(values = c(sign = 1.5, ns = .5), guide = F) +
  scale_x_discrete(labels = xlabels, drop = F) +
  scale_alpha_continuous(range = c(.1, .8), guide = F) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("Functional diversity (Rao's Q)") 

# extract y-axis range from plot (used to place the significance box)
yrange <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range

# add box for significances period/tree/interaction
p <- p + annotate(geom = "tile", x = 10.5, width = 1.5,
                  y = yrange[2] - diff(yrange) * .03, height = diff(yrange) * .1, fill = "grey80")

# add significance letters. either bold or transparent, depending on significanace
if (l.raoplot$pval.object[l.raoplot$pval.object$var == "Tree_sp", "sign"] == "sign") {
  p <- p + annotate(geom = "text", label = "T", fontface = "bold",
                    x = 10, y = yrange[2] - diff(yrange) * .03, hjust = .5, 
                    vjust = .5, col = "#7570b3", size = 6, alpha = 1)
} else {
  p <- p + annotate(geom = "text", label = "T",
                    x = 10, y = yrange[2] - diff(yrange) * .03, hjust = .5, 
                    vjust = .5, col = "#7570b3", size = 6, alpha = .5)
}
if (l.raoplot$pval.object[l.raoplot$pval.object$var == "Period:Tree_sp", "sign"] == "sign") {
  p <- p + annotate(geom = "text", label = "×", fontface = "bold",
                    x = 10.5, y = yrange[2] - diff(yrange) * .03, hjust = .5, 
                    vjust = .5, col = "#7570b3", size = 6, alpha = 1)
} else {
  p <- p + annotate(geom = "text", label = "×",
                    x = 10.5, y = yrange[2] - diff(yrange) * .03, hjust = .5, 
                    vjust = .5, col = "#7570b3", size = 6, alpha = .5)
}
if (l.raoplot$pval.object[l.raoplot$pval.object$var == "Period", "sign"] == "sign") {
  p <- p + annotate(geom = "text", label = "D", fontface = "bold",
                    x = 11, y = yrange[2] - diff(yrange) * .03, hjust = .5, 
                    vjust = .5, col = "#7570b3", size = 6, alpha = 1)
} else {
  p <- p + annotate(geom = "text", label = "D",
                    x = 11, y = yrange[2] - diff(yrange) * .03, hjust = .5, 
                    vjust = .5, col = "#7570b3", size = 6, alpha = .5)
}


p +
  theme(plot.margin = margin(0, 0, 0, 3, "mm"),
        legend.position = "top",
        legend.box.margin = margin(5, 0, 5, 0, "mm")) # to prevent that x-axis labels are cut

# CWM plot ---------------------------------------------------------------------.

l.combplot <- lapply(v.sel.traits.ext,
                     f_levelplot) # apply the function to all trait CWMs

# only three parts needed. convert to data frames
d.plotting <- do.call(rbind, lapply(l.combplot, function(x) x$d.plotting))
intercept_mns <- do.call(rbind, lapply(l.combplot, function(x) x$intercept_mns))
pval.object <- do.call(rbind, lapply(l.combplot, function(x) x$pval.object))

# define labels for traits
trait.labels <- c(wing.length.res.CWM = "Wing length",
                  wing.loading.rev.CWM = "Wing loading",
                  wing.aspect.CWM = "Wing aspect",
                  body.length.CWM = "Body length",
                  leg.length.res.CWM = "Leg length",
                  lightness.CWM = "Lightness",
                  antenna.length.res.CWM = "Antenna length",
                  eye.length.res.CWM = "Eye length",
                  hairiness.CWM = "Hairiness",
                  body.width.res.CWM = "Body width",
                  body.roundness.CWM = "Body roundness",
                  head.length.res.CWM = "Head length",
                  mandibular.aspect.CWM = "Mandibular aspect")

# extract y-axis range for each trait. This is needed to place the significance box
range_y <- d.plotting %>% 
  mutate(upper = ifelse(grepl("Period", var), Estimate, Estimate + 1.96 * se),
         lower = ifelse(grepl("Period", var), Estimate, Estimate - 1.96 * se)) %>% 
  group_by(trait) %>% 
  summarise(ymax = max(upper),
            ymin = min(lower),
            diff = diff(c(ymin, ymax)),
            .groups = "drop") %>% 
  mutate(trait = factor(trait, level = names(trait.labels)))

# edit p values of object-level variables
pval.object <-
  pval.object %>% 
  mutate(trait = factor(trait, level = names(trait.labels))) %>% # adjust factor levels (order!)
  filter(var != "(Intercept)") %>% 
  left_join(range_y, by = "trait") %>% # add y-range information (needed for plotting)
  mutate(Label = ifelse(var == "Period", "D", # Label that will be plotted
                        ifelse(var == "Tree_sp", "T", "×")),
         x_pos = ifelse(var == "Period", 11, # x-axis position at which will be plotted
                        ifelse(var == "Tree_sp", 10, 10.5)))


d.plotting <- d.plotting %>% 
  mutate(trait = factor(trait, level = names(trait.labels))) # adjust factor levels (order!)

intercept_mns <- intercept_mns %>% 
  mutate(trait = factor(trait, level = names(trait.labels))) # adjust factor levels (order!)


d.plotting %>%
  ggplot() +
  geom_point(aes(x = varlevel, y = Estimate), col = NA)  + # dummy points (to being able to add numeric stuff to x axis)
  geom_segment(data = intercept_mns[intercept_mns$level == "Site", ], # intercept site
               aes(x = -Inf,     
                   xend = length(vars.site)  + .5,
                   y = mean + intercept,
                   yend = mean + intercept), col = "#1b9e77", lty = 1) +
  geom_segment(data = intercept_mns[intercept_mns$level == "Patch", ],  # intercept patch
               aes(x = length(vars.site) + .5,
                   xend = length(vars.site) + length(vars.patch) + .5,
                   y = mean + intercept,
                   yend = mean + intercept), col = "#d95f02", lty = 1) +
  geom_segment(data = intercept_mns[intercept_mns$level == "Object", ], # intercept object
               aes(x = length(vars.site) + length(vars.patch) + .5,
                   xend = Inf,
                   y = mean + intercept,
                   yend = mean + intercept), col = "#7570B3", lty = 1) +
  geom_segment(data = . %>% filter(is.na(Tree_sp)), # confidence intercals for site and patch variables
               aes(x = varlevel, y = Estimate, col = level,
                   yend = Estimate - 1.96 * se, xend = varlevel,
                   size = sign)) +
  geom_segment(data = . %>% filter(is.na(Tree_sp)), # confidence intercals for site and patch variables
               aes(x = varlevel, y = Estimate, col = level,
                   yend = Estimate + 1.96 * se, xend = varlevel,
                   size = sign)) +
  geom_point(aes(x = varlevel, y = Estimate, shape = level), # white background of points
             col = "white", size = 2.75) +
  geom_point(aes(x = varlevel, y = Estimate, col = level, shape = level), # actual points
             size = 2) +
  geom_line(data = . %>% filter(!is.na(Tree_sp)), # connect two points per tree species
            aes(x = varlevel, y = Estimate, group = Tree_sp),
            alpha = .5, col = "#7570b3") +
  geom_vline(xintercept = seq(.5, 8.5, 1), lty = 3, size = .25) + # add vertical separation lines
  geom_vline(xintercept = c(10.5, 11.5), lty = 3, size = .25) + # add vertical separation lines
  geom_vline(xintercept = c(6.5, 9.5), lty = 1, size = .5) + # add vertical separation lines
  scale_colour_manual(values = c(Site = "#1b9e77", Patch = "#d95f02", 
                                 Object = rgb(117/255,112/255,179/255, .5)),
                      breaks = c("Site", "Patch", "Object"),
                      labels = c("Site                                   ", # space to being able to place the icons
                                 "Patch                                   ",  # space to being able to place the icons
                                 "Object                                 "), # space to being able to place the icons
                      name = "Filtering step") +
  scale_shape_manual(values = c(Site = 19, Patch = 15, Object = 17), 
                     labels = c("Site                                   ", # space to being able to place the icons
                                "Patch                                   ",  # space to being able to place the icons
                                "Object                                 "), # space to being able to place the icons
                     name = "Filtering step") +
  scale_size_manual(values = c(sign = 1.5, ns = .5), guide = F) + # significant lines bold
  scale_x_discrete(labels = xlabels, drop = F) +
  scale_alpha_continuous(range = c(.1, .8), guide = F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.x = element_blank(),
        legend.position = "top",
        legend.box.margin = margin(5, 0, 5, 0, "mm")) +
  facet_wrap(~ trait, scales = "free_y", ncol = 3, # wrap the facets per trait
             labeller = labeller(trait = trait.labels)) + # correct trait labels
  geom_tile(data = range_y, # add the significance box
            aes(x = 10.5, y = ymax + 0.1 * diff, height = diff * 0.12), # height & position relative to range (diff)
            width = 1.5, fill = "grey80") +
  geom_text(data = pval.object %>% filter(sign == "sign"), aes(x = x_pos, y = ymax + 0.1 * diff, # add significant letters (bold)
                                                               label = Label),
            hjust = .5, vjust = .5, fontface = "bold", col = rgb(117/255,112/255,179/255)) + 
  geom_text(data = pval.object %>% filter(sign == "ns"), aes(x = x_pos, y = ymax + 0.1 * diff, # add non-significant letters (transparent)
                                                             label = Label),
            hjust = .5, vjust = .5, col = rgb(117/255,112/255,179/255, .5)) +
  ylab("Predicted CWM")


# check whether intercepts are significant #####################################
################################################################################.

f_intercept <- function(trait_i){
  
  # Model site-level -----------------------------------------------------------.
  
  # subset for the respective trait
  d.target.site <- d.cwm.site.nulldiff %>% 
    filter(!is.na(!! sym(trait_i)))
  
  formula_site <- 
    as.formula(paste(trait_i, " ~ NMDS1 + NMDS2 + prop_forest2000 + DWV +
          Canopy_closure + ssc_mean"))
  
  # linear model
  mod.site <- lm(formula_site, data = d.target.site, na.action = na.fail)
  
  p_site <- summary(mod.site)$coefficients["(Intercept)", c("Estimate", "Pr(>|t|)")]
  
  # Model patch-level ----------------------------------------------------------.
  
  # subset for the respective trait
  d.target.patch <- d.cwm.patch.nulldiff %>% 
    filter(!is.na(!! sym(trait_i)))
  
  formula_patch <- 
    as.formula(paste(trait_i, " ~ DWV + Cover_above2m + ssc_mean"))
  
  # linear model
  mod.patch <- lm(formula_patch, data = d.target.patch, na.action = na.fail)
  
  p_patch <- summary(mod.patch)$coefficients["(Intercept)",c("Estimate", "Pr(>|t|)")]
  
  
  # Model object-level ---------------------------------------------------------.
  
  # subset for the respective trait
  d.target.object <- d.cwm.object.nulldiff %>% 
    filter(!is.na(!! sym(trait_i)))
  
  formula_object <- 
    as.formula(paste(trait_i, " ~ Period * Tree_sp  + (1 | PlotID/TreeID)"))
  
  # change contrasts to sum (necessary for later plotting)
  options(contrasts = c("contr.sum", "contr.poly"))
  
  # linear mixed effects model
  mod.object <- glmmTMB(formula_object,
                        na.action = na.fail, 
                        data = d.target.object)
  
  p_object <- summary(mod.object)$coefficients$cond["(Intercept)", c("Estimate", "Pr(>|z|)")]
  
  c(p_site = ifelse(p_site[2]<=0.05, sign(p_site[1]), 0), 
    p_patch = ifelse(p_patch[2]<=0.05, sign(p_patch[1]), 0), 
    p_object = ifelse(p_object[2]<=0.05, sign(p_object[1]), 0))
  
}


l.intercept <- lapply(v.sel.traits.ext,
                      f_intercept) 
l.intercept <- do.call(rbind, l.intercept)
rownames(l.intercept) <- v.sel.traits.ext
l.intercept


# detailed plots (supplement) ##################################################
################################################################################.

# plotting order of filtering variables
var.order <- c("prop_forest2000", "DWV", "NMDS1", "NMDS2", "ssc_mean",
               "Canopy_closure", "Cover_above2m", "Period1")

# filtering variables at site level
vars.site <- c("prop_forest2000", "DWV", "NMDS1", "NMDS2", "ssc_mean",
               "Canopy_closure")

# filtering variables at patch level
vars.patch <- c("DWV", "ssc_mean", "Cover_above2m")

vars.object <- c("Period1", 
                 paste0(c("", "Period1:"), rep(tree.order[-length(tree.order)], each = 2)))


var.labels <- vars.object
names(var.labels) <- vars.object
var.labels <- gsub("Period1", "1st period", var.labels)
var.labels <- gsub(":", " × ", var.labels)

var.labels <- c(prop_forest2000 = "Forest surroundings",
                DWV = "Deadwood volume",
                NMDS1 = "Broadleaf share",
                NMDS2 = "Conifer type",
                ssc_mean = "Structural complexity",
                Canopy_closure = "Canopy closure site",
                Cover_above2m = "Canopy closure patch",
                var.labels,
                Pseudotsuga = "Pseudotsuga",
                "Period1:Pseudotsuga" = "1st period × Pseudotsuga")

trait.labels <- c(wing.length.res= "Wing length",
                  wing.loading.rev= "Wing loading",
                  wing.aspect = "Wing aspect",
                  body.length = "Body length",
                  leg.length.res = "Leg length",
                  lightness = "Lightness",
                  antenna.length.res = "Antenna length",
                  eye.length.res = "Eye length",
                  hairiness = "Hairiness",
                  body.width.res = "Body width",
                  body.roundness = "Body roundness",
                  head.length.res= "Head length",
                  mandibular.aspect = "Mandibular aspect")

trait.labels.ext <- c(paste(trait.labels, "(CWM)"), paste(trait.labels, "(FD)"))
names(trait.labels.ext) <- c(paste0(names(trait.labels), ".CWM"), paste0(names(trait.labels), ".FD"))

trait.labels.ext <- c(trait.labels.ext, FDRao = "Functional diversity")

f_levelplot_detail <- function(trait_i){
  
  # Model site-level -----------------------------------------------------------.
  
  # subset for the respective trait
  d.target.site <- d.cwm.site.nulldiff %>% 
    filter(!is.na(!! sym(trait_i)))
  
  formula_site <- 
    as.formula(paste(trait_i, " ~ NMDS1 + NMDS2 + prop_forest2000 + DWV +
          Canopy_closure + ssc_mean"))
  
  # linear model
  mod.site <- lm(formula_site, data = d.target.site, na.action = na.fail)
  
  # Model patch-level ----------------------------------------------------------.
  
  # subset for the respective trait
  d.target.patch <- d.cwm.patch.nulldiff %>% 
    filter(!is.na(!! sym(trait_i)))
  
  formula_patch <- 
    as.formula(paste(trait_i, " ~ DWV + Cover_above2m + ssc_mean"))
  
  # linear model
  mod.patch <- lm(formula_patch, data = d.target.patch, na.action = na.fail)
  
  # Model object-level ---------------------------------------------------------.
  
  # subset for the respective trait
  d.target.object <- d.cwm.object.nulldiff %>% 
    filter(!is.na(!! sym(trait_i)))
  
  formula_object <- 
    as.formula(paste(trait_i, " ~ Period * Tree_sp  + (1 | PlotID/TreeID)"))
  
  # change contrasts to sum (necessary for later plotting)
  options(contrasts = c("contr.sum", "contr.poly"))
  
  # linear mixed effects model
  mod.object <- glmmTMB(formula_object,
                        na.action = na.fail, 
                        data = d.target.object)
  
  # extract coefficients
  d.plotting <- summary(mod.site)$coefficients %>% 
    as.data.frame() %>% 
    rownames_to_column("var") %>% 
    mutate(level = "Site") %>% 
    bind_rows(summary(mod.patch)$coefficients %>% 
                as.data.frame() %>% 
                rownames_to_column("var") %>% 
                mutate(level = "Patch")) %>% 
    rename(p_val = "Pr(>|t|)",
           se = "Std. Error") %>% 
    bind_rows(summary(mod.object)$coefficients$cond %>% 
                as.data.frame() %>% 
                rownames_to_column("var") %>% 
                mutate(level = "Object") %>% 
                rename(p_val = "Pr(>|z|)",
                       se = "Std. Error"))
  
  
  # change tree numbers (from sum contrasts) to names
  for (i in 1:nlevels(d.cwm.object.nulldiff$Tree_sp)){
    d.plotting <- d.plotting %>% 
      mutate(var = ifelse(var %in% c(paste0("Tree_sp", i), paste0("Period1:Tree_sp", i)),
                          gsub(paste0("Tree_sp", i), 
                               tree.order[levels(d.cwm.object.nulldiff$Tree_sp)[i]], 
                               var),
                          var))
  }
  
  # extract intercepts
  intercept_i <- d.plotting %>% 
    filter(var == "(Intercept)") %>% 
    select(Estimate, level) %>% 
    rename(intercept = Estimate)
  
  # add values for missing species
  d.plotting$lastlevel <- NA    
  
  # main effect
  d.plotting <- d.plotting %>% 
    filter(var %in% tree.order) %>% 
    summarise(Estimate = -sum(Estimate),
              se = mean(se), # approximation
              var = tree.order[length(tree.order)],
              level = "Object",
              lastlevel = T) %>% 
    mutate(p_val = ifelse(Estimate - 1.96 * se > 0 |
                            Estimate + 1.96 * se < 0,
                          0.05, 1)) %>% # p value approximation (confidence intervals crossing zero?)
    bind_rows(d.plotting, .) 
  
  # interaction effect
  d.plotting <- d.plotting %>% 
    filter(var %in% paste0("Period1:",tree.order)) %>% 
    summarise(Estimate = -sum(Estimate),
              se = mean(se), # approximation
              var = paste0("Period1:", tree.order[length(tree.order)]),
              level = "Object",
              lastlevel = T) %>% 
    mutate(p_val = ifelse(Estimate - 1.96 * se > 0 |
                            Estimate + 1.96 * se < 0,
                          0.05, 1)) %>% # p value approximation (confidence intervals crossing zero?)
    bind_rows(d.plotting, .)
  
  
  
  
  mns_i <- data.frame(mean = c(d.nullmeans.site[, paste0(trait_i, "_mn")],
                               d.nullmeans.patch[, paste0(trait_i, "_mn")],
                               d.nullmeans.object[, paste0(trait_i, "_mn")]),
                      level = c("Site", "Patch", "Object"))
  
  
  density.null.site <- data.frame(density(d.cwm.site.null[, trait_i])[c("x", "y")])
  
  range_hdi_site <- bayestestR::hdi(d.cwm.site.null[, trait_i], ci = .95)
  
  density.null.site <- density.null.site %>% 
    filter(x >= range_hdi_site$CI_low,
           x <= range_hdi_site$CI_high) %>% 
    mutate(y = y / max(y), # scale to max 1
           x_diff = median(diff(x))) # add differences for later plotting (width of rects)
  
  density.null.patch <- data.frame(density(d.cwm.patch.null[, trait_i])[c("x", "y")])
  
  range_hdi_patch <- bayestestR::hdi(d.cwm.patch.null[, trait_i], ci = .95)
  
  density.null.patch <- density.null.patch %>% 
    filter(x >= range_hdi_patch$CI_low,
           x <= range_hdi_patch$CI_high) %>% 
    mutate(y = y / max(y), # scale to max 1
           x_diff = median(diff(x))) # add differences for later plotting (width of rects)
  
  density.null.object <- data.frame(density(d.cwm.object.null[, trait_i])[c("x", "y")])
  
  range_hdi_object <- bayestestR::hdi(d.cwm.object.null[, trait_i], ci = .95)
  
  density.null.object <- density.null.object %>% 
    filter(x >= range_hdi_object$CI_low,
           x <= range_hdi_object$CI_high) %>% 
    mutate(y = y / max(y), # scale to max 1
           x_diff = median(diff(x))) # add differences for later plotting (width of rects)
  
  
  
  
  
  d.plotting <- d.plotting %>% 
    filter(var != "(Intercept)") %>% 
    left_join(mns_i, by = "level") %>% 
    left_join(intercept_i, by = "level") %>% 
    mutate(Estimate = mean + intercept + Estimate,
           level = ifelse(lastlevel & !is.na(lastlevel), "Lastlevel", level),
           level = factor(level, levels = c("Site", "Patch", "Object", "Lastlevel")),
           var = factor(var, levels = c(var.order, paste0(c("", "Period1:"), rep(tree.order, each = 2)))),
           sign = ifelse(p_val <= .05, "sign", "ns"),
           sign = ifelse(is.na(sign), "ns", sign),
           varlevel = interaction(var, level),
           varlevel = factor(varlevel, levels = c(paste0(vars.site, ".Site"),
                                                  paste0(vars.patch, ".Patch"),
                                                  paste0(vars.object, ".Object"),
                                                  paste0(tree.order[length(tree.order)], ".Lastlevel"),
                                                  paste0("Period1:", tree.order[length(tree.order)], ".Lastlevel"))))
  
  for (i in tree.order){
    d.plotting <- d.plotting %>% 
      mutate(Estimate = ifelse(grepl(paste0("Period1:", i), var),
                               Estimate + d.plotting$Estimate[d.plotting$var == i] - intercept - mean,
                               Estimate))
  }
  
  
  xlabels <- var.labels[gsub(".Site|.Patch|.Object|.Lastlevel", "", levels(d.plotting$varlevel))]
  names(xlabels) <- levels(d.plotting$varlevel)
  
  d.plotting %>% 
    ggplot() +
    geom_point(aes(x = varlevel, y = Estimate, col = level))  +
    geom_rect(data = density.null.site,
              aes(ymin = x+x_diff/2, ymax = x-x_diff/2, 
                  xmin = 0, xmax = length(unique(vars.site)) + .5, alpha = y), 
              fill = "grey60", col = NA) +
    geom_segment(data = NULL, aes(x = 0, xend = length(vars.site) + .5, 
                                  y = mns_i$mean[1] + intercept_i$intercept[1], 
                                  yend = mns_i$mean[1] + intercept_i$intercept[1],
                                  col = "Site"), lty = 1) +
    geom_rect(data = density.null.patch,
              aes(ymin = x+x_diff/2, ymax = x-x_diff/2, 
                  xmin = length(vars.site) + .5, 
                  xmax = length(vars.site) + length(vars.patch) + .5, 
                  alpha = y), 
              fill = "grey60", col = NA) +
    geom_segment(data = NULL, aes(x = length(vars.site) + .5, 
                                  xend = length(vars.site) + length(vars.patch) + .5, 
                                  y = mns_i$mean[2] + intercept_i$intercept[2], 
                                  yend = mns_i$mean[2] + intercept_i$intercept[2],
                                  col = "Patch"), lty = 1) +
    geom_rect(data = density.null.object,
              aes(ymin = x+x_diff/2, ymax = x-x_diff/2, 
                  xmin = length(vars.site) + length(vars.patch) + .5, 
                  xmax = Inf, 
                  alpha = y), 
              fill = "grey60", col = NA) +
    geom_segment(data = NULL, aes(x = length(vars.site) + length(vars.patch) + .5, 
                                  xend = Inf, 
                                  y = mns_i$mean[3] + intercept_i$intercept[3], 
                                  yend = mns_i$mean[3] + intercept_i$intercept[3],
                                  col = "Object"), lty = 1) +
    geom_segment(aes(x = varlevel, y = Estimate, col = level,
                     yend = Estimate - 1.96 * se, xend = varlevel,
                     size = sign)) +
    geom_segment(aes(x = varlevel, y = Estimate, col = level,
                     yend = Estimate + 1.96 * se, xend = varlevel,
                     size = sign)) +
    geom_point(aes(x = varlevel, y = Estimate, shape = level),
               col = "white", size = 2.75) +
    geom_point(aes(x = varlevel, y = Estimate, col = level, shape = level),
               size = 2) +
    geom_vline(xintercept = seq(.5, 8.5, 1), lty = 3, size = .25) +
    geom_vline(xintercept = seq(12.5, 50.5, 2), lty = 3, size = .25) +
    geom_vline(xintercept = 10.5, lty = 2, size = .25) +
    geom_vline(xintercept = c(6.5, 9.5), lty = 1, size = .5) +
    scale_fill_manual(values = c(Site = "#1b9e77", Patch = "#d95f02", Object = "#7570b3", Lastlevel = "grey40"),
                      guide = F) +
    scale_colour_manual(values = c(Site = "#1b9e77", Patch = "#d95f02", Object = "#7570b3", Lastlevel = "grey40"),
                        breaks = c("Site", "Patch", "Object"), 
                        labels = c("Region-to-site", "Site-to-patch", "Patch-to-object"),
                        name = "Filtering step") +
    scale_shape_manual(values = c(Site = 19, Patch = 15, Object = 17), 
                       labels = c("Region-to-site", "Site-to-patch", "Patch-to-object"),
                       name = "Filtering step") +
    scale_size_manual(values = c(sign = 1.5, ns = .5), guide = F) +
    scale_x_discrete(labels = xlabels, drop = F) +
    scale_alpha_continuous(range = c(.1, .8), guide = F) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
          axis.title.x = element_blank()) +
    ylab(trait.labels.ext[trait_i]) +
    # add interaction lines
    geom_segment(data = . %>% filter(var %in% tree.order),
                 aes(x = as.numeric(var) + 1.5, y = Estimate, xend = as.numeric(var) + 3.5, yend = Estimate, col = level),
                 lty = 1)   
  
  
  
}

for (trait_i in c("FDRao", v.sel.traits.ext)){
  p <- f_levelplot_detail(trait_i)
}
