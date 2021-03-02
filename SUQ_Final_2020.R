
##########################
##########################
#### BBP SUQ ANALYSES ####
##########################
##########################


###########################
#### REQUIRED PACKAGES ####
###########################


library(stringr)
library(dplyr)
library(ggplot2)
library(magrittr)

###################
#### KEY FILES ####
###################

load("/Users/whitmanet/Downloads/BBP_SUQ.analyses/BBP_SUQ.RData")

SUQ.82 # 82 people under aged 12 with SUQ data in BBP
scv_longit_acghst # MAGeT volumes for all NVs and SCAs that passed QC
SUQ.82_maget  # 82 people under aged 12 with SUQ data in BBP appended with MAGeT volumes (primary data set)

# save.image('~/Desktop/BBP_SUQ.RData')

####################################
#### LOCATIONS OF CRITICAL DATA ####
####################################

# all MAGeT brain primary output files
/Volumes/bv1/scinet/all_maget_output # segmentations
/Volumes/bv3/AK/nih_subcort_sa/  # surface vertex files with seperate folder per structure named: amygdala_left | gp_left |  new_lefthc |striatum_right | amygdala_right |  gp_right | new_righthc |	thalamus_left	| striatum_left	| thalamus_right

#######################
#### KEY VARIABLES ####
#######################

# Alcohol
  
SUQ__2c # drinking days per month at peak use
SUQ__2d # drinks per drinking day at peak use
DrinksMonth # product of 2c and 2d, note outlier with value >300
DrinksMonth.nona # product of 2c and 2d, plugging in 0 values for individuals who have never consumed alcohol

##########################################################################################
#### DESCRIPTIVE CHARACTERISTICS AND DIFFERENCES BETWEEN MALES AND FEMALES (TALBLE 1) ####
##########################################################################################

# Demographic Characteristics 
table(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,]$SEX)

with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(AGESCAN, SEX, mean))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(AGESCAN, SEX, sd))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], t.test(AGESCAN~SEX))

with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(IQ, SEX, mean))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(IQ, SEX, sd))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], t.test(IQ~SEX))

with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(SES, SEX, mean))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(SES, SEX, sd))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], t.test(SES~SEX))

# Alcohol Use Characteristics


## Summary of total elapsed time
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], summary(AGE_AT_BBP - AGESCAN))

with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], summary(AGE_AT_BBP))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(AGE_AT_BBP, SEX, mean))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(AGE_AT_BBP, SEX, sd))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], t.test(AGE_AT_BBP~SEX))

with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], summary(AGE_AT_PEAK))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(AGE_AT_PEAK, SEX, mean, na.rm = TRUE))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], tapply(AGE_AT_PEAK, SEX, sd, na.rm = TRUE))
with(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], t.test(AGE_AT_PEAK~SEX, na.rm = TRUE))


################################################
#### SETTING UP MAIN GLIM FILE FOR ANALYSIS ####
################################################

# confirm all SUQ.82 scans passed MAGeT QC

sum(SUQ.82$MASKID %in% scv_longit_acghst$MASKID)

# combine SUQ.82 with MAGeT values

SUQ.82_maget <- merge(SUQ.82, scv_longit_acghst, by="MASKID")



##########################################################################################
#### ADDED NOV 2020 - REACALCULATING ARE AND NUCLEAR AREA VARIABLES FROM VERTEX FILES ####
##########################################################################################

SUQ.82_maget$DrinksMonth.nona <- SUQ.82_maget$DrinksMonth
SUQ.82_maget$DrinksMonth.nona[is.na(SUQ.82_maget$DrinksMonth.nona)] <- 0

SUQ.82_maget$age.cent <- scale(SUQ.82_maget$AGESCAN, center=T, scale=F) %>% as.vector

SUQ.82_maget$bil.amygdala.a <- colSums(rbind(amygdala_vertex.sa_l, amygdala_vertex.sa_r))
SUQ.82_maget$bil.hippocampus.a <- colSums(rbind(hippocampus_vertex.sa_l, hippocampus_vertex.sa_r))

SUQ.82_maget$bil.amygdala.a_resid <- with(SUQ.82_maget, lm(bil.amygdala.a~ AGESCAN + SEX)) %>% use_series(residuals)
SUQ.82_maget$bil.hippocampus.a_resid <- with(SUQ.82_maget, lm(bil.hippocampus.a~ AGESCAN + SEX)) %>% use_series(residuals)

SUQ.82_maget$bil.amygdala.v_resid <- with(SUQ.82_maget, lm(bil.amygdala.v~ AGESCAN + SEX)) %>% use_series(residuals)
SUQ.82_maget$bil.hipp.v_resid <- with(SUQ.82_maget, lm(bil.hipp.v~ AGESCAN + SEX)) %>% use_series(residuals)

#############################################################
#### TOTAL SUBCORTICAL VOLUME AND SURFACE AREA (Table 2) ####
#############################################################

#### Amygdala ####

#volume
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.v~age.cent + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.v~age.cent + SEX*alcohol.cent))


#surface area
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.a~age.cent + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.a~age.cent + SEX*alcohol.cent))


#### Hippocampus ####

#volume
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hipp.v~age.cent + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hipp.v~age.cent + SEX*alcohol.cent))

#surface area
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hippocampus.a~age.cent + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hippocampus.a~age.cent + SEX*alcohol.cent))

#### Striatum ####

#volume
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.v~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.v~age.cent + SEX+DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.v~age.cent + SEX*alcohol.cent))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.v~age.cent + SEX+alcohol.cent))

#surface area 
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.a~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.a~age.cent + SEX+DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.a~age.cent + SEX*alcohol.cent))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.a~age.cent + SEX+alcohol.cent))

###################################################################
#### CONTROLLING FOR TOTAL BRAIN VOLUME - ADDED FEB 2021 ##########
###################################################################

liv_vols <- read.csv("/Users/whitmanet/Documents/long_alcohol/CatherineM_etoh_n82_civet10_csf_gmv_wmv.csv")
liv_vols[,5] <- liv_vols$GMV + liv_vols$WMV
liv_vols[,6] <- liv_vols$GMV + liv_vols$WMV + liv_vols$CSF

colnames(liv_vols) <- c("MASKID", "CSF", "GMV", "WMV", "GMV_WMV", "GMV_WMV_CSF")

SUQ.82_maget <- merge(SUQ.82_maget, liv_vols, by="MASKID")

# sanity check
cor.test(SUQ.82_maget$GMV, SUQ.82_maget$WMV)

cor.test(SUQ.82_maget$bil.amygdala.v, SUQ.82_maget$GMV_WMV)
cor.test(SUQ.82_maget$bil.hipp.v, SUQ.82_maget$GMV_WMV)
cor.test(SUQ.82_maget$bil.striatum.v, SUQ.82_maget$GMV_WMV)

cor.test(SUQ.82_maget$bil.amygdala.a, SUQ.82_maget$GMV_WMV)
cor.test(SUQ.82_maget$bil.hippocampus.a, SUQ.82_maget$GMV_WMV)
cor.test(SUQ.82_maget$bil.striatum.a, SUQ.82_maget$GMV_WMV)

cor.test(SUQ.82_maget$bil.amygdala.v_resid, SUQ.82_maget$GMV_WMV)
cor.test(SUQ.82_maget$bil.hipp.v_resid, SUQ.82_maget$GMV_WMV)
# cor.test(SUQ.82_maget$bil.striatum.v_resid, SUQ.82_maget$GMV_WMV)

cor.test(SUQ.82_maget$bil.amygdala.a_resid, SUQ.82_maget$GMV_WMV)
cor.test(SUQ.82_maget$bil.hippocampus.a_resid, SUQ.82_maget$GMV_WMV)
cor.test(SUQ.82_maget$bil.striatum.a_resid, SUQ.82_maget$GMV_WMV)


#### Amygdala ####

#volume
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.v~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.v~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.v~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.v~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.v~age.cent + GMV_WMV + SEX*alcohol.cent))


#surface area
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.a~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.a~age.cent + GMV +SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.a~age.cent + GMV_WMV +SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.a~age.cent + GMV_WMV_CSF +SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.amygdala.a~age.cent + GMV_WMV +  SEX*alcohol.cent))


#### Hippocampus ####

#volume
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hipp.v~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hipp.v~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hipp.v~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hipp.v~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hipp.v~age.cent + GMV_WMV + SEX*alcohol.cent))

#surface area
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hippocampus.a~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hippocampus.a~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hippocampus.a~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hippocampus.a~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.hippocampus.a~age.cent + GMV_WMV + SEX*alcohol.cent))

#### Striatum ####

#volume
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.v~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.v~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.v~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.v~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))

#surface area 
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.a~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.a~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.a~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], bil.striatum.a~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))



#### while excluding outliers



#### Amygdala ####

#volume
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.v~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.v~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.v~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.v~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.v~age.cent + GMV_WMV + SEX*alcohol.cent))


#surface area
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.a~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.a~age.cent + GMV +SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.a~age.cent + GMV_WMV +SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.a~age.cent + GMV_WMV_CSF +SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.a~age.cent + GMV_WMV +  SEX*alcohol.cent))


#### Hippocampus ####

#volume
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hipp.v~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hipp.v~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hipp.v~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hipp.v~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hipp.v~age.cent + GMV_WMV + SEX*alcohol.cent))

#surface area
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hippocampus.a~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hippocampus.a~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hippocampus.a~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hippocampus.a~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))

summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hippocampus.a~age.cent + GMV_WMV + SEX*alcohol.cent))

#### Striatum ####

#volume
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.v~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.v~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.v~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.v~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))

#surface area 
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.a~age.cent + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.a~age.cent + GMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.a~age.cent + GMV_WMV + SEX*DrinksMonth.nona))
summary(lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.a~age.cent + GMV_WMV_CSF + SEX*DrinksMonth.nona))


###########################################
#### VERTEX-LEVEL ANALYSES OF PEAK USE ####
###########################################

#### Amygdala ####

# making left and right vertex tables 

path_to_dir <- "/Volumes/bv3/AK/nih_subcort_sa/amygdala_left"
p <- data.frame(path=armin.path.to.subcort.text.files(SUQ.82_maget$MASKID, path_to_dir))
amygdala_vertex.sa_l <- mni.build.data.table(p)

p$path <- gsub("amygdala_left",  "amygdala_right", p$path)
amygdala_vertex.sa_r <- mni.build.data.table(p)

# calculating total amygdala surface area and plotting against total amygdala volume

SUQ.82_maget$bil.amygdala.a <- colSums(rbind(amygdala_vertex.sa_l, amygdala_vertex.sa_r))
with(SUQ.82_maget, plot(bil.amygdala.v, bil.amygdala.a))


# modeling vertex values by substance use metrics


### with sex interaction

SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX*DrinksMonth.nona', amygdala_vertex.sa_l[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l.vertstats')

SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX*DrinksMonth.nona', amygdala_vertex.sa_r[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r.vertstats')

### without sex interaction

SUQ.82_age.cent_sex_DrinksMonth.nona_amygdala.sa_l <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX + DrinksMonth.nona', amygdala_vertex.sa_l[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex_DrinksMonth.nona_amygdala.sa_l, '~/Desktop/SUQ.82_age.cent_sex_DrinksMonth.nona_amygdala.sa_l.vertstats')

SUQ.82_age.cent_sex_DrinksMonth.nona_amygdala.sa_r <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX + DrinksMonth.nona', amygdala_vertex.sa_r[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex_DrinksMonth.nona_amygdala.sa_r, '~/Desktop/SUQ.82_age.cent_sex_DrinksMonth.nona_amygdala.sa_r.vertstats')

## FDR for the amygdala

# seeing the t statistic column positions so you can select which ones you want to fdr over
attributes(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic)

# are any fdr corrected p values for the term of interest (here column #5) below 0.05? yes
c(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic[,5], SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic[,5]) %>% abs %>% pt(., 76, lower.tail=F) %>% p.adjust(., method="fdr") %>% summary() ## 76 = degrees of freedom = 81 - 5

## find vertices with significant sex interaction with FDR correction

SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic <- cbind(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic, c(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic[,5], SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic[,5]) %>% abs %>% pt(., 76, lower.tail=F) %>% p.adjust(., method="fdr") %>% extract(1:1473))
colnames(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic)[6] <- "fdr.corr_p"

SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic <- cbind(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic, c(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic[,5], SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic[,5]) %>% abs %>% pt(., 76, lower.tail=F) %>% p.adjust(., method="fdr") %>% extract(1474:2878))
colnames(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic)[6] <- "fdr.corr_p"

SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic <- cbind(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic, -log10(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic[,6]))
colnames(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l$tstatistic)[7] <- "neg.log10.p"

SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic <- cbind(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic, -log10(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic[,6]))
colnames(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r$tstatistic)[7] <- "neg.log10.p"

mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_l.vertstats')
mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_amygdala.sa_r.vertstats')

# dropping down to simpler model 

# seeing the t statistic column positions so you can select which ones you want to fdr over
attributes(SUQ.82_age.cent_sex_DrinksMonth.nona_amygdala.sa_r$tstatistic)

# are any fdr corrected p values for the term of interest (here column #5) below 0.05? no

c(SUQ.82_age.cent_sex_DrinksMonth.nona_amygdala.sa_l$tstatistic[,4], SUQ.82_age.cent_sex_DrinksMonth.nona_amygdala.sa_r$tstatistic[,4]) %>% abs %>% pt(., 77, lower.tail=F) %>% p.adjust(., method="fdr") %>% summary() ## 77 = degrees of freedom = 81 - 4

#### HIPPOCAMPUS ####

# making left and right vertex tables 

path_to_dir <- "/Volumes/bv3/AK/nih_subcort_sa/new_lefthc"
p <- data.frame(path=armin.path.to.subcort.text.files(SUQ.82_maget$MASKID, path_to_dir))
hippocampus_vertex.sa_l <- mni.build.data.table(p)

p$path <- gsub("new_lefthc",  "new_righthc", p$path)
hippocampus_vertex.sa_r <- mni.build.data.table(p)

# calculating total hippocampal surface area and plotting against total hippocampal volume

SUQ.82_maget$bil.hippocampus.a <- colSums(rbind(hippocampus_vertex.sa_l, hippocampus_vertex.sa_r))
with(SUQ.82_maget, plot(bil.hipp.v, bil.hippocampus.a))

### with sex interaction

SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX*DrinksMonth.nona', hippocampus_vertex.sa_l[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l.vertstats')

SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX*DrinksMonth.nona', hippocampus_vertex.sa_r[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r.vertstats')

### without sex interaction

SUQ.82_age.cent_sex_DrinksMonth.nona_hippocampus.sa_l <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX + DrinksMonth.nona', hippocampus_vertex.sa_l[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex_DrinksMonth.nona_hippocampus.sa_l, '~/Desktop/SUQ.82_age.cent_sex_DrinksMonth.nona_hippocampus.sa_l.vertstats')

SUQ.82_age.cent_sex_DrinksMonth.nona_hippocampus.sa_r <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX + DrinksMonth.nona', hippocampus_vertex.sa_r[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex_DrinksMonth.nona_hippocampus.sa_r, '~/Desktop/SUQ.82_age.cent_sex_DrinksMonth.nona_hippocampus.sa_r.vertstats')

## FDR for the hippocampus

# seeing the t statistic column positions so you can select which ones you want to fdr over
attributes(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic)

# are any fdr corrected p values for the term of interest (here column #5) below 0.05? yes
c(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic[,5], SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic[,5]) %>% abs %>% pt(., 76, lower.tail=F) %>% p.adjust(., method="fdr") %>% summary() ## 76 = degrees of freedom = 81 - 5

## find vertices with significant sex interaction with FDR correction

SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic <- cbind(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic, c(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic[,5], SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic[,5]) %>% abs %>% pt(., 76, lower.tail=F) %>% p.adjust(., method="fdr") %>% extract(1:1152))
colnames(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic)[6] <- "fdr.corr_p"

SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic <- cbind(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic, c(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic[,5], SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic[,5]) %>% abs %>% pt(., 76, lower.tail=F) %>% p.adjust(., method="fdr") %>% extract(1153:2367))
colnames(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic)[6] <- "fdr.corr_p"

SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic <- cbind(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic, -log10(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic[,6]))
colnames(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l$tstatistic)[7] <- "neg.log10.p"

SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic <- cbind(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic, -log10(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic[,6]))
colnames(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r$tstatistic)[7] <- "neg.log10.p"

mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_l.vertstats')
mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_hippocampus.sa_r.vertstats')


# dropping down to simpler model 

# seeing the t statistic column positions so you can select which ones you want to fdr over
attributes(SUQ.82_age.cent_sex_DrinksMonth.nona_hippocampus.sa_r$tstatistic)

# are any fdr corrected p values for the term of interest (here column #5) below 0.05? no

c(SUQ.82_age.cent_sex_DrinksMonth.nona_hippocampus.sa_l$tstatistic[,4], SUQ.82_age.cent_sex_DrinksMonth.nona_hippocampus.sa_r$tstatistic[,4]) %>% abs %>% pt(., 77, lower.tail=F) %>% p.adjust(., method="fdr") %>% summary() ## 77 = degrees of freedom = 81 - 4


#### STRIATUM ####

# making left and right vertex tables 

path_to_dir <- "/Volumes/bv3/AK/nih_subcort_sa/striatum_left"
p <- data.frame(path=armin.path.to.subcort.text.files(SUQ.82_maget$MASKID, path_to_dir))
striatum_vertex.sa_l <- mni.build.data.table(p)

p$path <- gsub("striatum_left",  "striatum_right", p$path)
striatum_vertex.sa_r <- mni.build.data.table(p)

# calculating total hippocampal surface area and plotting against total hippocampal volume

SUQ.82_maget$bil.striatum.a <- colSums(rbind(striatum_vertex.sa_l, striatum_vertex.sa_r))
with(SUQ.82_maget, plot(bil.hipp.v, bil.striatum.a))


# modeling vertex values by substance use metrics


### with sex interaction

SUQ.82_age.cent_sex.DrinksMonth.nona_striatum.sa_l <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX*DrinksMonth.nona', striatum_vertex.sa_l[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_striatum.sa_l, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_striatum.sa_l.vertstats')

SUQ.82_age.cent_sex.DrinksMonth.nona_striatum.sa_r <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX*DrinksMonth.nona', striatum_vertex.sa_r[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex.DrinksMonth.nona_striatum.sa_r, '~/Desktop/SUQ.82_age.cent_sex.DrinksMonth.nona_striatum.sa_r.vertstats')

### without sex interaction

SUQ.82_age.cent_sex_DrinksMonth.nona_striatum.sa_l <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX + DrinksMonth.nona', striatum_vertex.sa_l[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex_DrinksMonth.nona_striatum.sa_l, '~/Desktop/SUQ.82_age.cent_sex_DrinksMonth.nona_striatum.sa_l.vertstats')

SUQ.82_age.cent_sex_DrinksMonth.nona_striatum.sa_r <- mni.vertex.statistics(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], 'y~ age.cent + SEX + DrinksMonth.nona', striatum_vertex.sa_r[,SUQ.82_maget$DrinksMonth.nona < 300])
mni.write.vertex.stats(SUQ.82_age.cent_sex_DrinksMonth.nona_striatum.sa_r, '~/Desktop/SUQ.82_age.cent_sex_DrinksMonth.nona_striatum.sa_r.vertstats')

### FDR for striatum (NOTE: NO MODEL SURVIVES FDR CORRECTION)

# seeing the t statistic column positions so you can select which ones you want to fdr over
attributes(SUQ.82_age.cent_sex.DrinksMonth.nona_striatum.sa_r$tstatistic)

# are any fdr corrected p values for the term of interest (here column #5) below 0.05? no
c(SUQ.82_age.cent_sex.DrinksMonth.nona_striatum.sa_l$tstatistic[,5], SUQ.82_age.cent_sex.DrinksMonth.nona_striatum.sa_r$tstatistic[,5]) %>% abs %>% pt(., 76, lower.tail=F) %>% p.adjust(., method="fdr") %>% summary() ## 76 = degrees of freedom = 81 - 5

# dropping down to simpler model 

# seeing the t statistic column positions so you can select which ones you want to fdr over
attributes(SUQ.82_age.cent_sex_DrinksMonth.nona_striatum.sa_r$tstatistic)

# are any fdr corrected p values for the term of interest (here column #5) below 0.05? no

c(SUQ.82_age.cent_sex_DrinksMonth.nona_striatum.sa_l$tstatistic[,4], SUQ.82_age.cent_sex_DrinksMonth.nona_striatum.sa_r$tstatistic[,4]) %>% abs %>% pt(., 77, lower.tail=F) %>% p.adjust(., method="fdr") %>% summary() ## 77 = degrees of freedom = 81 - 4

############################################################
#### VISUALIZATION OF MOST IMPACTED VERTICES (Figure 3) ####
############################################################

### For Hippocampus
SUQ.82_maget$tmp <- hippocampus_vertex.sa_l[1071+1,]
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x=DrinksMonth.nona, y=tmp, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "LEFT HIPPOCAMPAL VOLUME AT VERTEX 1071") + theme(text = element_text(size=15))

SUQ.82_maget$tmp <- hippocampus_vertex.sa_r[1159+1,]
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x=DrinksMonth.nona, y=tmp, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "RIGHT HIPPOCAMPAL VOLUME AT VERTEX 1159") + theme(text = element_text(size=15))

### For Amygdala
SUQ.82_maget$tmp <- amygdala_vertex.sa_l[59+1,]
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x=DrinksMonth.nona, y=tmp, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "LEFT AMYGDALA VOLUME AT VERTEX 59") + theme(text = element_text(size=15))

SUQ.82_maget$tmp <- amygdala_vertex.sa_r[472+1,]
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x=DrinksMonth.nona, y=tmp, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "RIGHT AMYGDALA VOLUME AT VERTEX 472") + theme(text = element_text(size=15))

### Sensitivity Analysis: limitting drinks to < 100

### For Hippocampus
SUQ.82_maget$tmp <- hippocampus_vertex.sa_l[1071+1,]
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 75,], aes(x=DrinksMonth.nona, y=tmp, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "LEFT HIPPOCAMPAL VOLUME AT VERTEX 1071") + theme(text = element_text(size=15))

SUQ.82_maget$tmp <- hippocampus_vertex.sa_r[1159+1,]
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 75,], aes(x=DrinksMonth.nona, y=tmp, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "RIGHT HIPPOCAMPAL VOLUME AT VERTEX 1159") + theme(text = element_text(size=15))

### For Amygdala
SUQ.82_maget$tmp <- amygdala_vertex.sa_l[59+1,]
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 75,], aes(x=DrinksMonth.nona, y=tmp, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "LEFT AMYGDALA VOLUME AT VERTEX 59") + theme(text = element_text(size=15))

SUQ.82_maget$tmp <- amygdala_vertex.sa_r[472+1,]
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 75,], aes(x=DrinksMonth.nona, y=tmp, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "RIGHT AMYGDALA VOLUME AT VERTEX 472") + theme(text = element_text(size=15))


#### MAKING AND PLOTTING RESIDUALIZED VALUES ####

## e.g. for amygdala

# make variable
SUQ.82_maget$bil.amygdala.a_resid <- with(SUQ.82_maget, lm(bil.amygdala.a~ AGESCAN + SEX)) %>% use_series(residuals)

# plot
filter(SUQ.82_maget, DrinksMonth.nona<300) %>% ggplot(., aes(x=DrinksMonth.nona, y=bil.amygdala.a_resid, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") + facet_grid(~SEX)


#### AGE GRAPH #####

ggplot() +
  geom_segment(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x = AGESCAN, xend = AGE_AT_BBP, y = reorder(count, AGESCAN), yend = reorder(count, AGESCAN), color = SEX, size = .1)) +
  geom_point(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x = AGE_AT_PEAK, y = reorder(count, AGESCAN), color = SEX, size = DrinksMonth.nona)) +
  labs(x="AGE", y="INDIVIDUAL", color = "SEX", size = "DRINKS PER MONTH") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position ="right", axis.text.x=element_text(size=30), axis.title.x=element_text(size=30, vjust=-0.3), axis.title.y=element_text(size=30, vjust=0.3)) +
  theme(text = element_text(size=30))

ggplot() +
  geom_segment(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x = AGESCAN, xend = AGE_AT_BBP, y = reorder(count, AGESCAN), yend = reorder(count, AGESCAN), color = SEX)) +
  geom_point(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x = AGE_AT_PEAK, y = reorder(count, AGESCAN), color = SEX, size = DrinksMonth.nona)) +
  labs(x="AGE", y="INDIVIDUAL", color = "SEX", size = "DRINKS PER MONTH") + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position ="right", axis.text.x=element_text(size=18), axis.title.x=element_text(size=18, vjust=-0.3), axis.title.y=element_text(size=18, vjust=0.3)) +
  theme(text = element_text(size=30))


### GRAPHS WITH AGE RESIDUALIZED SURFACE AREAS #### 


ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x=DrinksMonth.nona, y=bil.amygdala.a_resid, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm", size = 2) +
  labs(x = "DRINKS PER MONTH", y = "BILATERIAL AMYGDALAR SURFACE AREA") + theme(text = element_text(size=30))

ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], aes(x=DrinksMonth.nona, y=bil.hippocampus.a_resid, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm", size = 2) +
  labs(x = "DRINKS PER MONTH", y = "BILATERIAL HIPPCAMPAL SURFACE AREA") + theme(text = element_text(size=30))

### For hippocampus
SUQ.82_maget$tmp <- hippocampus_vertex.sa_l[1071+1,]
SUQ.82_maget$peak.hippocampus.l.a_resid <- with(SUQ.82_maget, lm(tmp~ AGESCAN + SEX)$residuals) 
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], aes(x=DrinksMonth.nona, y=peak.hippocampus.l.a_resid, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm", size = 3, lty = 1) +
  labs(x = "DRINKS PER MONTH", y = "LEFT HIPPOCAMPAL VOLUME AT VERTEX 1071") + theme(text = element_text(size=30))

SUQ.82_maget$tmp <- hippocampus_vertex.sa_r[1159+1,]
SUQ.82_maget$peak.hippocampus.r.a_resid <- with(SUQ.82_maget, lm(tmp~ AGESCAN + SEX)$residuals) 
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], aes(x=DrinksMonth.nona, y=peak.hippocampus.r.a_resid, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm", size = 3) +
  labs(x = "DRINKS PER MONTH", y = "RIGHT HIPPOCAMPAL VOLUME AT VERTEX 1159") + theme(text = element_text(size=30))

### For Amygdala
SUQ.82_maget$tmp <- amygdala_vertex.sa_l[59+1,]
SUQ.82_maget$peak.amygdala.l.a_resid <- with(SUQ.82_maget, lm(tmp~ AGESCAN + SEX)$residuals) 
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], aes(x=DrinksMonth.nona, y=peak.amygdala.l.a_resid, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm", size = 3) +
  labs(x = "DRINKS PER MONTH", y = "LEFT AMYGDALA VOLUME AT VERTEX 59") + theme(text = element_text(size=30))

SUQ.82_maget$tmp <- amygdala_vertex.sa_r[472+1,]
SUQ.82_maget$peak.amygdala.r.a_resid <- with(SUQ.82_maget, lm(tmp~ AGESCAN + SEX)$residuals) 
ggplot(SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], aes(x=DrinksMonth.nona, y=peak.amygdala.r.a_resid, group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm", size = 3) +
  labs(x = "DRINKS PER MONTH", y = "RIGHT AMYGDALA VOLUME AT VERTEX 472") + theme(text = element_text(size=30))



## UPDATED GRAPHS WITH AGE RESIDUALIZED SURFACE AREAS - WITH FITS +/- OUTLIERS REMOVED - AREA

ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.amygdala.a_resid, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.amygdala.a_resid, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=bil.amygdala.a_resid, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = .5, se=F) +
  labs(x = "DRINKS PER MONTH", y = "BILATERIAL AMYGDALAR SURFACE AREA") + theme(text = element_text(size=24))

ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.hippocampus.a_resid, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.hippocampus.a_resid, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=bil.hippocampus.a_resid, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = .5, se=F) +
  labs(x = "DRINKS PER MONTH", y = "BILATERIAL HIPPOCAMPAL SURFACE AREA") + theme(text = element_text(size=24))



## UPDATED GRAPHS WITH AGE RESIDUALIZED SURFACE AREAS - WITH FITS +/- OUTLIERS REMOVED - VOLUME

ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.amygdala.v_resid, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.amygdala.v_resid, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=bil.amygdala.v_resid, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = 1, se=F) +
  labs(x = "PEAK DRINKS PER MONTH IN ADULTHOOD", y = "AMYGDALA VOLUME IN CHILDHOOD") + theme(text = element_text(size=20))

ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.hipp.v_resid, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.hipp.v_resid, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=bil.hipp.v_resid, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = 1, se=F) +
  labs(x = "PEAK DRINKS PER MONTH IN ADULTHOOD", y = "HIPPOCAMPUS VOLUME IN CHILDHOOD") + theme(text = element_text(size=20))

## GRAPHS WITH AGE, SEX, TOTAL BRAIN VOLUME RESIDUALIZED

# hippocampus

SUQ.82_maget_temp <- SUQ.82_maget[2:82,]

# total brain volume residualized
s_temp <- lm(data = SUQ.82_maget_temp[SUQ.82_maget_temp$DrinksMonth.nona < 300,], bil.hipp.v~age.cent + GMV_WMV + SEX)
ggplot(SUQ.82_maget_temp[SUQ.82_maget_temp$DrinksMonth.nona < 300,], aes(x=DrinksMonth.nona, y=resid(s_temp), group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "bil.hippocampus.a_age_sex_tbv_resid") + theme(text = element_text(size=15))

# total brain volume residualized + outliers excluded
s_temp <- lm(data = SUQ.82_maget_temp[SUQ.82_maget_temp$DrinksMonth.nona < 150,], bil.hipp.v~age.cent + GMV_WMV + SEX)
ggplot(SUQ.82_maget_temp[SUQ.82_maget_temp$DrinksMonth.nona < 150,], aes(x=DrinksMonth.nona, y=resid(s_temp), group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "bil.hippocampus.a_age_sex_tbv_resid") + theme(text = element_text(size=15))

ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.hipp.v_resid, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=bil.hipp.v_resid, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=bil.hipp.v_resid, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = 1, se=F) +
  labs(x = "PEAK DRINKS PER MONTH IN ADULTHOOD", y = "HIPPOCAMPUS VOLUME IN CHILDHOOD") + theme(text = element_text(size=20))

# amygdala

# total brain volume residualized
s_temp <- lm(data = SUQ.82_maget_temp[SUQ.82_maget_temp$DrinksMonth.nona < 300,], bil.amygdala.v~age.cent + GMV_WMV + SEX)
ggplot(SUQ.82_maget_temp[SUQ.82_maget_temp$DrinksMonth.nona < 300,], aes(x=DrinksMonth.nona, y=resid(s_temp), group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "bil.hippocampus.a_age_sex_tbv_resid") + theme(text = element_text(size=15))

# total brain volume residualized + outliers excluded
s_temp <- lm(data = SUQ.82_maget_temp[SUQ.82_maget_temp$DrinksMonth.nona < 150,], bil.amygdala.v~age.cent + GMV_WMV + SEX)
ggplot(SUQ.82_maget_temp[SUQ.82_maget_temp$DrinksMonth.nona < 150,], aes(x=DrinksMonth.nona, y=resid(s_temp), group=SEX, color=SEX)) + geom_point() + geom_smooth(method="lm") +
  labs(x = "DRINKS PER MONTH", y = "bil.hippocampus.a_age_sex_tbv_resid") + theme(text = element_text(size=15))






#### AMYGDALA, HIPPOCAMPUS, STRIATUM VOLUME/SA CONTROLLING FOR TBV - ADDED MARCH 2021 ####
#### SUPPLEMENTARY FIGURE 1


## UPDATED GRAPHS WITH AGE RESIDUALIZED SURFACE AREAS - WITH FITS +/- OUTLIERS REMOVED - AREA

y_temp <- lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.sa~age.cent + GMV_WMV, SEX)
ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = .5, se=F) +
  labs(x = "DRINKS PER MONTH", y = "BILATERIAL AMYGDALAR SURFACE AREA") + theme(text = element_text(size=24))

y_temp <- lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hippocampus.a~age.cent + GMV_WMV, SEX)
ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = .5, se=F) +
  labs(x = "DRINKS PER MONTH", y = "BILATERIAL HIPPOCAMPAL SURFACE AREA") + theme(text = element_text(size=24))

y_temp <- lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.a~age.cent + GMV_WMV, SEX)
ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = .5, se=F) +
  labs(x = "DRINKS PER MONTH", y = "BILATERIAL STRIATUM SURFACE AREA") + theme(text = element_text(size=24))

## UPDATED GRAPHS WITH AGE RESIDUALIZED SURFACE AREAS - WITH FITS +/- OUTLIERS REMOVED - VOLUME

y_temp <- lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.amygdala.v~age.cent + GMV_WMV, SEX)
ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = 1, se=F) +
  labs(x = "PEAK DRINKS PER MONTH IN ADULTHOOD", y = "AMYGDALA VOLUME IN CHILDHOOD") + theme(text = element_text(size=20))

y_temp <- lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.hipp.v~age.cent + GMV_WMV, SEX)
ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = 1, se=F) +
  labs(x = "PEAK DRINKS PER MONTH IN ADULTHOOD", y = "HIPPOCAMPUS VOLUME IN CHILDHOOD") + theme(text = element_text(size=20))

y_temp <- lm(data = SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 150,], bil.striatum.v~age.cent + GMV_WMV, SEX)
ggplot() +
  geom_point(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX)) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 300,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX, color=SEX), method="lm", size = 2) +
  geom_smooth(data=SUQ.82_maget[SUQ.82_maget$DrinksMonth.nona < 125,], mapping=aes(x=DrinksMonth.nona, y=y_temp, group=SEX), color="black", method="lm", alpha=0.5, lty=2, size = 1, se=F) +
  labs(x = "PEAK DRINKS PER MONTH IN ADULTHOOD", y = "STRIATUM VOLUME IN CHILDHOOD") + theme(text = element_text(size=20))

