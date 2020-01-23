# ----------------------------------------------------------------------------------------------------- #
# -------------- Brief tutorial & example script on some basic community analyses --------------------- # 
# ----------------------------------------------------------------------------------------------------- #
# Ecological analyses using metabarcoding mOTU assignments
# Uses data from an arthropod metabarcoding survey conducted over 3 consecutive years in the UK
# Aimed at assessing the impact of land cover & farm management on biodiversity & community composition

# Load packages for analyses & plotting 
required_pkgs <- c("tidyverse","ggridges","viridis","vegan","directlabels","reshape2","fso","MuMIn","lme4","rich")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(required_pkgs, character.only=TRUE)

# ------ Set up data for analyses

# Read in mOTU table 
all_samples <- readRDS("mOTU_samples.rds")
str(all_samples) # Data frame of 
# Sample = sample label ID 
# mOTU = unique mOTU identifier
# Domain ... Species = taxonomic asignments to species level
# Domain_pc ... Species_pc = LCA1 score for each assignment level
# Read_n = read number for that mOTU
# Run = run index (samples split over multiple runs)

# Read in ecological covariate data
lcm_data <- read.csv("LCM_2015.csv") %>% dplyr::select(1:16) %>% filter(Resolution==1) # CEH landcover map  - % cover of different land uses
cc_data  <- read.csv("CC_2016.csv") %>% dplyr::select(4:15) %>% filter(Resoultion==1) # CEH crop cover map  - % cover of different crops
cc_data$shannon <- diversity(cc_data[,-11],index="shannon") # Calculate a shannon DI for cc
lcm_data[,4:15] <- lcm_data[,4:15]/rowSums(lcm_data[,4:15]) # LCM as proportion 
lcm_data$prop_non <- 1-(lcm_data$Arable_horticulture+lcm_data$Improved_grassland) # Proportion of non-arable for LCM
all_meta <- cbind(lcm_data,cc_data) # Combine


# ------------------------------------------------------------------------------ #
# ------------------        Set up & filter data        ------------------------ # 
# ------------------------------------------------------------------------------ #

# Some basic data house keeping & indexing 
contam_filter <- all_samples %>% dplyr::filter(!str_detect(Sample, "UV" )) %>%  # Remove un-pooled validation samples 
                                     mutate(Farm   = str_sub(Sample, 1, 2), # Add in farm identifer
                                            Season = str_extract(Sample, pattern = "S[0-9]"), # Season 
                                            Trap   = str_extract(Sample, pattern ="T[0-9]"), # Pan trap number
                                            FS     = interaction(Farm, Season), # Farm * Season
                                            FT     = interaction(Farm, Trap), # Farm * Trap
                                            FTS    = interaction(Farm, Trap, Season)) %>%  # Farm * Trap * Season
                              dplyr::filter(!Farm %in% c("NA","EB","PC","TJ")) # Remove various extraction & sequencing blanks 

# ------------------------------------------------------------------------------ #
# ------------------      Filter out contamination     ------------------------ # 
# ------------------------------------------------------------------------------ #

# An index of organisms that appeared in extraction blanks, which we might want to exclude from the analyses
contam_index <- readRDS("contam_index.rds")
exclude_sp <- c("Chrysoperla nipponensis","Ectopsocus californicus") # Index of two non-uk species that appeared in samples - likely contamination from lab

# Filter out likely contaminants
all_filter<- left_join(contam_filter,contam_index,by="mOTU") %>%
                dplyr::select(-Species.y) 
                filter(Read_n > 30,   # Filter low read frequencies
                       Class_pc==100, # Specify that all mOTUs to be included in analyses should have high LCA1 scores up to classs
                       ! Species.x %in% exclude_sp ,  # Exclude non-uk/study site species
                       ! mOTU %in% contam_index$mOTU) # remove likely contaminant mOTUs

# --------------------------------------------------------------------------------- #
# ------------------      Filter out Singleton species     ------------------------ # 
# --------------------------------------------------------------------------------- #

# Common practice for many ordinations is to remove species that only occur once 
# Or sites that only have a single species. 
                
# Look at singleton Species & Sites with only one Species
singleton_species <- all_filter %>% group_by(mOTU) %>% filter(n()==1) %>% ungroup() # Singleton species
singleton_sites   <- all_filter %>% group_by(FTS)  %>% filter(n()==1) %>% ungroup() # Singleton sites

# Remove singleton species & Sites with only one species sequentially
no_singles_all <- all_filter %>% group_by(mOTU) %>% filter(n()>1) %>% ungroup() %>% 
                                 group_by(FTS)   %>% filter(n()>1) %>% ungroup() 


#-----------------------------------------------------------------------------#
#------------------ Non-metric multidimensional scaling ----------------------#
#-----------------------------------------------------------------------------#

# Ordination technique that allows visualisation of community compositions relative to one another
# Works well with absence presence data that MBC generates

# Setup ------------------
# Make an absence / presence matrix by Species
miniabs_pres <- dcast(unique(no_singles_all), FS ~ mOTU, length, value.var="mOTU") # Get a matrix of incidence (number of detections by mOTU)
miniabs_pres[,2:ncol(miniabs_pres)] <- 1*(miniabs_pres[,2:ncol(miniabs_pres)]>0) # Convert to absence / presence
miniabs_pres <- miniabs_pres %>% mutate(Farm=str_sub(FS,1,2)) # Quick restructure of farm level indexes
str(miniabs_pres) # DF where column = mOTU & row = sample

ncol(miniabs_pres) # Check number of species in analysis
no_qual <- miniabs_pres[,-ncol(miniabs_pres)] # Remove qualitative variables from data
ab_dist <- vegdist(no_qual[,-1],method="jaccard",binary=TRUE)  # Make a dissimilarity matrix binary for presence/absence
# might be worth looking at different method arguments depending on your data

# Analysis ------------------
nmds_1 <- metaMDS(ab_dist,k=5 ,wascores=FALSE, autotransform=FALSE,noshare=FALSE,trymax = 1000) # Run the NMDS 
# Number of dimensions , k , specifies how many dimensions to run the ordination in, the higher the number the more flexible the ordination, but becomes more difficult to interpret
# Sometimes the ordination will take a long time to converge, and it might be worth reinitialising with the `previous.best=` argument
nmds_2 <- metaMDS(ab_dist,k=5 ,wascores=FALSE, autotransform=FALSE,noshare=FALSE,trymax = 1000, previous.best = nmds_1) 

# Ordination diagnostics -----
# Stress is the major measure of pathological ordinations
# A good rule of thumb:
# >= 0.3 - Ordination is arbitrary tells you very little about community composition
# 0.2-0.3 - Ordination is suspect and should be interpreted with caution
# 0.1 - Fair fit
# 0.5 - Good fit

# Shepard plot
stressplot(nmds_1)
# Shows the original dissimilarity calculated in the object ab_dist to the  distances calculated in the ordination
# Larger scatter means a larger observed distance and a poorer fit
# This example does pretty well
# Should also look at a scree plot to determine the optimal stress for the NMDS - but can't remember right now how to get this...

# ------ Plot ordination against environmental variables
meta_lcm <- left_join(miniabs_pres,all_meta,by="Farm") # bind LCM & CC daata to abs pres matrix


# Retrieve data & format for plot 
NMDS1 <- nmds_1$points[,1] # First axis
NMDS2 <- nmds_1$points[,2] # Second axis
ggNMDS <- data.frame(NMDS1 = NMDS1, NMDS2 = NMDS2, 
                     FS   = miniabs_pres$FS) %>% 
  mutate(Season=str_sub(FS,4,5),
         Farm=str_sub(FS,1,2)) %>% left_join(meta_lcm,by="Farm")

# Plot ordination of community composition in each field
# Points represent a single farm in a single year, the closer the points the more similar communities are.
ggplot(ggNMDS,aes(NMDS1,NMDS2,colour=Region))+
  geom_point(size=3)+
  stat_ellipse()+ # Plots an ellipse for each sample region - clearly no discernable difference between UK counties
  scale_colour_viridis(discrete = T,end=0.8)+
  theme_classic()
  

# These next three lines calculate contours for three environmental variables w/ respect to ordination space
ordi.prop <- ordisurf(nmds_1,meta_lcm$prop_non,plot = FALSE,bs="ds",elect = FALSE, method = "GCV.Cp") # Proportion of non-arable surrounding a farm
ordi.blw <- ordisurf(nmds_1,meta_lcm$Broadleaf_woodland,plot = FALSE,bs="ds",elect = FALSE, method = "GCV.Cp")  # Proportion of Broadleaf woodland
ordi.shannon <- ordisurf(nmds_1,meta_lcm$shannon,plot = FALSE,bs="ds",elect = FALSE, method = "GCV.Cp") # Diversity index of crop cover

# Retrieve some summary information to display above the plots
shan.sum <- summary(ordi.shannon)
l_shannon <- paste(expression(R^{2}~(adj)),"==",round(shan.sum$r.sq,2)) # Rsquare
l_shannon2 <- paste("Deviance explained","=",round(shan.sum$dev.expl*100,2),"%") # Deviance explained
l_shannon3 <- paste(expression(P),"==",round(shan.sum$s.pv,2)) # P value

blw.sum <- summary(ordi.blw)
l_blw <- paste(expression(R^{2}~(adj)),"==",round(blw.sum$r.sq,2))
l_blw2 <- paste("Deviance explained","=",round(blw.sum$dev.expl*100,2),"%")
l_blw3 <- paste(expression(P),"==",round(blw.sum$s.pv,2))


prop.sum <- summary(ordi.prop)
l_prop <- paste(expression(R^{2}~(adj)),"==",round(prop.sum$r.sq,2))
l_prop2 <- paste("Deviance explained","=",round(prop.sum$dev.expl*100,2),"%")
l_prop3 <- paste(expression(P),"==",round(prop.sum$s.pv,2))

# Function to retrieve the contours to allow easy plotting --------------- 
retrieve_contours <- function(ordi){
## Retrieve contour info ## 
ordi.grid <- ordi$grid #extracts the ordisurf object
ordi.mite <- expand.grid(x = ordi.grid$x, y = ordi.grid$y) #get x and ys
ordi.mite$z <- as.vector(ordi.grid$z) #unravel the matrix for the z scores
ordi.mite.na <- data.frame(na.omit(ordi.mite)) #gets rid of the nas
}
# Function to retrieve the contours to allow easy plotting --------------- 

ordi.plot.blw <- retrieve_contours(ordi.blw)
ordi.plot.prop <- retrieve_contours(ordi.prop) 
ordi.plot.shannon <- retrieve_contours(ordi.shannon)



# Plot   
p.blw <-   ggplot()+
  stat_contour(data = ordi.plot.blw, aes(x = x, y = y, z = z,colour = ..level..),size=2, alpha=0.4)+
  scale_colour_viridis(option="D",begin = 0,end=0.8)+
  annotate("text",label= l_blw,x=-0.4,y=0.35,parse=T,hjust=0)+
  annotate("text",label= l_blw2,x=-0.4,y=0.3,hjust=0)+
  #annotate("text",label= l_blw3,x=-0.4,y=0.25,parse=T,size=8,hjust=0)+
  geom_point(data=ggNMDS,aes(NMDS1,NMDS2),pch=21,size=2,colour="black",fill="black", alpha=0.5)+
  labs(x="NMDS1",y="NMDS2")+
  theme_classic()

p.prop <-   ggplot()+
  stat_contour(data = ordi.plot.prop, aes(x = x, y = y, z = z,colour = ..level..),size=2, alpha=0.4)+
  scale_colour_viridis(option="D",begin = 0,end=0.8)+
  annotate("text",label= l_prop,x=-0.4,y=0.35,parse=T,hjust=0)+
  annotate("text",label= l_prop2,x=-0.4,y=0.30,hjust=0)+
  #annotate("text",label= l_prop3,x=-0.4,y=0.25,parse=T,size=8,hjust=0)+ 
  geom_point(data=ggNMDS,aes(NMDS1,NMDS2),pch=21,size=2,colour="black",fill="black", alpha=0.5)+
  labs(x="NMDS1",y="NMDS2")+
  theme_classic()

p.shannon <-   ggplot()+
  stat_contour(data = ordi.plot.shannon, aes(x = x, y = y, z = z,colour = ..level..),size=2, alpha=0.4)+
  scale_colour_viridis(option="D",begin = 0,end=0.8)+
  annotate("text",label= l_shannon,x=-0.4,y=0.35,parse=T,hjust=0)+
  annotate("text",label= l_shannon2,x=-0.4,y=0.3,hjust=0)+
  #annotate("text",label= l_shannon3,x=-0.4,y=0.25,parse=T,size=8,hjust=0)+
  geom_point(data=ggNMDS,aes(NMDS1,NMDS2),pch=21,size=2,colour="black",fill="black", alpha=0.5)+
  labs(x="NMDS1",y="NMDS2")+
  theme_classic()



# These plots display the ordination points across an environmental gradient (alongside metrics of model fit)
direct.label(p.blw,    list("top.points",cex=1,col="black",fill="black"))
direct.label(p.prop,   list("top.points",cex=1,col="black",fill="black"))
direct.label(p.shannon, list("top.points",cex=1,col="black",fill="black"))

# Low Rsquared and deviance values suggest environmental variables aren't driving community composition. 
