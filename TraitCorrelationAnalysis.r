############### NOTES ###################
# This code utilizes qualitative growth data across yeast species. The growth data was originally designated by letters/symbols. We changed the letters/symbols to numbers in the following way: + = 1; - = 0; n = NA; w = 1, and v = 0.5. The "v" stands for variable growth and was either tested across multiple strains or tested multiple times for the same strain and growth was found to be both weak and also not occurring. We are treating these cases as a different class from growth to be conservative.

options(stringsAsFactors = FALSE)

library("grofit")
library("doBy")
library("RColorBrewer")
library("pvclust")
library("ggplot2")
library("pvclust")
library("ape")
library("pheatmap")
library("dendextend")
library("igraph")

#################################################
############## Load Data ########################
#################################################

GrowthQual_m <- read.csv("YeastBookGrowth_df.csv", header = TRUE, check.names = FALSE) # Qualitative growth data that have been modified from letters/symbols to numbers
Treatment <- read.csv("Treatments_2015-06-15.csv", header = TRUE, check.names = FALSE) # Data table of different treatments being analyzed
GeneraList = read.csv("GeneraList.csv", header = TRUE)
####################################################
############## Functions ###########################
####################################################

######### Standard error function #########
std <- function(x) sd(x, na.rm = TRUE)/sqrt(length(x)) # Function to calculate standard error

######### Number of traits tested #########
tested <- function(d){
  output = vector()
  for(i in 1:nrow(d)){
    output[i] = length(which(!is.na(d[i,4:ncol(d)])))
  }
  return(output)
}

#####################################################
############## Data table modification ##############
#####################################################
# Remove anything that is not a Ascomycetous fungus
GeneraDrop = GeneraList[which(GeneraList[,2] == "B" | GeneraList[,2] == "NY"),1]

GrowthQual_m = GrowthQual_m[-which(GrowthQual_m[,"Genus"] %in% GeneraDrop),]

TreatNAs = data.frame(Treatment = character(), Number = numeric()) # Creates an empty dataframe to quantify the number of times a treatment was NA for the yeast

TREATMENT = 1 # Designates number for column to avoid string matching
NUMBER = 2 # Designates number for column to avoid string matching

for(i in 1:ncol(GrowthQual_m)){ # Start a loop to calculate the number NAs for a treatment
  TreatNAs[i, TREATMENT] = names(GrowthQual_m)[i] # Adds treatment name to row i of dataframe
  TreatNAs[i, NUMBER] = length(which(is.na(GrowthQual_m[,i]) == TRUE)) # Quantifies and adds the number of NAs to row i of the dataframe
} # End loop

AverageTests = mean(TreatNAs[,2])
KeepTreats = TreatNAs[which(TreatNAs[,"Number"] < 155),"Treatment"] # Determines which treatments have less than 155 NAs

GrowthQual_m = GrowthQual_m[,which(colnames(GrowthQual_m) %in% KeepTreats)] # Keeps treatments with less than 155 NAs
write.csv(GrowthQual_m, file = "GrowthQual_m.csv")


#YeastNAs = data.frame(Yeast = character(), Number = numeric()) # Creates an empty dataframe to quantify the number of NAs for a yeast

#YEAST = 1 # Designates number for column to avoid string matching
#NUMBER = 2 # Designates number for column to avoid string matching

#for(i in 1:nrow(GrowthQual_m)){ # Start a loop to calculate the number of NAs for a yeast
#  YeastNAs[i, YEAST] = GrowthQual_m[i,"Species"] # Adds species name to row i of dataframe
#  YeastNAs[i, NUMBER] = length(which(is.na(GrowthQual_m[i,]) == TRUE)) # Quantifies and adds the number of NAs for a yeast to row i of the dataframe
#} # End loop

#AverageSpecies = mean(YeastNAs[,2])
#KeepYeast = YeastNAs[which(YeastNAs[,"Number"] < 16),"Yeast"] # Determines which yeasts have less than 16 NAs

#GrowthQual_m = GrowthQual_m[which(GrowthQual_m[,"Species"] %in% KeepYeast),] # Keeps yeasts with less than 16 NAs

######### Create vectors of different treatments #########

Treatment_split = splitBy("Identifier", Treatment)
Fermentation_T = as.vector(Treatment_split[[which(names(Treatment_split) == "F")]][,1]) # Split out fermentation environments and create a vector
Carbon_T = as.vector(Treatment_split[[which(names(Treatment_split) == "C")]][,1]) # Split out carbon environments and creates a list
Nitrogen_T = as.vector(Treatment_split[[which(names(Treatment_split) == "N")]][,1]) # Splits out nitrogen environments and creates a list
Temperature_T = as.vector(Treatment_split[[which(names(Treatment_split) == "T")]][,1]) # Splits out temperature environments and creates a list

######### Treatment and Info column designation #########
TESTCOLSTART <- 3 # Insert the column that actual data starts here
INFOCOLENDS <- 2 # Insert the column that info data ends at

######### Data table of information columns in CBSQual_m #########
InfoColumns <- GrowthQual_m[,c(1:INFOCOLENDS)] # Columns that contain strain information
TestColumns <- colnames(GrowthQual_m[TESTCOLSTART:ncol(GrowthQual_m)]) # Data columns

################################################
################ Trait Analyses ################
################################################

########### Quantify the number of environments a a species can grow in ###########

TotalGrowthSpecies_df <- data.frame(Species = character(), RawGrowth = numeric(), TotalTested = numeric(), GrowthRatio = numeric()) # Creates empty dataframe to quantify growth in

# Designates columns for loop - Removes string matching
SPECIES = 1
RAWGROWTH = 2
TOTALTESTED = 3
GROWTHRATIO = 4

for(i in 1:nrow(GrowthQual_m)){ #Starts loop to quantify growth
  Values = as.numeric(GrowthQual_m[i,TESTCOLSTART:ncol(GrowthQual_m)])
  Values[which(Values == 0.5)] = 1
  Growth = sum(Values, na.rm = TRUE) # Calculates raw growth for a species(sum of 1's and 0.5 in the column)
  Tested = length(which(is.na(GrowthQual_m[i,TESTCOLSTART:ncol(GrowthQual_m)]) == FALSE)) # Calculates the number of environments tested for the species (removes NAs)
  Ratio = Growth/Tested
  TotalGrowthSpecies_df[i,SPECIES] = GrowthQual_m[i, "Species"] # Adds the species name to row i of the dataframe
  TotalGrowthSpecies_df[i, RAWGROWTH] = Growth # Adds the number of environments the species could grow in to row i of the dataframe
  TotalGrowthSpecies_df[i, TOTALTESTED] = Tested # Adds the number of environments that were tested for the species to row i of the dataframe
  TotalGrowthSpecies_df[i, GROWTHRATIO] = Ratio # Calculates the ratio of growth environments a species could grow in to row i of the dataframe
} # Ends loop

# Creates an ordered scatter plot of growth ratios across species
#SpeciesGrowthPlot = ggplot(TotalGrowthSpecies_df, aes(x = reorder(TotalGrowthSpecies_df$Species, TotalGrowthSpecies_df$GrowthRatio, order = is.ordered(TotalGrowthSpecies_df$Species)), y = GrowthRatio))+geom_point(stat = "identity")
#SpeciesGrowthPlot = SpeciesGrowthPlot +theme_bw()+xlab("Species")+ylab("Growth Ratio")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 5), strip.text.y = element_text(size = 12))
#ggsave(paste("SpeciesGrowthPlot", Sys.Date(), ".pdf", sep = ""))

# Creates a histogram of growth ratios across species
a = ggplot(TotalGrowthSpecies_df, aes(GrowthRatio))+geom_histogram(binwidth = 0.05, fill = "gray", colour = "black")
a = a + theme_bw()+xlab("Growth Ratio")+ylab("Frequency")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 8), strip.text.y = element_text(size = 12))
ggsave(paste("SpeciesGrowthBreadth_histogram_", Sys.Date(), ".pdf", sep = ""))

########### Quantify the number of fermentation a a species can grow in ###########

#FermentGrowthSpecies_df <- data.frame(Species = character(), RawGrowth = numeric(), TotalTested = numeric(), GrowthRatio = numeric()) # Creates empty dataframe to quantify growth in different fermentation environments

# Designates columns for loop - Removes string matching
#SPECIES = 1
#RAWGROWTH = 2
#TOTALTESTED = 3
#GROWTHRATIO = 4

#FermentData = GrowthQual_m[,which(colnames(GrowthQual_m) %in% Fermentation_T)] # Creates a dataframe of growth data for all fermentation treatments

#for(i in 1:nrow(FermentData)){ #Starts loop to quantify growth
#   Growth = sum(FermentData[i,], na.rm = TRUE) # Calculates raw growth for a species(sum of 1's and 0.5 in the column)
#   Tested = length(which(is.na(FermentData[i,]) == FALSE)) # Calculates the number of fermentation environments tested for the species (removes NAs)
#   Ratio = Growth/Tested
#   FermentGrowthSpecies_df[i,SPECIES] = GrowthQual_m[i,"Species"] # Adds the species name to row i of the dataframe
#   FermentGrowthSpecies_df[i, RAWGROWTH] = Growth # Adds the number of fermentation environments the species could grow in to row i of the dataframe
#   FermentGrowthSpecies_df[i, TOTALTESTED] = Tested # Adds the number of fermentation environments that were tested for the species to row i of the dataframe
#   FermentGrowthSpecies_df[i, GROWTHRATIO] = Ratio # Calculates the ratio of fermentation growth environments a species could grow in to row i of the dataframe
# } # Ends loop

## Creates an ordered scatter plot of fermentation growth ratios across species
#a = ggplot(FermentGrowthSpecies_df, aes(x = reorder(FermentGrowthSpecies_df$Species, FermentGrowthSpecies_df$GrowthRatio, order = is.ordered(FermentGrowthSpecies_df$Species)), y = GrowthRatio))+geom_point(stat = "identity")
#a = a +theme_bw()+xlab("Species")+ylab("Fermentation Growth Ratio")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 5), strip.text.y = element_text(size = 12))
#ggsave(paste("Fermentation_SpeciesGrowthPlot", Sys.Date(), ".pdf", sep = ""))

## Creates a histogram of fermentation growth ratios across species
#a = ggplot(FermentGrowthSpecies_df, aes(GrowthRatio))+geom_histogram(binwidth = 0.1, fill = "gray", colour = "black")
#a = a + theme_bw()+xlab("Growth Ratio")+ylab("Frequency")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 8), strip.text.y = element_text(size = 12))
#ggsave(paste("Fermentation_SpeciesGrowth_histogram_", Sys.Date(), ".pdf", sep = ""))

############ Quantify the number of carbon a a species can grow in ###########

#CarbonGrowthSpecies_df <- data.frame(Species = character(), RawGrowth = numeric(), TotalTested = numeric(), GrowthRatio = numeric()) # Creates empty dataframe to quantify carbon growth in

## Designates columns for loop - Removes string matching
#SPECIES = 1
#RAWGROWTH = 2
#TOTALTESTED = 3
#GROWTHRATIO = 4

#CarbonData = GrowthQual_m[which(colnames(GrowthQual_m) %in% Carbon_T),] # Creates a dataframe of growth data for all carbon treatments


#for(i in 1:nrow(CarbonData)){ #Starts loop to quantify growth
#  Growth = sum(CarbonData[i,], na.rm = TRUE) # Calculates raw carbon growth for a species(sum of 1's and 0.5 in the column)
#  Tested = length(which(is.na(CarbonData[i]) == FALSE)) # Calculates the number of carbon environments tested for the species (removes NAs)
#  Ratio = Growth/Tested
#  CarbonGrowthSpecies_df[i,SPECIES] = colnames(CarbonData)[i] # Adds the species name to row i of the dataframe
#  CarbonGrowthSpecies_df[i, RAWGROWTH] = Growth # Adds the number of carbon environments the species could grow in to row i of the dataframe
#  CarbonGrowthSpecies_df[i, TOTALTESTED] = Tested # Adds the number of carbon environments that were tested for the species to row i of the dataframe
#  CarbonGrowthSpecies_df[i, GROWTHRATIO] = Ratio # Calculates the ratio of carbon growth environments a species could grow in to row i of the dataframe
#} # Ends loop

# # Creates an ordered scatter plot of carbon growth ratios across species
# a = ggplot(CarbonGrowthSpecies_df, aes(x = reorder(CarbonGrowthSpecies_df$Species, CarbonGrowthSpecies_df$GrowthRatio, order = is.ordered(CarbonGrowthSpecies_df$Species)), y = GrowthRatio))+geom_point(stat = "identity")
# a = a +theme_bw()+xlab("Species")+ylab("Carbon Growth Ratio")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 5), strip.text.y = element_text(size = 12))
# ggsave(paste("Carbon_SpeciesGrowthPlot", Sys.Date(), ".pdf", sep = ""))
# 
# # Creates a histogram of carbon growth ratios across species
# a = ggplot(CarbonGrowthSpecies_df, aes(GrowthRatio))+geom_histogram(binwidth = 0.1, fill = "gray", colour = "black")
# a = a + theme_bw()+xlab("Carbon Growth Ratio")+ylab("Frequency")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 8), strip.text.y = element_text(size = 12))
# ggsave(paste("Carbon_SpeciesGrowth_histogram_", Sys.Date(), ".pdf", sep = ""))
# 
# ########### Quantify the number of nitrogen environments a species can grow in ###########
# 
# NitrogenGrowthSpecies_df <- data.frame(Species = character(), RawGrowth = numeric(), TotalTested = numeric(), GrowthRatio = numeric()) # Creates empty dataframe to quantify nitrogen growth in
# 
# # Designates columns for loop - Removes string matching
# SPECIES = 1
# RAWGROWTH = 2
# TOTALTESTED = 3
# GROWTHRATIO = 4
# 
# NitrogenData = GrowthQualSpecies[which(rownames(GrowthQualSpecies) %in% Nitrogen_T),] # Creates a dataframe of growth data for all nitrogen treatments
# 
# 
# for(i in 1:ncol(NitrogenData)){ #Starts loop to quantify nitrogen growth
#   Growth = sum(NitrogenData[i], na.rm = TRUE) # Calculates raw nitrogen growth for a species(sum of 1's and 0.5 in the column)
#   Tested = length(which(is.na(NitrogenData[i]) == FALSE)) # Calculates the number of nitrogen environments tested for the species (removes NAs)
#   Ratio = Growth/Tested
#   NitrogenGrowthSpecies_df[i,SPECIES] = colnames(NitrogenData)[i] # Adds the species name to row i of the dataframe
#   NitrogenGrowthSpecies_df[i, RAWGROWTH] = Growth # Adds the number of nitrogen environments the species could grow in to row i of the dataframe
#   NitrogenGrowthSpecies_df[i, TOTALTESTED] = Tested # Adds the number of nitrogen environments that were tested for the species to row i of the dataframe
#   NitrogenGrowthSpecies_df[i, GROWTHRATIO] = Ratio # Calculates the ratio of nitrogen growth environments a species could grow in to row i of the dataframe
# } # Ends loop
# 
# # Creates an ordered scatter plot of nitrogen growth ratios across species
# a = ggplot(NitrogenGrowthSpecies_df, aes(x = reorder(NitrogenGrowthSpecies_df$Species, NitrogenGrowthSpecies_df$GrowthRatio, order = is.ordered(NitrogenGrowthSpecies_df$Species)), y = GrowthRatio))+geom_point(stat = "identity")
# a = a +theme_bw()+xlab("Species")+ylab("Nitrogen Growth Ratio")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 5), strip.text.y = element_text(size = 12))
# ggsave(paste("Nitrogen_SpeciesGrowthPlot", Sys.Date(), ".pdf", sep = ""))
# 
# # Creates a histogram of nitrogen growth ratios across species
# a = ggplot(NitrogenGrowthSpecies_df, aes(GrowthRatio))+geom_histogram(binwidth = 0.1, fill = "gray", colour = "black")
# a = a + theme_bw()+xlab("Nitrogen Growth Ratio")+ylab("Frequency")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 8), strip.text.y = element_text(size = 12))
# ggsave(paste("Nitrogen_SpeciesGrowth_histogram_", Sys.Date(), ".pdf", sep = ""))
# 
# ########### Quantify the number of temperature a a species can grow in ###########
# 
# TemperatureGrowthSpecies_df <- data.frame(Species = character(), RawGrowth = numeric(), TotalTested = numeric(), GrowthRatio = numeric()) # Creates empty dataframe to quantify temperature growth in
# 
# # Designates columns for loop - Removes string matching
# SPECIES = 1
# RAWGROWTH = 2
# TOTALTESTED = 3
# GROWTHRATIO = 4
# 
# TemperatureData = GrowthQualSpecies[which(rownames(GrowthQualSpecies) %in% Temperature_T),] # Creates a dataframe of growth data for all temperature treatments
# 
# 
# for(i in 1:ncol(TemperatureData)){ #Starts loop to quantify temperature growth
#   Growth = sum(TemperatureData[i], na.rm = TRUE) # Calculates raw temperature growth for a species(sum of 1's and 0.5 in the column)
#   Tested = length(which(is.na(TemperatureData[i]) == FALSE)) # Calculates the number of environments tested for the species (removes NAs)
#   Ratio = Growth/Tested
#   TemperatureGrowthSpecies_df[i,SPECIES] = colnames(TemperatureData)[i] # Adds the species name to row i of the dataframe
#   TemperatureGrowthSpecies_df[i, RAWGROWTH] = Growth # Adds the number of temperature environments the species could grow in to row i of the dataframe
#   TemperatureGrowthSpecies_df[i, TOTALTESTED] = Tested # Adds the number of temperature environments that were tested for the species to row i of the dataframe
#   TemperatureGrowthSpecies_df[i, GROWTHRATIO] = Ratio # Calculates the ratio of temperature growth environments a species could grow in to row i of the dataframe
# } # Ends loop
# 
# # Creates an ordered scatter plot of temperature growth ratios across species
# a = ggplot(TemperatureGrowthSpecies_df, aes(x = reorder(TemperatureGrowthSpecies_df$Species, TemperatureGrowthSpecies_df$GrowthRatio, order = is.ordered(TemperatureGrowthSpecies_df$Species)), y = GrowthRatio))+geom_point(stat = "identity")
# a = a +theme_bw()+xlab("Species")+ylab("Temperature Growth Ratio")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 5), strip.text.y = element_text(size = 12))
# ggsave(paste("Temperature_SpeciesGrowthPlot", Sys.Date(), ".pdf", sep = ""))
# 
# # Creates a histogram of temperature growth ratios across species
# a = ggplot(TemperatureGrowthSpecies_df, aes(GrowthRatio))+geom_histogram(binwidth = 0.1, fill = "gray", colour = "black")
# a = a + theme_bw()+xlab("Temperature Growth Ratio")+ylab("Frequency")+theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 8), strip.text.y = element_text(size = 12))
# ggsave(paste("Temperature_SpeciesGrowth_histogram_", Sys.Date(), ".pdf", sep = ""))
# 
# ########### Combined treatment growth counts ###########
# TreatmentGrowthSpecies_df = cbind(FermentGrowthSpecies_df, CarbonGrowthSpecies_df, NitrogenGrowthSpecies_df, TemperatureGrowthSpecies_df) # Combines all growth treatment dataframes
# colnames(TreatmentGrowthSpecies_df) = c("Species", "Ferment_Raw", "Ferment_Tested", "Ferment_GrowthRatio", "Drop", "Carbon_Raw", "Carbon_Tested", "Carbon_GrowthRatio", "Drop", "Nitrogen_Raw", "Nitrogen_Tested", "Nitrogen_GrowthRatio", "Drop", "Temperature_Raw", "Temperature_Tested", "Temperature_GrowthRatio") # Manually set column name of treatment growth dataframe
# TreatmentGrowthSpecies_df = TreatmentGrowthSpecies_df[,-which(colnames(TreatmentGrowthSpecies_df) == "Drop")] # Drops redundant species columns
# 
# ########### Correlations between environmental growth conditions ###########
# Carb_Ferment_lm = lm(Carbon_GrowthRatio ~ Ferment_GrowthRatio, TreatmentGrowthSpecies_df) # linear model of carbon and fermentation growth
# Temp_Ferment_lm = lm(Temperature_GrowthRatio ~ Ferment_GrowthRatio, TreatmentGrowthSpecies_df) # linear model of temperature and fermentation growth
# Nitrogen_Ferment_lm = lm(Nitrogen_GrowthRatio ~ Ferment_GrowthRatio, TreatmentGrowthSpecies_df) # linear model of nitrogen and fermentation growth
# Temp_Carbon_lm = lm(Temperature_GrowthRatio ~ Carbon_GrowthRatio, TreatmentGrowthSpecies_df) # linear model of carbon and temperature growth
# Nitrogen_Carbon_lm = lm(Nitrogen_GrowthRatio ~ Carbon_GrowthRatio, TreatmentGrowthSpecies_df) # linear model of carbon and nitrogen growth
# Nitrogen_Temp_lm = lm(Nitrogen_GrowthRatio ~ Temperature_GrowthRatio, TreatmentGrowthSpecies_df) # linear model of temperature and nitrogen growth
# 
########### Create a dendrogram based on traits for all species data ###############

GrowthQualSpecies_m = read.csv("GrowthQual_mSpecies.csv", header = TRUE, row.names = 1)
SpeciesDendro <- pvclust(GrowthQualSpecies_m, method.dist = "euclidean", method.hclust = "ward", nboot = 1000) # Creates a species tree based on trait presence

pdf(paste("Species_Dendrogram_", Sys.Date(), ".pdf", sep = ""), height = 25, width = 25) # Creates a pdf
SpeciesDend_results = as.dendrogram(SpeciesDendro)
SpeciesDend_results %>% pvclust_show_signif(SpeciesDendro, show_type = "lwd") %>% 
	circlize_dendrogram()
SpeciesDendro %>% text
SpeciesDend_results %>% pvrect()
dev.off() # Turns off PDF

SpeciesDendro.Clusters <- pvpick(SpeciesDendro) # Writes out significant clusters of species
write.csv(capture.output(SpeciesDendro.Clusters[[1]]), file = paste("SpeciesDendro_Clusters_", Sys.Date(), ".csv", sep = "")) # Saves a csv of all significant clusters

#pdf(paste("Species_Long_", Sys.Date(), ".pdf", sep = ""), height = 50, width = 200) # Creates a pdf
#SpeciesDend_results = as.dendrogram(SpeciesDendro)
#SpeciesDend_results %>% pvclust_show_signif(SpeciesDendro, show_type = "lwd") %>% 
#	plot()
#SpeciesDendro %>% text
#SpeciesDendro %>% pvrect()
#dev.off()

#library(dynamicTreeCut)
#hc <- GrowthQualSpecies %>% dist %>% hclust
#dend <- hc %>% as.dendrogram

#clusters <- cutreeDynamic(hc, distM = as.matrix(dist(GrowthQualSpecies)), method = "tree")
#clusters <- clusters[order.dendrogram(dend)]
#cluster_numbers <- unique(clusters) - (0 %in% clusters)
#n_clusters <- length(cluster_numbers)
#library("colorspace")
#cols <- rainbow_hcl(n_clusters)
#true_species_cols <- rainbow_hcl(10)[as.numeric(colnames(GrowthQualSpecies)[[order.dendrogram(dend)]))

########### Create a dendrogram based on traits for all genera data ###############
GrowthQualGeneraT_m = read.csv("GrowthQual_mGenusT.csv", header = TRUE, row.names = TRUE)

GeneraDendro <- pvclust(GrowthQualGeneraT_m, method.dist = "euclidean", method.hclust = "ward", nboot = 1000) # Creates a genera tree based on trait presence

pdf(paste("Genera_Dendrogram_", Sys.Date(), ".pdf", sep = ""), height = 25, width = 25) # Creates a pdf
GeneraDend_results = as.dendrogram(GeneraDendro)
GeneraDend_results %>% pvclust_show_signif(GeneraDendro, show_type = "lwd") %>% 
	circlize_dendrogram()
GeneraDendro %>% text
GeneraDend_results %>% pvrect()
dev.off() # Turns off PDF

GeneraDendro.Clusters <- pvpick(GeneraDendro) # Writes out significant clusters of genera
write.csv(capture.output(GeneraDendro.Clusters[[1]]), file = paste("GeneraDendro_Clusters_", Sys.Date(), ".csv", sep = "")) # Saves a csv of all signficant clusters

######### Calculates Trait Correlations #########
#GrowthQual_m = GrowthQual_m[,-j]

treatments = colnames(GrowthQual_m)[3:ncol(GrowthQual_m)]

CorrelationMatrix = matrix(0, nrow = length(treatments), ncol = length(treatments))
colnames(CorrelationMatrix) = treatments
rownames(CorrelationMatrix) = treatments

Correlation_df = data.frame(Treatment_A = character(), Treatment_B = character(), TValue = numeric(), Df = numeric(), CorrValue = numeric(), Pvalue = numeric())
k = 1
#GrowthQual_m = GrowthQual_m[,-j]
for(i in 3:length(GrowthQual_m)){
  for(j in 3:length(GrowthQual_m)){
    temp = cor.test(as.numeric(GrowthQual_m[,i]), as.numeric(GrowthQual_m[,j]))
    CorrelationMatrix[which(colnames(CorrelationMatrix) == colnames(GrowthQual_m)[i]),which(colnames(CorrelationMatrix) == colnames(GrowthQual_m)[j])] = temp$estimate[[1]]
    Correlation_df[k, 1] = colnames(GrowthQual_m)[i]
    Correlation_df[k, 2] = colnames(GrowthQual_m)[j]
    Correlation_df[k, 3] = temp[[1]][1]
    Correlation_df[k, 4] = temp[[2]][1]
    Correlation_df[k, 5] = temp[[4]][1]
    Correlation_df[k, 6] = temp[[3]][1]
    k = 1+k
  }
}

BHCorrection = p.adjust(Correlation_df[,6], method = "BH")
Correlation_df=cbind(Correlation_df, BHCorrection)

SigCorr_df = Correlation_df[which(Correlation_df$BHCorrection < 0.05),]
PositiveCorr = SigCorr_df[which(SigCorr_df$CorrValue > 0),c(1,2)]
NegativeCorr = SigCorr_df[which(SigCorr_df$CorrValue < 0),c(1,2)]
write.csv(PositiveCorr, file = "PosCorrelation_padj.csv", row.names = FALSE)
write.csv(NegativeCorr, file = "NegCorrelation_padj.csv", row.names = FALSE)

CorrelationMatrix= CorrelationMatrix[-which(rownames(CorrelationMatrix) == "Glucose"), -which(colnames(CorrelationMatrix) == "Glucose")]

#CORRMAT[lower.tri(CORRMAT)] = NA
COLORS = c("#8EA6CC","#325A99","white", "#EEFEFF","#FFC4AE", "#CCAEAA")
pdf("Correlation_hmap.pdf")
pheatmap(CorrelationMatrix, clustering_distance_rows = "euclidean", cluster_cols = FALSE, clustering_method = "ward", border_color = "black", cellheight = 8, cellwidth = 8,color = colorRampPalette(COLORS)(100))
dev.off()

pvClustMatrix = GrowthQual_m[,3:ncol(GrowthQual_m)]
Treatment_dendro = pvclust(pvClustMatrix, nboot = 1000, method.dist = "euclidean", method.hclust = "ward")

######### Calculate the number of correlated and anti-correlated traits per species #########

# Calculates the number correlated traits within a species
Species = GrowthQual_m$Species

NegCorrCounts_df = data.frame(Species = character(), NegCorrTrait = numeric())
for(i in 1:length(Species)){
  tempData = GrowthQual_m[which(GrowthQual_m$Species == Species[[i]]),]
  tempCount_list = list()
  for(j in 1:nrow(NegativeCorr)){
    tempCount_list[[j]] = list()
    TraitA = which(colnames(tempData) == NegativeCorr[j,1])
    TraitB = which(colnames(tempData) == NegativeCorr[j,2])
    temp_df=tempData[,c(TraitA, TraitB)]
    if(is.na(temp_df[,1]) == FALSE & is.na(temp_df[,2]) == FALSE){
      if(temp_df[,1] == temp_df[,2]){
        tempValue = temp_df[,1]
      }
      if(tempValue == 0.5){
        tempValue = 1
      }
      tempCount_list[[j]] = tempValue
    }
  }
  NegCorrCounts_df[i,1] = Species[[i]]
  VALUES = as.numeric(unlist(tempCount_list))
  NegCorrCounts_df[i,2] = sum(VALUES)/2
}

BreadthCounts_df = merge(TotalGrowthSpecies_df, NegCorrCounts_df, by = "Species")

# Creates a scatterplot of the number of correlated and anti-correlated traits per species
a = ggplot(BreadthCounts_df, aes(x = Correlated, y = AntCorrelated))+geom_point(shape = 19)+geom_smooth(method = lm, se = FALSE)
a = a+ theme_bw()+xlab("Correlation Trait Counts")+ylab("Anti-Correlated Trait Counts")+scale_fill_grey()+ theme(axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 8, angle = 90),axis.text.y = element_text(size = 8), strip.text.y = element_text(size = 12))
ggsave(paste("CorrAntCorrTraitCountScatter_", Sys.Date(), ".pdf", sep = "")) # Saves scatterplot of correlated and anti-correlated traits

######### Determine whether there is a trade-off between the number of anti-correlated traits and the number of traits a species can use #########

#############################################################
############## Trait Co-occurrence among Genera ##############
#############################################################
GrowthQualGenera_m = read.csv("GrowthQual_mGenus.csv", header = TRUE)
GeneraSplit = splitBy("X", GrowthQualGenera_m)

GeneraCorr = data.frame()
