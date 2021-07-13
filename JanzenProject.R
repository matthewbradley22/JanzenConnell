#Working on pathogen tradeoff
#Load packages
library(tidyverse) 
library(gridExtra) 


### Initialize params ###

# Set up a grid to represent land

#Place individuals on plot
par(mfrow = c(1, 1))
xsize = ysize = 100

numInd = 100
numSpecies = 2
data <-  tibble("ID" = 1:numInd, "age" = 3, 
                species = sample(numSpecies, numInd, replace = TRUE), "xlocation" = runif(numInd, 0, 100), 
                "ylocation" = runif(numInd,0,100), "parentalDistance" = 0, "Pathogen" = 0, "Resistance" = 0,
                "pathStrength" = 0)

#Choose an initial 25% of population to be infected
infected <- sample(numInd, numInd/4)
data$Pathogen[which(data$ID %in% infected)] <- 1


plot(data$xlocation, data$ylocation, col = data$species, pch = 20, xlim = c(1,100), ylim = c(1,100),
     xlab = "X-Coord", ylab = "Y-Coord", main = "Initial Trees")


#Set up rates of birth/death/growth etc...
deltaT = 1 #1 year periods of growth
ageChange = deltaT
dispParam = 8 #dispersalParam
uniqueSp = unique(data$species)

numYears = 5000

#Carrying Capacity
K= 2000

# Can be anywhere from ~0.2 to ~3. Once its above ~1 you get oscillation above and below carrying capacity and program
#gets a fair bit slower as you near 3
growthrates = c(1, 1)

#Pathogen initials

#This is how far away a tree can be affected. But its a torus so setting the value to 100 does not man 
#every tree is infected. best values seem to be between 2.5 and 3 (for infected death rate around 0.03)
#Too low and pathogens go extinct, too high and all trees quickly pathogenized
PathDensParam = 2.7



#Propotion of trees within pathogen-infecting distance that get infected each timestep
infectionRates = c(0.7, 0.7)
resistanceStrength = c(0.5,0.5)
for(i in 1:numSpecies){
  data[data$species==i,]$Resistance = resistanceStrength[i]
}

data[data$Pathogen==1,]$pathStrength = data[data$Pathogen==1,]$Resistance + 0.06
#Has to stay lower than infected death rates. Slows growth to carrying capacity as it gets larger but
# tree dynamics seem to be similar between 0.01 and 0.05, but pathogens end up growing much slower
#and have larger swings as death rates increase
healthyDeathRate = 0.02

#Dynamics get much more variable when value is ~0.04 or more. 0.02-0.04 pretty good for 0.7 infection Rate.
infectedDeathRates = c(0.06,0.06)

# Used for summary statistics
speciesPopulation = NULL
pathPop = NULL
agePops = NULL
pathAgePop = NULL
denseSpecies <- NULL
pathStrengths <- NULL
avgResistance <- NULL
### MAIN LOOP ####


for(t in seq(0, numYears, by = deltaT)){ 
  if(nrow(data[data$Pathogen == 1, ]) == 0 ){
    break
  }
  #disperse seed
  adults <- data[data$age > 3, ]
  seedsPerSpecies <- NULL
  #Get number of seeds per species (logarithmic growth)
  for (i in 1:numSpecies){
    speciesPop = nrow(data[data$species == i,])
    numSeeds <- carryingCapacity(i, speciesPop, nrow(data), K)
    seedsPerSpecies <- c(seedsPerSpecies, numSeeds)
  }
  
  Ind = data[1, ]
  babies <- dispersal(adults, numSpecies, seedsPerSpecies, dispParam, xsize, ysize, data)
  data <- bind_rows(data, babies)
  
  #Background death of adults. Death rate decreases with age
  data <- backgroundDeath(data, infectedDeathRates, healthyDeathRate)
  
  #Count number of each species
  countSpecies <- speciesSize(data$species, numSpecies)
  speciesPopulation <-  bind_rows(speciesPopulation, countSpecies)
  
  #pathogens reproduce
  newPathogens <- pathGrowth(data,PathDensParam, xsize, infectionRates)
  data$Pathogen[which(data$ID %in% newPathogens$ID)] <- 1
  
  #Here we add our pathogen strength values to the new pathogens. We make sure that
  #Pathogen strength is greater than tree resistance for infections
  for(i in 1:nrow(data)){
    if(data[i,]$Pathogen==1 & data[i,]$pathStrength == 0){
      newlyInfected <- data[i,]$ID
      data[i,]$pathStrength <- newPathogens[newPathogens$ID == newlyInfected,]$pathStrength
    
    }
  }
  for(i in 1:nrow(data)){
    if(data[i,]$Resistance > (data[i,]$pathStrength - 0.05)){
      data[i,]$Pathogen = 0
      data[i,]$pathStrength = 0
    }
  }
  
  pathStrengths <- c(pathStrengths, mean(data[data$Pathogen==1,]$pathStrength))
  #Count number of pathogens
  countPaths <- data[data$Pathogen==1,]
  pathAges <- speciesSize(countPaths$species, numSpecies)
  pathPop <- bind_rows(pathPop, pathAges)
  
  
  #Supp  functions/Statistics
  ages <- tibble("gen" = t, "babies" = nrow(subset(data, age < 6)), "young" = nrow(subset(data, age >= 6 & age <=20)), "middle aged" = nrow(subset(data, age > 20 &
                                                                                                                                                     age < 50)), 
                 "old" = nrow(subset(data, age >= 50)))
  agePops <- bind_rows(agePops, ages)

  avgResistance <- c(avgResistance, mean(data$Resistance))
  
  #These fucntions plot the trees every generation. I mainly use these when trying to
  #paramterize something so I can see the effects on several generations. Otherwise comment them out
  
  plot(data$xlocation, data$ylocation, col = data$species, pch = 20, xlim = c(1,100), ylim = c(1,100),
       main = "Tree Population", ylab = "Count", xlab = "Generation")
 
  print(t)
  
  #Grow
  data$age = data$age + deltaT
  
}


##### MAIN FUNCTIONS #######

#Carrying Capacity limits birth rates

carryingCapacity = function(species,speciesN, totalN, K){
  growthrate = growthrates[species]
  s = speciesN * growthrate
  popGrowth <- (s*(K-totalN)/K)
  popGrowth <-  ifelse(popGrowth > 0, popGrowth, 0)
  return(popGrowth)
  
}

#Disperse x number of seeds per parent. 

dispersal <- function(Trees, numSpecies, seedsPerSpecies, dispParam, xsize, ysize, data){
  babies = tibble()
  for(i in 1:numSpecies){
    speciesTree <- Trees[Trees$species == i,]
    if(nrow(speciesTree)> 0){
      #Before running this, we calculate seeds per species using the carrying capacity function
      numBabies <- round(seedsPerSpecies[i])
      #take a random sample of conspecific adult trees to disperse the seeds
      #weight_by = 1/Resistance, 
      parents <- slice_sample(speciesTree, n = numBabies, replace = TRUE)
      if (nrow(parents)>0){
        for (i in 1:nrow(parents)){
          parent = parents[i,]
          baby = parent
          baby$age = as.integer(0)
          baby$ID  = max(max(babies$ID),max(data$ID)) + 1
          baby$species  =  parent$species
          babyDist = (rnorm(1, mean=0, sd=dispParam))
          baby$Pathogen = 0
          baby$Resistance = min(max(rnorm(n = 1, mean = parent$Resistance, sd = 0.05), 0),1)
          baby$pathStrength = 0
          r = babyDist * sqrt(runif(1))
          theta = runif(1)*2*pi
          baby$xlocation = ((r*cos(theta)) + parent$xlocation) %% xsize
          baby$ylocation = ((r*sin(theta)) + parent$ylocation) %% ysize
          baby$parentalDistance = abs(r)
          babies <- bind_rows(babies, baby)
        }
      }
      
    }
  }
  return(babies)
}




#Natural Background death
backgroundDeath <- function(data, infectedDeathRates, healthyDeathRate){
  living = c()
  for (i in 1:nrow(data)){
    if (data[i,]$Pathogen == 1){
      pathEffect <- data[i,]$pathStrength
      chanceDeath <- max(((infectedDeathRates[data[i,]$species])*pathEffect), (healthyDeathRate + 0.01))
      alive <- runif(1) > (chanceDeath)
      living[i] = alive
    }
    if(data[i,]$Pathogen == 0){
      alive <- (runif(1) > healthyDeathRate)
      living[i] = alive
    }
    
  }
  return(data[living,])
}




#### virulence FUNCTIONS ######

pathGrowth <- function(data,PathDensParam, xsize, infectionRates){
  newPathogens <- NULL
  infectedTrees <- data[data$Pathogen == 1,]
  healthyTrees <- data[data$Pathogen == 0,]
  if(nrow(infectedTrees) > 0){
    for (i in 1:nrow(infectedTrees)){
      indTree = infectedTrees[i, ]
      newInfections = subset(healthyTrees, xlocation > (indTree$xlocation - (PathDensParam))%% xsize
                             & xlocation < (indTree$xlocation + (PathDensParam))%% xsize
                             & ylocation > (indTree$ylocation - (PathDensParam))%% xsize
                             & ylocation < (indTree$ylocation + (PathDensParam))%% xsize
                             & species == indTree$species)
      newInfections <- slice_sample(newInfections, n = round(nrow(newInfections)* infectionRates[indTree$species]))
      newInfections$pathStrength = rnorm(n = nrow(newInfections), mean = indTree$pathStrength, sd = 0.05)
      newInfections$pathStrength[newInfections$pathStrength>1] <- 1
      newPathogens <- bind_rows(newPathogens, newInfections)
      
    }
  }
  if (nrow(infectedTrees) == 0){
    newPathogens <- NULL
  }
  
  newPathogens <- newPathogens[!duplicated(newPathogens$ID),]
  return(newPathogens)
}



### SUPPLEMENTAL FUNCTIONS###

  
#Plot trees 
plotSpecies(speciesPopulation)
totals <- speciesPopulation %>% mutate(total = rowSums(across(everything())))
plot(totals$total, pch = 20, col = "cadetblue3")

#Plot Pathogens
plotPaths(pathPop)
totals <- pathPop %>% mutate(total = rowSums(across(everything())))
plot(totals$total, pch = 20, col = "cadetblue3")


#Plot both
grid.arrange(plotSpecies(speciesPopulation), plotPaths(pathPop))

#Plot the age distribution over time of trees
agePops1 <-  pivot_longer(agePops, cols = !gen, names_to = "Group", values_to = "Count")
agePopGraph <- ggplot(agePops1, aes(gen, Count, color = Group))+
  geom_line()+
  scale_x_continuous(n.breaks = 10)
agePopGraph


#Plot pathogen strengths over time
par(mfrow = c(1, 2))

PSplot <- plot(pathStrengths, xlab = "Gen", ylab = "Pathogen Strength")

#Plot resistance over time
TRplot <- plot(avgResistance, xlab = "Gen", ylab = "Resistance")


#Plot amount of each species
plotSpecies <- function(speciesPopulation){
  for(i in 1:numSpecies){
    colnames(speciesPopulation)[i] <- paste0("Species", i)
  }
  speciesPopulation <- add_column(speciesPopulation, gen = 1:nrow(speciesPopulation))
  speciesPopulation <- pivot_longer(speciesPopulation, cols = starts_with("Species"), names_to = "Species",
                                    names_prefix = "Species", values_to = "Count")
  ggplot(data = speciesPopulation, aes(gen, Count, col = Species))+
    geom_point()+
    geom_line()+
    ylim(0, max(speciesPopulation$Count))+
    ylab("Tree Density (Inds)")+
    xlab("Generation (years)")+
    scale_color_manual(values=c("#bf5700", "#333f48"))+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          plot.title = element_text(size=16))+
    ggtitle("Neutral Model")
  
}

#Plot amount of pathogenss
plotPaths <- function(speciesPopulation){
  for(i in 1:numSpecies){
    colnames(speciesPopulation)[i] <- paste0("Species", i)
  }
  speciesPopulation <- add_column(speciesPopulation, gen = 1:nrow(speciesPopulation))
  speciesPopulation <- pivot_longer(speciesPopulation, cols = starts_with("Species"), names_to = "Species",
                                    names_prefix = "Species", values_to = "Count")
  ggplot(data = speciesPopulation, aes(gen, Count, col = Species))+
    geom_point()+
    geom_line()+
    ylim(0, max(speciesPopulation$Count))+
    ylab("Pathogen Density (Inds)")+
    xlab("Generation (years)")+
    scale_color_manual(values=c("#bf5700", "#333f48"))+
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          plot.title = element_text(size=23))
  
}

#Called from main for loop to track species
speciesSize <- function(speciesTypes, numSpecies){
  df = tibble()
  for (i in 1:numSpecies){
    df[1, i] = length(speciesTypes[speciesTypes == i])
  }
  return(df)
}
