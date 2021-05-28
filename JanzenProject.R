#Load packages
library(tidyverse) 
library(gridExtra) 


### Initialize params ###

# Set up a grid to represent land
#Place individuals on plot
xsize = ysize = 100

numInd = 100
numSpecies = 2
data <-  tibble("ID" = 1:numInd, "age" = 3, 
                species = sample(numSpecies, numInd, replace = TRUE), "xlocation" = runif(numInd, 0, 100), 
                "ylocation" = runif(numInd,0,100), "parentalDistance" = 0, "Pathogen" = 0)

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
numYears = 1000



#Carrying Capacity
K= 2000

# Can be anywhere from ~0.2 to ~3. Once its above ~1 you get oscillation above and below carrying capacity and program
#gets a fair bit slower as you near 3
growthRate1 = 1
growthRate2= 1

#Pathogen initials

#This is how far away a tree can be affected. But its a torus so setting the value to 100 does not man 
#every tree is infected. best values seem to be between 2.5 and 3 (for infected death rate around 0.03)
#Too low and pathogens go extinct, too high and all trees quickly pathogenized
PathDensParam = 2.7

infectionRate1 = 0.7
infectionRate2 = 0.7

#Has to stay lower than infected death rates. Slows growth to carrying capacity as it gets larger but
# tree dynamics seem to be similar between 0.01 and 0.05, but pathogens end up growing much slower
#and have larger swings as death rates increase
healthyDeathRate = 0.02

#Dynamics get much more variable when value is ~0.04 or more. 0.02-0.04 pretty good
infectedDeathRate1 = 0.05
infectedDeathRate2 = 0.05
# Used for summary statistics
speciesPopulation = NULL
pathPop = NULL
agePops = NULL
pathAgePop = NULL
denseSpecies <- NULL

### MAIN LOOP ####


for(t in seq(0, numYears, by = deltaT)){ 
  
  #disperse seed
  adults <- data[data$age > 3, ]
  seedsPerSpecies <- NULL
  for (i in 1:numSpecies){
    speciesPop = nrow(data[data$species == i,])
    numSeeds <- carryingCapacity(i, speciesPop, nrow(data), K)
    seedsPerSpecies <- c(seedsPerSpecies, numSeeds)
  }
  Ind = data[1, ]
  babies <- dispersal(adults)
  data <- bind_rows(data, babies)
  
  #Background death of adults. Death rate decreases with age
  data <- backgroundDeath(data)
  
  #Count number of each species
  countSpecies <- speciesSize(data$species)
  speciesPopulation <-  bind_rows(speciesPopulation, countSpecies)
  
  #pathogens reproduce
  newPathogens <- pathGrowth(data)
  data$Pathogen[which(data$ID %in% newPathogens$ID)] <- 1
  
  #Supplemental
  countPaths <- data[data$Pathogen==1,]
  pathAges <- speciesSize(countPaths$species)
  pathPop <- bind_rows(pathPop, pathAges)
  
  
  #Supp  functions/Statistics
  ages <- tibble("gen" = t, "babies" = nrow(subset(data, age < 6)), "young" = nrow(subset(data, age >= 6 & age <=20)), "middle aged" = nrow(subset(data, age > 20 &
                                                                                                                                                     age < 50)), 
                 "old" = nrow(subset(data, age >= 50)))
  agePops <- bind_rows(agePops, ages)

  
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
  if(species==1){
    s = speciesN * growthRate1
  }else{
    s = speciesN * growthRate2
  }
  popGrowth <- (s*(K-totalN)/K)
  popGrowth <-  ifelse(popGrowth > 0, popGrowth, 0)
  return(popGrowth)
  
}

#Disperse x number of seeds per parent. 

dispersal <- function(Trees){
  babies = NULL
  for(i in 1:numSpecies){
    speciesTree <- Trees[Trees$species == i,]
    if(nrow(speciesTree)> 0){
      #Before running this, we calculate seeds per species using the carrying capacity function
      numBabies <- ceiling(seedsPerSpecies[i])
      #take a random sample of conspecific adult trees to disperse the seeds
      parents <- slice_sample(speciesTree, n = numBabies, replace = TRUE)
      if (nrow(parents)>0){
        for (i in 1:nrow(parents)){
          parent = parents[i,]
          baby = parent[1, ]
          baby$age = as.integer(0)
          if(parent$species == 1){
            baby$ID  = max(data$ID) + i 
          }
          else{
            baby$ID  = max(babies$ID) + i 
          }
          baby$species  =  parent$species
          babyDist = (rnorm(1, mean=0, sd=dispParam))
          baby$Pathogen = 0
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
backgroundDeath <- function(population){
  living = c()
  probDeath = c()
  for (i in 1:nrow(population)){
    if (data[i,]$Pathogen == 1 & data[i,]$species ==1){
      alive <- (runif(1) > infectedDeathRate1)
      living[i] = alive
    }else if (data[i,]$Pathogen == 1 & data[i,]$species ==2) {
      alive <- (runif(1) > infectedDeathRate2)
      living[i] = alive
    }
    else{
      alive <- (runif(1) > healthyDeathRate)
      living[i] = alive
    }
    
  }
  return(population[living,])
}




#### virulence FUNCTIONS ######

pathGrowth <- function(trees){
  newPathogens <- NULL
  infectedTrees <- data[data$Pathogen == 1,]
  healthyTrees <- data[data$Pathogen == 0,]
  for (i in 1:nrow(infectedTrees)){
    indTree = infectedTrees[i, ]
    newInfections = subset(healthyTrees, xlocation > (indTree$xlocation - (PathDensParam))%% xsize
                           & xlocation < (indTree$xlocation + (PathDensParam))%% xsize
                           & ylocation > (indTree$ylocation - (PathDensParam))%% xsize
                           & ylocation < (indTree$ylocation + (PathDensParam))%% xsize
                           & species == indTree$species)
    if (indTree$species==1){
      newInfections <- slice_sample(newInfections, n = round(nrow(newInfections)* infectionRate1))
      newPathogens <- bind_rows(newPathogens, newInfections)
    }else{
      newInfections <- slice_sample(newInfections, n = round(nrow(newInfections)* infectionRate2))
      newPathogens <- bind_rows(newPathogens, newInfections)
    }
    
    
  }
  newPathogens <- newPathogens[!duplicated(newPathogens),]
  return(newPathogens)
}



### SUPPLEMENTAL FUNCTIONS###


#Plot trees 
plotSpecies(speciesPopulation) 
totals <- speciesPopulation %>% mutate(total = (`Species 1` + `Species 2` + `Species 3`))
plot(totals$total, pch = 20, col = "cadetblue3")


#Plot Pathogens
plotSpecies(pathPop)
totals <- pathPop %>% mutate(total = (`Species 1` + `Species 2` + `Species 3`))
plot(totals$total, pch = 20, col = "cadetblue3")

#Plot the age distribution over time of trees
agePops1 <-  pivot_longer(agePops, cols = !gen, names_to = "Group", values_to = "Count")
agePopGraph <- ggplot(agePops1, aes(gen, Count, color = Group))+
  geom_line()
agePopGraph



#Plot amount of each species
plotSpecies <- function(speciesPopulation){
  speciesPopulation <- add_column(speciesPopulation, gen = 1:nrow(speciesPopulation))
  speciesPopulation <- pivot_longer(speciesPopulation, cols = starts_with("Species"), names_to = "Species",
                                    names_prefix = "Species ", values_to = "Count")
  ggplot(data = speciesPopulation, aes(gen, Count, col = Species))+
    geom_point()+
    geom_line()+
    ylim(0, max(speciesPopulation$Count))+
    ylab("Tree Count")
  
}

#Called from main for loop to track species
speciesSize <- function(speciesTypes){
  a = speciesTypes[speciesTypes == 1]
  b = speciesTypes[speciesTypes == 2]
  c = speciesTypes[speciesTypes == 3]
  
  df <- tibble("Species 1" = length(a), "Species 2" = length(b), "Species 3" = length(c))
  return(df)
}




