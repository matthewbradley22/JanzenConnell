### Updated with the new carrying capacity equation. Still parameterizing a bit
#Load packages

#Load packages

#Load packages

library(tidyverse) 
library(gridExtra)


### Initialize params ###

# Set up a grid to represent land
#Place individuals on plot
xsize = ysize = 100

numInd = 100
numSpecies = 3
data <-  tibble("ID" = 1:numInd, "age" = 3, 
                 species = sample(numSpecies, numInd, replace = TRUE), "xlocation" = runif(numInd, 0, 100), 
                 "ylocation" = runif(numInd,0,100), "parentalDistance" = 0)

#Ignore:
#sample(150, numInd, replace = TRUE)

herbivores <- data %>%  slice_sample(n = 25) %>% select(-c(ID, age, parentalDistance))
herbivores <- herbivores %>% mutate(ID = 1:nrow(herbivores), age = 0)
plot(data$xlocation, data$ylocation, col = data$species, pch = 20, xlim = c(1,100), ylim = c(1,100),
     xlab = "X-Coord", ylab = "Y-Coord", main = "Initial Trees")



#Set up rates of birth/death/growth etc...
deltaT = 1 #1 year periods of growth
ageChange = deltaT
dispParam = 8 #dispersalParam
DensParam = 0.3 #DensityParamx
# DensParam1 = 2
uniqueSp = unique(data$species)
numYears = 10000



#Carrying Capacity
K= 3500
growthRate = 0.3

#Herbivore initials
HerbDensParam = 4
HerbEffectiveness = 0.84

# Used for summary statistics
speciesPopulation = NULL
herbPop = NULL
agePops = NULL
herbAgePop = NULL
denseSpecies <- NULL

### MAIN LOOP ####


for(t in seq(0, numYears, by = deltaT)){ 
  #grow
  data$age = data$age + deltaT
  herbivores$age = herbivores$age + 1
  
  #disperse seed
  adults <- data[data$age > 3, ]
  seedsPerSpecies <- NULL
  for (i in 1:numSpecies){
    speciesPop = nrow(data[data$species == i,])
    numSeeds <- carryingCapacity(speciesPop, nrow(data), K)
    seedsPerSpecies <- c(seedsPerSpecies, numSeeds)
  }
  Ind = data[1, ]
  babies <- dispersal(adults)
  #survival of seeds
  #babies <- seedPredation(babies, babies$parentalDistance)
  data <- bind_rows(data, babies)
  
  #Background death of adults. Death rate decreases with age
  data <- backgroundDeath(data)
  
  #trees over max age die
  
  countSpecies <- speciesSize(data$species)
  speciesPopulation <-  bind_rows(speciesPopulation, countSpecies)
  
  #Herbivores die
  herbivores <- herbivoreDeath(data, herbivores)
  
  
  #Herbivores reproduce
  newHerbivores <- herbivoreGrowth(data, herbivores)
  herbivores <- bind_rows(herbivores, newHerbivores)
  herbivores <- herbivores[!duplicated(herbivores$xlocation) & !duplicated(herbivores$ylocation),]
  
  #Supplemental
  countHerbs <- speciesSize(herbivores$species)
  herbPop <- bind_rows(herbPop, countHerbs)
  
  #Herbivores damage trees
  data <-  herbivory(data, herbivores)
  
  
  #Supp  functions/Statistics
  ages <- tibble("gen" = t, "babies" = nrow(subset(data, age < 6)), "young" = nrow(subset(data, age >= 6 & age <=20)), "middle aged" = nrow(subset(data, age > 20 &
                                                                                                   age < 50)), 
                 "old" = nrow(subset(data, age >= 50)))
  
  herbAges <- tibble("gen" = t, "babies" = nrow(subset(herbivores, age < 6)), 
                     "young" = nrow(subset(herbivores, age >= 6 & age <=20)),
                     "middle aged" = nrow(subset(herbivores, age > 20 &                                                                                                                                 age < 50)), 
                     "old" = nrow(subset(herbivores, age >= 50)))
  currentDens <- tibble("gen" = t, "num" = speciesDense(data))
  
  denseSpecies <- bind_rows(denseSpecies, currentDens)
  agePops <- bind_rows(agePops, ages)
  herbAgePop <- bind_rows(herbAgePop, herbAges)
  
  #These fucntions plot the trees every generation. I mainly use these when trying to
  #paramterize something so I can see the effects on several generations. Otherwise comment them out
  
  plot(data$xlocation, data$ylocation, col = data$species, pch = 20, xlim = c(1,100), ylim = c(1,100),
       main = "Tree Population Size", ylab = "Count", xlab = "Generation")
  #plot(herbivores$xlocation, herbivores$ylocation, col = herbivores$species, pch = 20, xlim = c(1,100), ylim = c(1,100))
  
  print(t)

}


##### MAIN FUNCTIONS #######

#Carrying Capacity limitis birth rates

carryingCapacity = function(speciesN, totalN, K){
  s = speciesN * growthRate
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
          baby$ID  = max(data$ID) + i 
          baby$species  =  parent$species
          babyDist = (rnorm(1, mean=0, sd=dispParam))
          r = babyDist * sqrt(runif(1))
          theta = runif(1)*2*pi
          baby$xlocation = ((r*cos(theta)) + parent$xlocation) %% xsize
          baby$ylocation = (((r*sin(theta)) * babyDist) + parent$ylocation) %% ysize
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
    alive <- (runif(1) > 0.01)
    living[i] = alive
  }
  return(population[living,])
}




#### HERBIVORY FUNCTIONS ######


#Herbivory function: herbivores kill trees
herbivory <- function(trees, herbivores){
  deadTrees = NULL
  for (i in 1:nrow(herbivores)){
    indHerbivore = herbivores[i, ]
    infected = data[data$xlocation == indHerbivore$xlocation  & data$ylocation == indHerbivore$ylocation, ]
    probSurv = (((infected$age+1)) / (((indHerbivore$age)+1)*HerbEffectiveness))
    if(probSurv < abs(rnorm(1,mean = 1, sd = 0.2))){
      deadTrees <- bind_rows(deadTrees, infected)
    }
  }
  deadID = deadTrees$ID
  return(subset(trees, !(ID %in% deadID)))
}

herbivoreGrowth <- function(trees, herbivores){
  infectedTrees <- NULL
  for (i in 1:nrow(herbivores)){
    indHerbivore <- herbivores[i, ]
    infectedTree <- subset(trees, xlocation > indHerbivore$xlocation - (HerbDensParam)
                           & xlocation < indHerbivore$xlocation + (HerbDensParam )
                           & ylocation > indHerbivore$ylocation - (HerbDensParam)
                           & ylocation < indHerbivore$ylocation + (HerbDensParam )
                           & species == indHerbivore$species)
    infectedTree <- subset(infectedTree, (xlocation != indHerbivore$xlocation | ylocation != indHerbivore$ylocation))
    infectedTrees <- bind_rows(infectedTrees, infectedTree)
  }
  if(nrow(infectedTrees) > 0){
    newHerbivores <- infectedTrees %>%  select(-c(parentalDistance, ID, age)) %>% 
      mutate(ID = seq(max(herbivores$ID)+1, 
      max(herbivores$ID)+nrow(infectedTrees)), age = 0)
  }else{
    newHerbivores <- NULL
  }
  return(newHerbivores)
}

herbivoreDeath <- function(trees, herbivores){
  dead = NULL
  totalInf = NULL
  for (i in 1:nrow(herbivores)){
    indHerbivore = herbivores[i, ]
    infected = data[data$xlocation == indHerbivore$xlocation  & data$ylocation == indHerbivore$ylocation, ]
    if (nrow(infected) < 1){
      dead <- bind_rows(dead, indHerbivore)
    }
  }
  herbID = dead$ID
  return(subset(herbivores, !(ID %in% herbID)))
}


### SUPPLEMENTAL FUNCTIONS###


#Plot pathogen vs trees
grid.arrange(plotSpecies(speciesPopulation), plotHerbSpecies(herbPop), ncol = 1)

totals <- speciesPopulation %>% mutate(total = (`Species 1` + `Species 2` + `Species 3`))
plot(totals$total, pch = 20, col = "cadetblue3")
#Plot the age distribution over time of trees
agePops1 <-  pivot_longer(agePops, cols = !gen, names_to = "Group", values_to = "Count")
agePopGraph <- ggplot(agePops1, aes(gen, Count, color = Group))+
  geom_line()

herbAgePop1 <-  pivot_longer(herbAgePop, cols = !gen, names_to = "Group", values_to = "Count")
herbAgeGraph <- ggplot(herbAgePop1, aes(gen, Count, color = Group))+
  geom_line()

grid.arrange(agePopGraph, herbAgeGraph, ncol = 1)


#plot adult pop
treesOld <- subset(agePops1, Group == "old")

agePopGraph <- ggplot(treesOld, aes(gen, Count))+
  geom_line()
agePopGraph
#Plot species density
ggplot(data = denseSpecies, aes(x = gen, y= num))+
  geom_line()+
  labs(y = "Population per 5x5",x = "generation", title = "Conspecific density")

#Used to plot population over time
plotSize <- function(popSize){
  popTime <- tibble("gen" = c(1:length(popSize)), "population" = popSize)
  plot(popTime$gen, popTime$population, type = "b", col = "blue")
}

#used to plot average density for a specific parameter
densSize <- function(xlocation, ylocation){
  Densities = NULL
  for (i in seq(0,95, by = 5)){
    popSegment = subset(xlocation, xlocation > i & xlocation < i+5)
    Densities = c(Densities, length(popSegment))
  }
  return(mean(Densities))
}

#Track population of each species
speciesSize <- function(speciesTypes){
  a = speciesTypes[speciesTypes == 1]
  b = speciesTypes[speciesTypes == 2]
  c = speciesTypes[speciesTypes == 3]
  
  df <- tibble("Species 1" = length(a), "Species 2" = length(b), "Species 3" = length(c))
  return(df)
}


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
plotHerbSpecies <- function(speciesPopulation){
  speciesPopulation <- add_column(speciesPopulation, gen = 1:nrow(speciesPopulation))
  speciesPopulation <- pivot_longer(speciesPopulation, cols = starts_with("Species"), names_to = "Species",
                             names_prefix = "Species ", values_to = "Count")
  ggplot(data = speciesPopulation, aes(gen, Count, col = Species))+
    geom_point()+
    geom_line()+
    ylim(0, max(speciesPopulation$Count))+
    ylab("Herbivore Count")
  
}

speciesDense <- function(population){
  totalNeighbors <- tibble("set" = NULL, "neighbors" = NULL)
  for (i in seq(0,100,5)){
      numNeighbors <- nrow(subset(population, species == 1 & xlocation > (i-5) & xlocation <= i
                                  & ylocation > (i-5) & ylocation < (i)))
      currentNeighbors <- tibble("set" = (i/5), "neighbors" = numNeighbors)
      totalNeighbors <- bind_rows(totalNeighbors, currentNeighbors)
  }
  avg = mean(totalNeighbors$neighbors)
  return(avg)
  }
#If want to track pop/dens put this in main loop and recreate variables popSize and avgDens
#Then use plotSize and densSize functions
popSize = c(popSize, nrow(data))
currentDens <- densSize(data$xlocation, data$ylocation)
avgDens <- c(avgDens, currentDens)





