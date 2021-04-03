### This version of the model still has a "phenomenological" density function that won't be in the final code.
#BetaCode branch has a version without this function, and with slightly more organized code, but not totally parameterized yet



#Load packages
library(tidyverse) 
library(gridExtra)


### Initialize params ###

# Set up a grid to represent land
#Place individuals on plot
xsize = ysize = 100
numInd = 50
numSpecies = 3
data <-  tibble("ID" = 1:numInd, "age" = sample(150, numInd, replace = TRUE), 
                 species = sample(numSpecies, numInd, replace = TRUE), "xlocation" = runif(numInd, 0, 100), 
                 "ylocation" = runif(numInd,0,100), "parentalDistance" = 0)

herbivores <- data %>%  slice_sample(n = 25) %>% select(-c(ID, age, parentalDistance))
herbivores <- herbivores %>% mutate(ID = 1:nrow(herbivores), age = 0)
plot(data$xlocation, data$ylocation, col = data$species, pch = 20, xlim = c(1,100), ylim = c(1,100))

#Set up rates of birth/death/growth etc...
deltaT = 1 #1 year periods of growth
ageChange = deltaT
maxAge = 150
seeds = 2
dispParam = 10 #dispersalParam
DensParam = 4 #DensityParam
DensParam1 = 2
uniqueSp = unique(data$species)
numYears = 10000

#Herbivore initials
HerbDensParam = 4.8
HerbEffectiveness = 2.3

speciesPop = NULL
herbPop = NULL
agePops = NULL
herbAgePop = NULL





### MAIN LOOP ####


for(t in seq(0, numYears, by = deltaT)){ 
  #grow
  data$age = data$age + deltaT
  herbivores$age = herbivores$age + 1
  
  #disperse seeds/
  Ind = data[1, ]
  babies <- dispersal(data)
  #survival of seeds
  babies <- seedPredation(babies, babies$parentalDistance)
  data <- bind_rows(data, babies)
  #Background death of adults. Death rate decreases with age
  data <- backgroundDeath(data, data$age)
  
  #Conspecific density death
  livePastDens = NULL
  for (k in 1:length(uniqueSp)){
    sp = uniqueSp[k]
    subPop = subset(data, species == sp)
    survived <- densityPredation(subPop)
    livePastDens <-  bind_rows(livePastDens, survived)
    livePastDens <- livePastDens[order(livePastDens$ID), ]
  }
  data <- livePastDens
  #trees over max age die
  
  countSpecies <- speciesSize(data$species)
  speciesPop <-  bind_rows(speciesPop, countSpecies)
  
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
  
  
  #Supp  functions
  ages <- tibble("gen" = t, "babies" = nrow(subset(data, age < 6)), "young" = nrow(subset(data, age >= 6 & age <=20)), "middle aged" = nrow(subset(data, age > 20 &
                                                                                                   age < 50)), 
                 "old" = nrow(subset(data, age >= 50)))
  
  herbAges <- tibble("gen" = t, "babies" = nrow(subset(herbivores, age < 6)), 
                     "young" = nrow(subset(herbivores, age >= 6 & age <=20)),
                     "middle aged" = nrow(subset(herbivores, age > 20 &                                                                                                                                 age < 50)), 
                     "old" = nrow(subset(herbivores, age >= 50)))
  
  agePops <- bind_rows(agePops, ages)
  herbAgePop <- bind_rows(herbAgePop, herbAges)
  #plot(data$xlocation, data$ylocation, col = data$species, pch = 20, xlim = c(1,100), ylim = c(1,100))
  #plot(herbivores$xlocation, herbivores$ylocation, col = herbivores$species, pch = 20, xlim = c(1,100), ylim = c(1,100))
  
  print(t)

}



##### MAIN FUNCTIONS #######

#Disperse 10 seeds per parent
dispersal <- function(Trees){
    babies = NULL
    for (i in 1:nrow(Trees)){
      eachAdult = Trees[i,]
      for (j in 1:seeds){
        baby = Trees[1, ]
        baby$age = as.integer(0)
        baby$ID  = max(data$ID) + (10*(i-1)) + j
        baby$species  =  eachAdult$species
        babyDist = (rnorm(1, mean=0, sd=dispParam))
        xDirection = runif(1)
        baby$xlocation = ((xDirection*babyDist) + eachAdult$xlocation) %% xsize
        baby$ylocation = (((1-xDirection) * babyDist) + eachAdult$ylocation) %% ysize
        baby$parentalDistance = abs(babyDist)
        babies <- bind_rows(babies, baby)
      }
    }
    return(babies)
    
}



#Seeds predated based on distance from parent
seedPredation <- function(Babies, parentalDist){
  babySurv = NULL
  for (i in 1:nrow(Babies)){
    dist = parentalDist[i]
    danger = runif(1)*(1/(dist))
    indLive = danger < .07
    babySurv[i] = indLive
  }
  return(Babies[babySurv, ])
}

#predation based on conspecifics nearby
densityPredation <- function(subPop){
  survived = NULL
  #dead = NULL
  for (i in 1:nrow(subPop)){
    #subPop <- subset(subPop, !(ID %in% dead$ID) )
    neighbors <-  NULL
    ind <- subPop[i, ]
    xNeighbors = data.frame(subset(subPop, xlocation <= ind$xlocation + DensParam & 
                                     xlocation >= ind$xlocation - DensParam
                                   & xlocation != ind$xlocation))
    yNeighbors = data.frame(subset(subPop, ylocation <= ind$ylocation + DensParam & 
                                     ylocation >= ind$ylocation - DensParam
                                   & ylocation != ind$ylocation))
    neighbors = intersect(xNeighbors, yNeighbors)
    neighbors = unique(neighbors)
    prob <- ifelse(nrow(neighbors)>0,  (DensParam1 * (ind$age+1))/(((max(neighbors$age))+1)*
                                                                     (nrow(neighbors))),
                   1.01)
    luck = runif(1)
    if(runif(1) < prob){
      survived <- bind_rows(ind , survived)
    } 
  }
  return(survived)
}

#Natural Background death
backgroundDeath <- function(population, age){
  living = c()
  probDeath = c()
  for (i in 1:nrow(population)){
    probDeath[i] <- ifelse(age[i] > 0,  0.3/(age[i]), 0.3)
    alive <- (runif(1) > probDeath[i])
    living[i] = alive
  }
  return(population[living,])
}






#### HERBIVORY STUFF ######


#Herbivory function: herbivores kill trees
herbivory <- function(trees, herbivores){
  deadTrees = NULL
  for (i in 1:nrow(herbivores)){
    indHerbivore = herbivores[i, ]
    infected = data[data$xlocation == indHerbivore$xlocation  & data$ylocation == indHerbivore$ylocation, ]
    probSurv = (((infected$age+3)) / (((indHerbivore$age)+1)*HerbEffectiveness))
    if(probSurv < 1){
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
    # for (j in 1:7){
    #   xDisplacement <- (indHerbivore$xlocation + rnorm(1,0, sd = 6))
    #   yDisplacement <- (indHerbivore$ylocation + rnorm(1,0, sd = 6))
    #   infectedTree <- subset(trees, xlocation > xDisplacement - 4 & xlocation < xDisplacement + 4 &
    #                            ylocation > yDisplacement - 4 & ylocation < yDisplacement + 4 &
    #                            xlocation != indHerbivore$xlocation & ylocation != indHerbivore$ylocation &
    #                            species == indHerbivore$species)
    #   infectedTrees <- bind_rows(infectedTrees, infectedTree)
    # }
    infectedTree <- subset(trees, xlocation > indHerbivore$xlocation - HerbDensParam & xlocation < indHerbivore$xlocation + 
                             HerbDensParam & ylocation > indHerbivore$ylocation - 
                             HerbDensParam & ylocation < indHerbivore$ylocation + HerbDensParam &
                            xlocation != indHerbivore$xlocation & ylocation != indHerbivore$ylocation &
                           species == indHerbivore$species)
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
grid.arrange(plotSpecies(speciesPop), plotHerbSpecies(herbPop), ncol = 1)

#Plot the age distribution over time of trees
agePops1 <-  pivot_longer(agePops, cols = !gen, names_to = "Group", values_to = "Count")
agePopGraph <- ggplot(agePops1, aes(gen, Count, color = Group))+
  geom_line()

herbAgePop1 <-  pivot_longer(herbAgePop, cols = !gen, names_to = "Group", values_to = "Count")
herbAgeGraph <- ggplot(herbAgePop1, aes(gen, Count, color = Group))+
  geom_line()

grid.arrange(agePopGraph, herbAgeGraph, ncol = 1)


#plot adult pop
herbivoreOld <- subset(herbivores, age > 50)
herbAgePopOld <-  pivot_longer(herbivoreOld, cols = !gen, names_to = "Group", values_to = "Count")
herbAgeGraph <- ggplot(herbAgePopOld, aes(gen, Count, color = Group))+
  geom_line()


agePopGraph <- ggplot(agePops, aes(gen, old))+
  geom_line()
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
plotSpecies <- function(speciesPop){
  speciesPop <- add_column(speciesPop, gen = 1:nrow(speciesPop))
  speciesPop <- pivot_longer(speciesPop, cols = starts_with("Species"), names_to = "Species",
                             names_prefix = "Species ", values_to = "Count")
  ggplot(data = speciesPop, aes(gen, Count, col = Species))+
    geom_point()+
    geom_line()+
    ylim(0, max(speciesPop$Count))+
    ylab("Tree Count")
  
}
plotHerbSpecies <- function(speciesPop){
  speciesPop <- add_column(speciesPop, gen = 1:nrow(speciesPop))
  speciesPop <- pivot_longer(speciesPop, cols = starts_with("Species"), names_to = "Species",
                             names_prefix = "Species ", values_to = "Count")
  ggplot(data = speciesPop, aes(gen, Count, col = Species))+
    geom_point()+
    geom_line()+
    ylim(0, max(speciesPop$Count))+
    ylab("Herbivore Count")
  
}
#If want to track pop/dens put this in main loop and recreate variables popSize and avgDens
#Then use plotSize and densSize functions
popSize = c(popSize, nrow(data))
currentDens <- densSize(data$xlocation, data$ylocation)
avgDens <- c(avgDens, currentDens)

