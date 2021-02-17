#Load packages
library(tidyverse) 
# Set up a grid to represent land
#Place individuals on plot
xsize = ysize = 100
numInd = 50
numSpecies = 3
data <-  tibble("ID" = 1:numInd, "age" = sample(150, 50, replace = TRUE), 
                 species = sample(numSpecies, 50, replace = TRUE), "xlocation" = runif(50, 0, 100), 
                 "ylocation" = runif(50,0,100), "parentalDistance" = 0)
plot(data$xlocation, data$ylocation, col = data$species, pch = 20, xlim = c(1,100), ylim = c(1,100))

#Set up rates of birth/death/growth etc...
deltaT = 1 #1 year periods of growth
ageChange = deltaT
maxAge = 150
seeds = 10
dispParam = 10 #dispersalParam
DensParam = 2.5 #DensityParam
DensParam1 = 200
uniqueSp = unique(data$species)
numYears = 50
popSize = NULL
avgDens = NULL
### MAIN LOOP ####


for(t in seq(0, numYears, by = deltaT)){ 
  #grow
  data$age = data$age + deltaT
  #disperse seeds/
  Ind = data[1, ]
  babies <- dispersal(Ind, data$species, data$xlocation, data$ylocation)
  #survival of seeds
  babies <- seedPredation(babies, babies$parentalDistance)
  data <- bind_rows(data, babies)
  #Background death of adults. Death rate decreases with age
  data <- backgroundDeath(data, data$age)
  #Conspecific density death
  livePastDens = NULL
  for (i in 1:length(uniqueSp)){
    sp = uniqueSp[i]
    subPop = subset(data, species == sp)
    survived <- densityPredation(subPop, subPop$age, subPop$ID, subPop$xlocation, subPop$ylocation)
    livePastDens <-  bind_rows(livePastDens, survived)
    livePastDens <- livePastDens[order(livePastDens$ID), ]
  }
  data <- livePastDens
  #trees over max age die
  data  <- filter(data, age <= 150)

  #Supplemental funs

  popSize = c(popSize, nrow(data))

  currentDens <- densSize(data$xlocation, data$ylocation)
  avgDens <- c(avgDens, currentDens)


  plot(data$xlocation, data$ylocation, col = data$species, pch = 20, xlim = c(1,100), ylim = c(1,100))

}

}






##### MAIN FUNCTIONS #######

#Disperse 10 seeds per parent
dispersal <- function(Ind, species, xlocation, ylocation){
    babies = NULL
    Adults = tibble("species" = species, "x" = xlocation, "y" = ylocation)
    for (i in 1:nrow(Adults)){
      eachAdult = Adults[i,]
      for (j in 1:seeds){
        baby = Ind
        baby$age = as.integer(0)
        baby$ID  = max(data$ID) + (10*(i-1)) + j
        baby$species  =  eachAdult$species
        babyDist = abs(rnorm(1, mean=0, sd=dispParam))
        xDirection = runif(1)
        baby$xlocation = ((xDirection*babyDist) + eachAdult$x) %% xsize
        baby$ylocation = (((1-xDirection) * babyDist) + eachAdult$y) %% ysize
        baby$parentalDistance = babyDist
        babies <- bind_rows(babies, baby)
      }
    }
    return(babies)
    #babies going over borders of ecosystem, fix  that later
}



#Seeds predated based on distance from parent
seedPredation <- function(Babies, parentalDist){
  babySurv = NULL
  for (i in 1:nrow(Babies)){
    dist = parentalDist[i]
    danger = 1/dist
    luck = runif(1)*danger
    indDeath = luck > 0.1
    babySurv[i] = indDeath
  }
  return(Babies[!babySurv, ])
}

#predation based on conspecifics nearby
densityPredation <- function(subPop, age, ID, xlocation, ylocation){
  survived = NULL
  dead = NULL
  indSpecies <-  data.frame(cbind("ID" = ID, "age" = age, "xlocation" = xlocation, "ylocation" = ylocation))
  for (i in 1:nrow(subPop)){
    neighbors <-  NULL
    ageInd = age[i]
    indX <-  xlocation[i]
    indY <-  ylocation[i]
    xNeighbors = data.frame(subset(indSpecies, xlocation <= indX + DensParam & xlocation >= indX 
                                   - DensParam
                                   & xlocation != indX))
    yNeighbors = data.frame(subset(indSpecies, ylocation <= indY + DensParam & 
                                     ylocation >= indY - DensParam
                                   & ylocation != indY))
    neighbors = bind_rows(neighbors,xNeighbors)
    neighbors = bind_rows(neighbors, yNeighbors)
    neighbors = unique(neighbors)
    prob <- ifelse(nrow(neighbors)>0,  (DensParam1/nrow(neighbors) * (ageInd+1))/((max(neighbors$age))+1),
           1.01)
    ifelse(runif(1) < prob, survived <-  bind_rows(subPop[i, ], survived), next)
  }
  return(survived)
}

#Natural Background death
backgroundDeath <- function(population, age){
  living = c()
  probDeath = c()
  for (i in 1:nrow(population)){
    probDeath[i] <- ifelse(age[i] > 0,  0.4/age[i], 0.7)
    alive <- (runif(1) > probDeath[i])
    living[i] = alive
  }
  return(population[living,])
}


### SUPPLEMENTAL FUNCTIONS###

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

plot(densSize)

#Only one tree in comp should die if multiple trees in area?
#Redo function driving density dependence



