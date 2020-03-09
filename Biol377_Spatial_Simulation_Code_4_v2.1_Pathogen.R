# Spatial simulation code for University of Canterbury Biol377.
# Pathogen Dynamics Simulation
# Intended to be used in the context of this document: https://docs.google.com/document/d/1rkp9M2OBmIHHmQu4y52eutd-0Ylll0utAj7guDvg8XI/edit?usp=sharing

# I. Dickie,  J. Tylianakis, D. Stouffer & M. Turnbull, School of Biological Sciences, University of Canterbury NZ
# CC-BY license: https://creativecommons.org/licenses/by/3.0/nz/

#################################
## Model 1: Predator-Prey
#################################

# The model assumes four states in order: empty, host (tree), pathogen infected host, resistant host. These are states 0, 1, 2 and 3

#################################
##  Load library (only need to do this once!)
#################################

## If the library is not already installed, you need to run the next 2 lines
chooseCRANmirror() #choose a location near NZ. Australia works well when it isn't on fire.
install.packages("simecol")

#load the library
library(simecol)

#################################
##  Clean up workspace
#################################

remove(list=ls())

#################################
##  Define transitions function
#################################

transitions <- function(...) 
    # Function used in all versions of the model -- don't modify.
    # Given a set of transition probabilities of length n, return a randomly determined
    # outcome (0 to n-1). Note that probabilities can sum to any value, but it is probably
    # most logical if they sum to 1 or 100.
    {
    probs <- c(...)
    return(sample(0:(length(probs)-1),1, prob=probs))
    }

#################################
## Pathogen spread model
#################################

## Tree recruitment is assumed to be spatially aggregated.  Note that "burnIn" in this model allows the tree population to
## show baseline change before infection.

## Tree recruitment is assumed to be spatially aggregated.  Note that "burnIn" in this model allows the tree population to
## show baseline change before infection.


gridSize <- 100        # Size of the spatial simulation in grid cells
burnIn <- 100          # Number of time steps to run before starting infection AND plotting and saving output
maxTime <- 2000        # Total time steps

mapsFreq <- 5           # How frequently to refresh map. 1 shows every time step 
                        # (slow but pretty), higher values will be faster but jumpier.
                        # set to any value greater than maxTime to suppress mapping.
colVec <- c("white", "#009E73", "#E69F00", "#0072B2")  # set colours 

## define growth, predation and death:
infect <- 1              #probability that a tree will become infected if all neighbours infected
die <- 0.01              #probability that an infected tree will die
resist <- 0.001          #probability that an infected tree will become resistant
recruit <- 0.005         #probability that an empty cell will have a new seedling grow
fitnessCost <- 0.4       #cost of being resistant to plant recruitment. Set to 0 for no cost.
mortality <- 0.002       #probability that a tree will die regardless of disease

initial.trees <- .5      #proportion of cells with tree at first time step
infection.points <- 20   #number of infection start points

# For no particular reason, I've defined the neighbourhood of infection as different from the
# neighbourhood for recruitment, with infection being a shorter-distance dispersal.

# define neighbourhood matrix for infection, in this case being the 4 touching cells:
infectdist <- matrix(c(1, 1, 1,
                  1, 1, 1,
                  1, 1, 1), nrow=3, byrow=TRUE)
# define neighbourhood matrix for recruitment, in this case being the 24 surrounding cells.
recruitdist <- matrix(rep(1, 25), nrow=5)
               
#################################
##  Define transition functions, based on initial state
#################################            

transitions0 <- function(N1, N3, maxNeigh, recruit, fitnessCost) # probability of an empty cell getting a new tree
    {
    if(length(N1)>0) 
        {
        prob1 <- N1*recruit / maxNeigh                          # probability of an empty cell becoming a host tree
        prob3 <- N3*(recruit*(1-fitnessCost)) / maxNeigh        # probability of becoming a resistant tree (note cost)
        prob0 <- 1 - (prob1 + prob3)                            # probability of staying the same
        prob2 <- rep(0, length(N1))                             # probabiliyt of becoming infected (always 0)
        return(mapply(transitions, prob0, prob1, prob2, prob3)) # Use the function transitions with the 3 probability
                                                                # vectors to determine outcomes. Note that the order
                                                                # is important (0, 1, 2)
        } else {
        return(NA)
        }
    }
transitions1 <- function(N2, maxNeigh, infect, mortality)        # probability of susceptible tree being infected
    {
    if(length(N2)>0) 
        {
        prob2 <- (1-mortality) * (N2*infect)/maxNeigh            # Infection. Note 1-mortality is the chance of random death                                                  
        prob0 <- rep(mortality, length(N2))                      # Random death
        prob1 <- 1-(prob2+prob0)                                 # Staying healthy
        prob3 <- rep(0, length(N2))                              # No chance of a healthy tree becoming resistant
        return(mapply(transitions, prob0, prob1, prob2, prob3))
        } else {
        return(NA)
        }
    }
transitions2 <- function(len, death, resist)                     # probability of an infected tree either dying or 
                                                                 # becoming resistant (or lingering)
                                              			         # note that "death" is die+mortality, as two ways to die...
    {
    if(len>0) {         							                   
    	prob1 <- rep(0, len)
    	prob0 <- rep(death, len)
    	prob3 <- rep(resist, len)
    	prob2 <- 1-(prob0+prob3)                                  
        return(mapply(transitions, prob0, prob1, prob2, prob3))
	    } else {
	    return(NA)
	    }
    } 
transitions3 <- function(len, mortality)                         # probability of resistant tree dying
    {
    if(len>0) {							                         # The if statement is for when there are no state 3 cells
    	prob2 <- rep(0,len)                                                            
    	prob0 <- rep(mortality, len)
    	prob1 <- rep(0, len)  
    	prob3 <- 1-prob0
    	return(mapply(transitions, prob0, prob1, prob2, prob3))
        } else {
        return(NA)
        }   
    }
                                    
#################################
##  Initialise map and define empty variables to hold results
#################################
              
map0 <- matrix(sample(c(0,1), gridSize^2, replace=TRUE, prob=c(1-initial.trees, initial.trees)),  nrow=gridSize)
image(map0, col=colVec[1:(max(map0)+1)])
susceptible <- infected <- resistant <- c()                                       # set up "recording" vectors for results
              
#################################
##  Run model
#################################             

for(time in 1:maxTime)                                  
    {
    N1r <- neighbours(map0, state=1, wdist=recruitdist, bounds=1)   #num neighbours of state 1 in recruitment neighbourhood
    N3r <- neighbours(map0, state=3, wdist=recruitdist, bounds=1)   #num neighbours of state 3 in recruitment neighbourhood
    N2i <- neighbours(map0, state=2, wdist=infectdist, bounds=1)    #num neighbours of state 2 in infection neighbourhood
    map1 <- map0                                        # 
    map1[map0 == 0] <- transitions0(N1r[map0==0],N3r[map0==0], maxNeigh = sum(recruitdist)-1, recruit, fitnessCost)
    map1[map0 == 1] <- transitions1(N2i[map0 == 1], maxNeigh=sum(infectdist)-1, infect, mortality)
    map1[map0 == 2] <- transitions2(length(map1[map0 == 2]), mortality+die, resist)
    map1[map0 == 3] <- transitions3(length(map1[map0 == 3]), mortality)
    #After burn-in period, produce visual images (if maps == TRUE) and save output
    if(time == burnIn)
        {
        map1[sample(1:gridSize^2, infection.points)] <- 2
        }
    if(time >= burnIn)
        {
        if(time %% mapsFreq == 0)
            {
            image(map1, add=TRUE,         			     # Draw a map, adding time step to label    
                col=colVec[sort(unique(c(map1)))+1])                # This colors the map correctly
            } #end if(time %% mapsFreq == 0)
        susceptible  <- append(susceptible , sum(map1 == 1))
        infected <- append(infected, sum(map1 == 2))
        resistant <- append(resistant, sum(map1 == 3))
        }  #end if(time > burnIn)
    map0 <- map1 
    }  #end for(time)

#################################
##  Plot model output
#################################
plot(susceptible, type="l", col=colVec[2], ylim=c(0, max(susceptible, infected, resistant)), xlab="Time", ylab="Populations")
points(infected, type="l", col=colVec[3])
points(resistant, type="l", col=colVec[4])
# Add a legend.
legend(1, 0.8 *max(susceptible, infected, resistant), legend=c("Susceptible", "Infected", "Resistant"), col=c(colVec[2], colVec[3]), lty=1:1, cex=0.8)

