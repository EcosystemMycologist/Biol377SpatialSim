# Spatial simulation code for University of Canterbury Biol377
# Predator Prey Dynamics Simulation
# Intended to be used in the context of this document: https://docs.google.com/document/d/1rkp9M2OBmIHHmQu4y52eutd-0Ylll0utAj7guDvg8XI/edit?usp=sharing

# I. Dickie,  J. Tylianakis, D. Stouffer & M. Turnbull, School of Biological Sciences, University of Canterbury NZ
# CC-BY license: https://creativecommons.org/licenses/by/3.0/nz/

#################################
## Model 1: Predator-Prey
#################################

# The model assumes three states in order: empty, prey, and predator. These are states 0, 1, 2.
# Two variants are given at bottom.

#################################
##  Load library (only need to do this once!)
#################################

## If the library is not already installed, you need to run this block:
chooseCRANmirror() #choose a location near NZ. Australia works well when it isn't on fire.
install.packages("simecol")
#load the library
library(simecol)

#################################
##  Clean up workspace
#################################

remove(list=ls())

#################################
##  Set general parameters. These are constants -- they don't change value during the simulation,
#################################

gridSize <- 150        # Size of the spatial simulation in grid cells
burnIn <-  0           # Number of time steps to run before saving output
maxTime <- 2000        # Total time steps

mapsFreq <- 5           # How frequently to refresh map. 1 shows every time step 
                        # (slow but pretty), higher values will be faster but jumpier.
                        # set to any value greater than maxTime to suppress mapping.
colVec <- c("white", "#009E73", "#D55E00")   # Set colors for empty, prey, and predator.
                                      
# define growth, predation and death:
a <- 0.17               # "a" is the probability that an empty cell with all neighbouring 
                        # cells being prey will become prey (transition 0 -> 1)
b <- 0.77               # "b" is the probability that a prey cell with all neighbouring  
                        # cells being predators will become predator (transition 1 -> 2)
c <- 0.06               # "c" is the death rate of predators (transition 2 -> 0)

# Note c = 0.06, a= 0.17 and b=0.77 is from figure 10 in the source paper and gives nice 
# oscillations representing a low-birth rate prey and a long-lived, highly efficient
# predator. The variables a, b and c can have any values, but the source paper suggests
# restricting to a + b + c = 1

# define neighbourhood, in this case being the 4 touching cells. You can change this.
wdist <- matrix(c(0, 1, 0,
                  1, 1, 1,
                  0, 1, 0), nrow=3, byrow=TRUE)
                  
initials <- c(.9,0.005, 0.095)   # Starting map will have this proportion of each 
                                 # state (in order), where states are 0 = empty, 1 = prey, 2 = predator.

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
##  Define three transition functions, based on initial state (0, 1, 2 or empty, prey, predator)
#################################

##  Define transition function for cells that are currently in state 0.
transitions0 <- function(N1, maxNeigh, a)
    {
    if(length(N1)>0)                                        # Test whether any output is required
        {
        prob1 <- (N1*a)/maxNeigh                            # Define the probability of 0 -> 1 (see handout)
        prob0 <- 1-prob1                                    # Define the probability of no change (see handout)
        prob2 <- rep(0, length(N1))                         # Define the probability of 0 -> 2 (see handout)
        return(mapply(transitions, prob0, prob1, prob2)) 	# Use the function transitions with the 3 probability
                                                            # vectors to determine outcomes. Note that the order
                                                            # is important (0, 1, 2)
        } else {  #if length(N1)                                          
        return(NA)                                          # Return NA if no output required
        }  #end if length(N1)                                              
    } #end function
    
##  Define transition function for cells that are currently in state 1.
transitions1 <- function(N2, maxNeigh, b)
    {                                                       # Similar to transitions0, but with different
    if(length(N2)>0)                                        # transitions (from handout)
        {                                                   # note that prob1 can only be calculated after prob 2
        prob2 <- (N2*b)/maxNeigh                            # but they are given in the correct order in the mapply
        prob1 <- 1-prob2                                    # function
        prob0 <- rep(0, length(N2))
        return(mapply(transitions, prob0, prob1, prob2))
        } else {
        return(NA)
        }
    }
##  Define transition function for cells that are currently in state 2.
transitions2 <- function(len, c)
    {
    if(len>0)                                               # Here only length as "len"
        {                                                   # is needed, as no effect of
        prob1 <- rep(0, len)                                # neighbours.
        prob0 <- rep(c, len)
        prob2 <- rep(1-c, len)                                    
        return(mapply(transitions, prob0, prob1, prob2))
        } else {
        return(NA)
        }
    } 

#################################
##  Initialise map and define two empty variables to hold results
#################################

map0 <- matrix(sample(c(0,1,2), gridSize^2, prob=initials, 
  replace=TRUE), nrow=gridSize)
image(map0, col=colVec)

pred <- c() 
prey <- c() 

#################################
##  Run model
#################################

for(time in 1:maxTime)                                                           # Loop through time, incrementing by 1 
    {                                                                            # time step each run
    #calculate the number of neighbors of each cell type for all points in the matrix
    N1 <- neighbours(map0, state=1, wdist=wdist, bounds=1)                       # Count number of neighbours of state 1
                                                                                 # for every cell in matrix
    N2 <- neighbours(map0, state=2, wdist=wdist, bounds=1)                       # Count number of neighbours of state 2
                                                                                 # for every cell in matri
    #determine map at t+1                                                                             
    map1 <-  map0                                                                # Initialise map at t+1
    map1[map0 == 0] <- transitions0(N1[map0==0], maxNeigh = sum(wdist)-1, a)     # Determine state at t+1 for all state 0
                                                                                 # cells at time t.
    map1[map0 == 1] <- transitions1(N2[map0 == 1], maxNeigh = sum(wdist)-1, b)   # Determine t+1 state for all state 1 cells
    map1[map0 == 2] <- transitions2(length(map1[map0 == 2]),  c)                 # Determine t+1 state for all state 2 cells
    #Mapping and reporting steps:
    #After burn-in period, produce visual images and save output
    if(time > burnIn)
        {
        if(time %% mapsFreq == 0)                                                # plot a map to screen at specified frequency
            {                                                                    # using the modulo (%%) function (or remainder) 
            image(map1, col=colVec[1:(max(map1)+1)], add=TRUE)
            } #end if(time %% mapsFreq == 0)
        pred <- append(pred, sum(map1 == 2))                                     #save results into recording vectors for predators
        prey <- append(prey, sum(map1 == 1))                                     #save results into recording vectors for prey
        }  #end if(time > burnIn)
    map0 <- map1                                                                 #set map at time 0 to the new time step
    }  #end for(time)

#################################
##  Plot model output
#################################

# Plot the prey population through time in green (we set the colours above when we defined "colVec"
plot(prey, type="l", col=colVec[2], ylim=c(0, max(pred, prey)), xlab="Time", ylab="Populations")
# Add the predator data
points(pred, type="l", col=colVec[3]) 
# Add a legend.
legend(1, 0.8 * max(pred, prey), legend=c("Prey", "Predator"), col=c(colVec[2], colVec[3]), lty=1:1, cex=0.8)
# Note that in the legend function, the first two arguments (1, 0.8 * max(pred, prey)) give the x and y coordinates for 
# placement of the legend, and you may need to move to avoid overlapping any lines.

## Make a plot of how predator and prey population sizes change as a function of each other.
plot(pred ~ prey, cex=0.5, col=rainbow(length(pred)), type="b", pch=16)

# Plot how often you get different numbers of predators in the landscape.
hist(pred, breaks=1000)


#################################
## Variation 1: Habitat fragmentation
#################################

## After running the above model until it is stable, bulldoze some roads and see what happens
## We do this by adding a state "3" that has no transition probability (it always stays a 3)
## Note that this variation is run AFTER the above without resetting the map or output vectors.

colVec[4] <- "#999999"                                                           # Add road colour to colVec

nRoads <- 4                                                                      # Number of roads to add in each dimension
map0[,sample(1:gridSize,nRoads) ] <- 3                                           # Turn some columns of the map into roads
map0[sample(1:gridSize,nRoads),] <- 3                                            # Turn some rows of the map into roads

image(map0, col=colVec)                                                          # Plot the fragmented map

## Now re-run the model (starting at step "Run model") and re-do your plots (step "Plot model output")

#################################
## Variation 1.1: Why did the chicken cross the road?
#################################

# If you wish to allow dispersal across roads, you could use a larger wdist with small values 
# in the outer neighbourhood. Try changing this setting, then run model again.
# (Again, starting at step "Run model" and re-doing your plots (step "Plot model output")

wdist <- matrix(c(0,  0, .02,  0,  0,
                  0, .5,  1, .5,  0, 
                  .02, 1,  1,  1, .02, 
                  0, .5,  1, .5,  0, 
                  0,  0, .02,  0,  0), nrow=5, byrow=TRUE)
image(map0, col=colVec)     