
# Spatial simulation code for University of Canterbury Biol377.
# Pathogen Dynamics Simulation : Two prey species variation
# Intended to be used in the context of this document: https://docs.google.com/document/d/1rkp9M2OBmIHHmQu4y52eutd-0Ylll0utAj7guDvg8XI/edit?usp=sharing

# I. Dickie,  J. Tylianakis, D. Stouffer & M. Turnbull, School of Biological Sciences, University of Canterbury NZ
# CC-BY license: https://creativecommons.org/licenses/by/3.0/nz/

# This variant is basically the same as the Predator Prey model, except that there are now two
# prey species. This requires rewriting the model with 4 possible states (empty, prey1, prey2, predator)

#################################
##  Load library
#################################

## If the library is not already installed, you need to run the next 2 lines
chooseCRANmirror() #choose a location near NZ. Australia works well when it isn't on fire.
install.packages("simecol")

require(simecol)

#################################
##  Clean up workspace
#################################

remove(list=ls())
     
#################################
##  Set general parameters. These are constants -- they don't change value during the simulation,
#################################

gridSize <- 100          # Size of the spatial simulation in grid cells
burnIn <-  1000           # Number of time steps to run before saving output
maxTime <- 2000          # Total time steps

mapsFreq <- 20           # How frequently to refresh map. 1 shows every time step 
                         # (slow but pretty), higher values will be faster but jumpier.
                         # set to any value greater than maxTime to suppress mapping.
                        
colVec <- c("white", "#009E73", "#0072B2", "#D55E00")

## define growth, predation and death:
c <- 0.06   #3 based death rate of predators
a1 <- 0.17  ## "a1" is the probability that an empty cell with all neighbouring prey1 will become prey1
b <- 0.2   ## "b" is the probability that a prey cell with all neigbouring predators will become predator

# define neighbourhood, in this case being the 4 touching cells. You can change this.
wdist <- matrix(c(0, 1, 0,
                  1, 1, 1,
                  0, 1, 0), nrow=3, byrow=TRUE)
                  
initials <- c(.9,0.0025, 0.0025, 0.095)   ## Starting map will have this proportion of each state (in order)


#################################
##  Set habitat quality parameters
#################################

goodHabitatSize <- round(gridSize / 4)  ## How big a patch of "good quality" habitat to include

qualFactor <- 4      ## This sets the relative quality of the two habitats for prey 2. Set to 1 for no difference.
a2 <- a1/qualFactor  ## "a2" is the probability that an empty cell with all neighbouring prey2 will become prey2
                     ## I've modified this to by less than a1 by 1/qualFactor
a2goodHabitat <- a1*qualFactor - a2  ## This is the probability that an empty cell with 4 neighbouring prey2 will 
                    ## become prey2 in a good quality patch. I've made this qualFactor times greater than a1.
                    ## the '- a2' is just because of the way the math works in the function.
                                            


#################################
## Define functions. Note that these are very similar to earlier, but we have added the second prey.
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
    
transitions0 <- function(N1, N2, habitat, maxNeigh, a1, a2, a2goodHabitat) 
    {
    if(length(N1)>0)             
        {  
        prob1 <- (N1*a1)/maxNeigh     
        prob2 <- (N2*(a2+habitat*a2goodHabitat))/maxNeigh                    
        prob0 <- 1-(prob1 + prob2)    ##                              
        prob3 <- rep(0, length(N1)) 
        return(mapply(transitions, prob0, prob1, prob2, prob3))
        } else {
        return(NA)
        }  
    }
transitions1 <- function(N3, maxNeigh, b)  # Function for prey 1 becoming predator is not changed much
    {
    if(length(N3)>0)             
        {  
        prob3 <- (N3*b)/maxNeigh                           
        prob1 <- 1-prob3                                    
        prob0 <- rep(0, length(N3))
        prob2 <- rep(0, length(N3))
        return(mapply(transitions, prob0, prob1, prob2, prob3))
        } else {
        return(NA)
        }     
    }
transitions2 <- function(N3, maxNeigh, b)    # Function for prey 2 becoming predator is almost the same as for prey 1
    {
    if(length(N3)>0)             
        {  
        prob3 <- (N3*b)/maxNeigh                           
        prob2 <- 1-prob3                                    
        prob0 <- rep(0, length(N3))
        prob1 <- rep(0, length(N3))
        return(mapply(transitions, prob0, prob1, prob2, prob3))
        } else {
        return(NA)
        }     
    }
transitions3 <- function(len, c) 		# Function for predator death is not changed much
    {
    if(length(N1)>0)             
        {   
        prob0 <- rep(c, len)
        prob1 <- rep(0, len)  
        prob2 <- rep(0, len)    
        prob3 <- rep(1-c, len)                                    
        return(mapply(transitions, prob0, prob1, prob2, prob3))
        } else {
        return(NA)
        }  
    } 
             
#################################
##  define habitat heterogeneity
#################################


habitatQuality <- matrix(0, nrow=gridSize, ncol=gridSize)
improved <- (gridSize/2 - goodHabitatSize/2):(gridSize/2 + goodHabitatSize/2)
habitatQuality[improved, improved] <- 1

#can visualise: 
image(habitatQuality, col=c("yellow", "green"), main = "Habitat quality map")
                                   
#################################
##  Initialise map and define empty variables to hold results
#################################
      
map0 <- matrix(sample(c(0,1,2, 3), gridSize^2, prob=initials, replace=TRUE), nrow=gridSize)
image(map0, col=colVec)
pred <- prey1 <- prey2 <- c()                             # set up "recording" vectors for results


#################################
##  Run model -- note that this is very similar to earlier, except we have added one more possibles state.    
#################################
          
for(time in 1:maxTime)                                  
    {
    N1 <- neighbours(map0, state=1, wdist=wdist, bounds=1)
    N2 <- neighbours(map0, state=2, wdist=wdist, bounds=1)
    N3 <- neighbours(map0, state=3, wdist=wdist, bounds=1)
    map1 <- map0                                        # 
    map1[map0 == 0] <- transitions0(N1[map0==0], N2[map0==0],habitatQuality[map0==0], 
        maxNeigh = sum(wdist), a1, a2, a2goodHabitat)
    map1[map0 == 1] <- transitions1(N3[map0 == 1], maxNeigh = sum(wdist)-1, b)
    map1[map0 == 2] <- transitions2(N3[map0 == 2], maxNeigh = sum(wdist)-1, b)
    map1[map0 == 3] <- transitions3(length(map1[map0 == 3]),  c)
    #After burn-in period, produce visual images (if maps == TRUE) and save output
    if(time > burnIn)
        {
        if(time %% mapsFreq == 0)
            {
            image(map1, col=colVec, add=TRUE)
            } #end if(time %% mapsFreq == 0)
        pred <- append(pred, sum(map1 == 3))
        prey1 <- append(prey1, sum(map1 == 1))
        prey2 <- append(prey2, sum(map1 == 2))
        }  #end if(time > burnIn)
    map0 <- map1 
    }  #end for(time) <- neighbours(map0, state=1, wdist=wdist, bounds=1)


#################################
##  Plot model output -- Adding a line for the second prey
#################################

plot(prey1, type="l", col=colVec[2], ylim=c(0, max(pred, prey1, prey2)), xlab="Time", ylab="Populations")
points(pred, type="l", col=colVec[4])
points(prey2, type="l", col=colVec[3])
legend(1, 0.8 * max(pred, prey1, prey2), legend=c("Prey 1", "Prey 2", "Predator"), col=c(colVec[2], colVec[3], colVec[4]), lty=1:1, cex=0.8)

#################################
##  Saving results
#################################

## For all the spatial models, you will want to do multiple runs and compare results with different settings.
## the "remove(list=ls())" line in the code removes all saved information before each run, so we need to output results
## to a file and then re-load them for comparison graphs.

## Say that I decide to save prey1, prey2, and predator population data.  I can do this by writing to a csv file
## I've commented out this line to avoid accidentally over-writing a file:
#
write.csv(file="modelOutput_qualFactor4_r1.csv", data.frame(time=1:length(prey1), prey1=prey1, prey2=prey2, pred=pred), row.names=FALSE)

## you should also save the results as an R object. This is VERY useful if you are trying to fix bugs, or forget
## some critical detail.  All you have to do is load('file name') and it will bring back the run:
save(list=ls(), file="Model_Run_qualFactor4_r1.rdat")


#################################
##  Comparing results
#################################

##I also ran the model with no difference in substrate quality and then saved results:
#write.csv(file="modelOutput_qualFactor1_r1.csv", data.frame(time=1:length(prey1), prey1=prey1, prey2=prey2, pred=pred), row.names=FALSE)

##Then I can read back in the file and compare outcomes:
output1 <- read.csv(file="modelOutput_qualFactor4_r1.csv", header=TRUE)
output1$runName <- "QualDiff4_r1"
output2 <- read.csv(file="modelOutput_qualFactor1_r1.csv", header=TRUE)
output2$runName <- "QualDiff1_r1"
output <- rbind(output1, output2)

output$totalPrey <- output$prey1 + output$prey2
plot(output$totalPrey[output$runName == "QualDiff4_r1"], type="l", lty=1)
points(output$totalPrey[output$runName == "QualDiff1_r1"], col="purple", lty=2)

plot(pred ~ time, type="l", col="darkgreen", data=output[output$runName == "QualDiff4_r1",])
points(output$pred[output$runName == "QualDiff1_r1"], col="purple", type="l")

## The graphs are pretty, but it might be just as effective to simply report the mean values.
## Remember that replication is important, even in a simulation. You can run each set of parameters multiple times and save multiple versions.