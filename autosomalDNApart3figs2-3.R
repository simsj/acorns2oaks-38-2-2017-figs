# Mon Mar 20 08:46:13 2017 ------------------------------
# code to generate figures 2 & 3 for Genetic Genealogy: Autosomal DNA Part 3
# by James Sims
# Acorns to Oaks 38(2), 2017 manuscript in preparation

# set local directory to one linked to GitHub account
setwd("~/git/acorns2oaks-38-2-2017-figs")

# load plotting and data manipulation libraries  
library(ggplot2)
library(dplyr)

# the simulation function
getSimulationList <- function(){
        # takes integer values of cMs as integer variable or vector of integers
        # in the parent frame
        # set.seed for reproducible results
        set.seed(76107)
        for(j in 1:length(cMs)){
                # paramater set up:
                zeros <- rep.int(0,100-cMs[j])
                ones <- rep.int(1,cMs[j])
                choices <- c(zeros,ones)
                rchoices <- sample(choices)
                howmanyTimes <- 10000
                results <- c(0)
                results <- as.data.frame(results)
                results$val <- c(0)
                #Simulation of Recombination events:
                for(i in (1:howmanyTimes)){
                        counter <- 0
                        thisPick <- 0
                        while(thisPick < 1){
                                counter <- counter + 1
                                thisPick <- sample(rchoices,1, replace = TRUE)
                                if(thisPick == 1){
                                        thisData <- c(i,counter) 
                                        results <- rbind(results, thisData)     
                                }
                        }
                }
                # remove top row that is not part of results
                results <- results[2:nrow(results),] 
                # push the datafame to the variable in the parent enviornment
                resultsList[[j]] <<- as.data.frame(results)  
        }
}

##----------------- Figure 2 Simulation of 10,000 recombination events  
##----------------- each with a cM value of 6
cMs <- c(6)
# create an empty list to hold simulation results
resultsList <- list()

# call the simulation function; values in the parent environment provide
# necessary variable names and values, namely cMs and resultsList
getSimulationList()

# convert to data frame for ease of use in ggplot
aSim <- as.data.frame(resultsList[1])

# sort the simulation distribution by number of generations to recombination
# slice the sorted simulation distribution at various genealogical time depths
# use this information for dashed lines and labels
aSimSored <- arrange(aSim,val)
ninetyfive <- aSimSored$val[as.integer(nrow(aSimSored) * .95)]
ninety <- aSimSored$val[as.integer(nrow(aSimSored) * .90)]
eightyfive <- aSimSored$val[as.integer(nrow(aSimSored) * .85)]
eighty <- aSimSored$val[as.integer(nrow(aSimSored) * .80)]
fifty <- aSimSored$val[as.integer(nrow(aSimSored) * .50)]
dashedlinelabels <- paste("Recombination distribution, dashed lines\n left to right: ",fifty,", ",
                          eighty,", ",eightyfive,", ",
                          ninety,", ",ninetyfive," generations",sep="")
solidlinelabels <- paste ("Average generations to \nrecombination, solid line  = ",
                          round(mean(aSim$val),2), sep="")
# set the width of bins in the historgram
chosenBinWidth <- 2

# plot data with labels
ggplot(data = aSim, aes(aSim$val)) +
        geom_histogram(binwidth=chosenBinWidth,
                       colour = "#000000", fill= "#FFFFFF") +
        scale_x_continuous(name="Generations to Recombination") +
        scale_y_continuous(name="Count") +
        theme_bw() +
        geom_vline(xintercept = ninetyfive, size = .5, colour = "#000000",
                   linetype = "dashed") +
        geom_vline(xintercept = ninety, size = .5, colour = "#000000",
                   linetype = "dashed") +
        geom_vline(xintercept = eightyfive, size = .5, colour = "#000000",
                   linetype = "dashed") +
        geom_vline(xintercept = eighty, size = .5, colour = "#000000",
                   linetype = "dashed") +
        geom_vline(xintercept = fifty, size = .5, colour = "#000000",
                   linetype = "dashed") +
        geom_vline(xintercept = mean(aSim$val), size = 1, colour = "#000000",
                   linetype = "solid") +
        annotate("text", label = dashedlinelabels, 
                 x = 120, y = 1500, color = "black", size = 3) +
        annotate("text", label = solidlinelabels, 
                 x = 120, y = 500, color = "black", size = 3) +
        annotate("text", label = "50%", 
         x = 7, y = 1500, color = "black", size = 3) + 
        annotate("text", label = "80%", 
                 x = 22, y = 1500, color = "black", size = 3) +
        annotate("text", label = "85%", 
                 x = 34, y = 1500, color = "black", size = 3) +
        annotate("text", label = "90%", 
         x = 42, y = 1500, color = "black", size = 3) +
        annotate("text", label = "95%", 
                 x = 53, y = 1500, color = "black", size = 3) 
ggsave("figure_2.png")


##----------------- Figure 3 Simulation of 10,000 recombination events for each
##----------------- different segment cM value at various genealogical time
##----------------- depths

# set the cMs values of the DNA fragment to values of interest; integers only!
cMs <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,35)

# create an empty list to hold a results list
resultsList <- list()

# run simulations; cMs values are from parent frame set above
getSimulationList()

# Set timedepth in generations to desired values; 
# max. is 6 for autoselection of plotting symbols
genealogicalTimeDepth <- c(8,10,15,20,30,50)

#Set up dataframe to contain depth, fraglenght and survival
curves <- as.data.frame(x = 0)
colnames(curves) <- "depth"
curves$fraglength <- 0
curves$percentsurvival <- 0

# populate curves datafame
for(k in 1:length(genealogicalTimeDepth)){
        for(m in 1:length(cMs)){
                depth <- genealogicalTimeDepth[k]
                frag <- cMs[m]
                intact <- round((nrow(filter(resultsList[[m]], val > genealogicalTimeDepth[k]))/nrow(resultsList[[m]]))*100,1)
                curves <- rbind(curves,c(depth, frag,intact))
        }
}
curves <- curves[2:nrow(curves),]  # delete first row, an artifact of coding

# plot results
ggplot(curves,aes(x=fraglength,y=percentsurvival, shape=factor(depth))) +
        theme_bw() +
        geom_hline(yintercept = 5, 
                   size = .5, 
                   colour = "black",
                   linetype = "dashed") +
        geom_hline(yintercept = 10, 
                   size = .5, 
                   colour = "black",
                   linetype = "dashed") +
        geom_hline(yintercept = 15, 
                   size = .5, 
                   colour = "black",
                   linetype = "dashed") +
        geom_hline(yintercept = 20, 
                   size = .5, 
                   colour = "black",
                   linetype = "dashed") +
        geom_point() + geom_line() + 
        scale_shape_discrete(name="Time\nDepth in\nGenerations") +
        labs(x="Segment cM Value",y="Percent Survival") +
        annotate("text", label = "Recombination distribution,\n dashed lines:", 
                 x = 25, y = 40, color = "black", size = 3) +
        annotate("text", label = "80%", 
                 x = 33, y = 22, color = "black", size = 3) +
        annotate("text", label = "85%", 
                 x = 33, y = 17, color = "black", size = 3) +
        annotate("text", label = "90%", 
                 x = 33, y = 12, color = "black", size = 3)  +
        annotate("text", label = "95%", 
                 x = 33, y = 7, color = "black", size = 3) +
        annotate("text", label = "100%", 
                 x = 33, y = 1, color = "black", size = 3)  
ggsave("figure_3.png")
        



        