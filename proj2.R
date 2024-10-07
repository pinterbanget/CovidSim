# Authors:
# Ryanson Jonathan (s2570340@ed.ac.uk)
# Joseph Gill (s1910643@ed.ac.uk)
# Fransiskus Budi Kurnia Agung (s2670828@ed.ac.uk)

# Contributions:


##################################
##################################

# Overview:
# This R code aims to implement a simple simulation method to infer fatal 
# incidence rates from Covid deaths in English hospitals

##################################
##################################
# ***
# Sets the working directories for the coders.
# setwd("/Users/rj/Documents/Codes/StatProg/covidsim") # Ryan's path
# setwd("/Users/josephgill/Documents/covidsim") # Joseph's path
setwd("/Users/fransiskusbudi/uoe/stat_prog/covidsim") # Frans' path

# Data Loading
data <- read.table("engcov.txt", header = TRUE)
data <- data[1:150,]
days <- data$julian
deaths <- data$nhs

meanlog <- 3.152; sdlog <- 0.451

infection_to_death <- dlnorm(1:80, meanlog, sdlog)
infection_to_death_normalized <- infection_to_death / sum(infection_to_death)

#how to normalize?

sum(infection_to_death)

plot(infection_to_death,type="l", col='red')

n <- sum(deaths) #29422
death_day <- rep(days,deaths)


infection_duration <- sample(1:80,n,prob=infection_to_death_normalized, replace = TRUE)
t0 <- death_day - infection_duration 
