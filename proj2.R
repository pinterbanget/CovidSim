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
setwd("/Users/josephgill/covidsim") # Joseph's path
# setwd("/Users/fransiskusbudi/uoe/stat_prog/covidsim") # Frans' path

# Data Loading
data <- read.table("engcov.txt", header = TRUE)
data <- data[1:150,]
days <- data$julian
deaths <- data$nhs

meanlog <- 3.152; sdlog <- 0.451

infection_to_death <- dlnorm(1:80, meanlog, sdlog)
infection_to_death_normalized <- infection_to_death / sum(infection_to_death)


sum(infection_to_death)

plot(infection_to_death,type="l", col='red')

n <- sum(deaths) #29422
death_day <- rep(days,deaths)


infection_duration <- sample(1:80,n,prob=infection_to_death_normalized, replace = TRUE)
t0 <- death_day - infection_duration 

print(tabulate(death_day)) #how many people died on each day

for (i in 1:n){
    infection_duration <- sample(1:80,n,prob=infection_to_death_normalized, replace = TRUE)
    t0 <- death_day - infection_duration 
}

print(t0)

deconv<- function(t,deaths,n.rep=100,bs=FALSE,t0=NULL){

}


# below is function to calculate p
pearson_eval <- function(actual_deaths, simulated_deaths) {
  # Evaluates the fitness of the model compared to the actual data.
  # Parameters:
  #   actual_deaths (vector): actual data containing actual deaths
  #   simulated_deaths (vector): simulated data
  # Returns:
  #   p (num): the final error calculated by the modified Pearson formula
  
  # Sets p to store the sum of values by elements.
  p <- 0
  
  # Loops through all data.
  for (i in 1:length(actual_deaths)) {
    di <- actual_deaths[i]
    di_s <- simulated_deaths[i]
    p <- p + ((di - di_s)**2 / max(c(1, di_s)))
  }
  
  # Returns the sum of it all
  return(p)
}

print(pearson_eval(deaths,t0))