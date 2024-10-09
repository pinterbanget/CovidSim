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
setwd("/Users/rj/Documents/Codes/StatProg/covidsim") # Ryan's path
# setwd("/Users/josephgill/covidsim") # Joseph's path
# setwd("/Users/fransiskusbudi/uoe/stat_prog/covidsim") # Frans' path

# Data Loading
data <- read.table("engcov.txt", header = TRUE)
data <- data[1:150,]
days <- data$julian
deaths <- data$nhs

meanlog <- 3.152
sdlog <- 0.451

infection_to_death <- dlnorm(1:80, meanlog, sdlog)
infection_to_death_normalized <- infection_to_death / sum(infection_to_death)

n <- sum(deaths) #29422
death_day <- rep(days,deaths)


infection_duration <- sample(1:80,n,prob=infection_to_death_normalized, replace = TRUE)
t0 <- death_day - infection_duration 

tabby <- tabulate(death_day, nbins = 310) #how many people died on each day


pearson_eval <- function(actual_deaths, simulated_deaths) {
  # Evaluates the fitness of the model compared to the actual data.
  # Parameters:
  #   actual_deaths (vector): actual data containing actual deaths
  #   simulated_deaths (vector): simulated data
  # Returns:
  #   p (num): the final error calculated by the modified Pearson formula

  # Sets p to store the sum of values by elements.
  di <- actual_deaths
  di_s <- simulated_deaths
  pearson_score <- sum((di - di_s)**2 / max(c(1, di_s)))

  # Returns the sum of it all
  return(pearson_score)
}

deconv <- function(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL,
                   meanlog = 3.152, sdlog = 0.451) {
  # blablal
  # Params:
  #     t (vec)         : 
  #     deaths (vec)    :
  #     n.rep (int)     :
  #     bs (bool)       :
  #     t0 (vec)        :
  #     meanlog (float) :
  #     sdlog (float)   :
  #
  # Returns:
  #     
  #     

  # Calculates the number of deaths occurred from the data.
  n <- sum(deaths)

  # Creates an n-spaced vector to store in which day of the year
  # each patient died. E.g.: 1 means the patient died on January 1,
  # 32 means the patient died on February 1.
  death_days <- rep(t, deaths)

  # Creates a 310-spaced vector to store the total deaths occurred on each day.
  total_deaths_per_day <- tabulate(death_days, nbins = 310)

  # Generates a probability of possible infection-to-death time for COVID-19
  # patients (the number of days between when a COVID-19 patient contracted 
  # COVID-19 until their death).
  infection_to_death <- dlnorm(1:80, meanlog, sdlog)
  normalized_inf_to_d <- infection_to_death / sum(infection_to_death)

  # Generates a distribution of possible duration that a COVID-19 patient
  # contracted COVID-19 until their death.
  infection_duration <- sample(1:80, n, prob = normalized_inf_to_d, 
                               replace = TRUE)

  # Generates t0, if not provided.
  if (is.null(t0)) {
    t0 <- death_days - infection_duration
  }

  inft <- matrix(data = NA, nrow = 310, ncol = n.rep)

  # Creates an n.rep-spaced vector to store the Pearson score after each epoch
  P <- rep(0, n.rep)

  to_pick_1 <- c(-8, -4, -2, -1, 1, 2, 4, 8)
  to_pick_2 <- c(-4, -2, -1, 1, 2, 4)
  to_pick_3 <- c(-2, -1, 1, 2)
  
  # sim_inf_days <- sample(1:80, n, prob=normalized_inf_to_d, replace = TRUE)
  # sim_death_days <- t0 + sim_inf_days
  # sim_death_days <- pmax(sim_death_days, 1) # done to filter out values < 1
  # sim_death_days <- pmin(sim_death_days, 310) # done to filter out values > 310
  # total_sim_deaths_per_day <- tabulate(sim_death_days, nbins = 310)
  # pearson_score <- pearson_eval(total_deaths_per_day, total_sim_deaths_per_day)

  for (i in 1:n.rep) {
    sim_inf_days <- sample(1:80, n, prob=normalized_inf_to_d, replace = TRUE)
    sim_death_days <- t0 + sim_inf_days
    sim_death_days <- pmax(sim_death_days, 1) # done to filter out values < 1
    sim_death_days <- pmin(sim_death_days, 310) # done to filter out values > 310
    total_sim_deaths_per_day <- tabulate(sim_death_days, nbins = 310)
    pearson_score <- pearson_eval(total_deaths_per_day, total_sim_deaths_per_day)

    if (i < round(n.rep*0.5)) {
      to_sample <- to_pick_1
    } else if (i < round(n.rep*0.75)) {
      to_sample <- to_pick_2
    } else {
      to_sample <- to_pick_3
    }

    random_indices <- sample(1:length(t0), n, replace = FALSE)
    random_additions <- sample(to_sample, n, replace = TRUE)

    for (j in 1:length(t0)) {
      current_index <- random_indices[j]

      old_death_day <- sim_death_days[current_index]
      new_death_day <- sim_death_days[current_index] + random_additions[current_index]

      new_death_day <- max(new_death_day, 1)
      new_death_day <- min(new_death_day, 310)

      total_sim_deaths_per_day[new_death_day] <- total_sim_deaths_per_day[new_death_day] + 1
      total_sim_deaths_per_day[old_death_day] <- total_sim_deaths_per_day[old_death_day] - 1

      current_ps <- pearson_eval(total_deaths_per_day, total_sim_deaths_per_day)

      if (current_ps < pearson_score) {
        pearson_score <- current_ps
      } else {
        sim_death_days[current_index] <- old_death_day
        total_sim_deaths_per_day[new_death_day] <- total_sim_deaths_per_day[new_death_day] - 1
        total_sim_deaths_per_day[old_death_day] <- total_sim_deaths_per_day[old_death_day] + 1
      }
    }

    P[i] <- pearson_score
    t0 <- sim_death_days - sim_inf_days
    t0_freq <- tabulate(t0, nbins = 310)
    inft[, i] <- t0_freq

    matplot(1:310, cbind(t0_freq, total_deaths_per_day, total_sim_deaths_per_day))
  }
  
  return(t0)
}


deconv(t = data$julian, deaths = data$nhs, n.rep = 100)