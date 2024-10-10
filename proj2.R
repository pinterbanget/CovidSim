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

# Loads the data.
data <- read.table("engcov.txt", header = TRUE)
data <- data[1:150,]

pearson_eval <- function(real_deaths, sim_deaths) {
  # Evaluates the fitness of the model compared to the actual data.
  # A modified Pearson value formula is used, where the difference
  # between each element, di and di_s, is squared and then divided by
  # di_s if di_s > 1, or 1 if di_s < 1; where
  # di indicates a point from the actual data, and di_s indicates a point
  # from the simulated data. The sum of all the values (the Pearson value)
  # is returned as p.
  # Note that the lower the p, the better the fitness of the model.
  #
  # Parameters:
  #   real_deaths (vec) : actual data of deaths by day
  #   sim_deaths (vec)  : simulated data of deaths by day
  #
  # Returns:
  #   p (float)         : the Pearson value calculated by the modified 
  #                       Pearson formula

  # di and di_s are declared as vectors, 
  # so vectorisation computation can be done.
  di <- real_deaths
  di_s <- sim_deaths
  p_vec <- (di - di_s)**2 / pmax(di_s, 1)
  
  # The sum of each element in the vectors is the Pearson value.
  p <- sum(p_vec)

  return(p)
}

deconv <- function(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL,
                   meanlog = 3.152, sdlog = 0.451) {
  # Computes the estimated COVID-19 infection date for patients, given
  # the death date for said patients.
  # ____
  # TODO: add more thorough comment
  #
  # Parameters:
  #     t (vec)         :
  #     deaths (vec)    : days of death
  #     n.rep (int)     : the number of iterations to be done.
  #     bs (bool)       : if bool
  #     t0 (vec)        : a vector of initial guesses for the day of infection
  #     meanlog (float) : the mean of the log normal distribution.
  #                       Defaults to 3.152 in line with the ISARIC study.
  #     sdlog (float)   : the deviation standard of the log normal distribution.
  #                       Defaults to 0.451 in line with the ISARIC study.
  #
  # Returns:
  #     P (vec)         : an n.rep-sized vector containing the history of the
  #                       Pearson values for each iteration
  #     inft (mat)      : a matrix of shape (310, n.rep) containing deaths
  #                       by t0 for each iteration
  #     t0 (vec)        : a vector of infection day estimation for each patient

  # Calculates the total number of deaths occurred from the data.
  n <- sum(deaths)

  # Creates an n-spaced vector, with each element of the vector representing
  # a patient, to store the day of the year said patient died, e.g.
  # number 1 on the 1st element means the 1st patient died on January 1;
  # number 32 on the 2nd element means the 2nd patient died on February 1, etc.
  death_days <- rep(t, deaths)
  
  # Creates a 310-spaced vector to store the total deaths occurred on each day.
  # For this simulation, a limit is set where the deaths can only occur
  # between day 1 (inclusive) and day 310 (inclusive) of the year.
  total_deaths_by_day <- tabulate(death_days, nbins = 310)
  
  # If bootstrapping is turned on, the total_deaths_by_day data is changed
  # to reflect Poisson distribution ____ 
  # TODO: complete this comment
  if (bs) {
    total_deaths_by_day <- rpois(total_deaths_by_day,lambda=total_deaths_by_day)
  }

  # Generates a probability of possible infection-to-death time for COVID-19
  # patients (the number of days between when a COVID-19 patient contracted 
  # COVID-19 and their death).
  infection_to_death <- dlnorm(1:80, meanlog, sdlog)
  inf_to_d_normalised <- infection_to_death / sum(infection_to_death)

  # If not provided, generates t0, an n-spaced vector to estimate the
  # day of COVID-19 infection for each patient.
  # Elements in t0 are gained from subtracting the days of patient deaths
  # by a randomised number from 1 to 80, picked using inf_to_d_normalised
  # as the weights for each number.
  if (is.null(t0)) {
    t0 <- death_days - sample(1:80, n, prob = inf_to_d_normalised, 
                              replace = TRUE)
  }
  
  # Creates a matrix of shape (310, n.rep) to store the estimated
  # 
  # for every iteration from 1 to n.rep.
  inft <- matrix(data = NA, nrow = 310, ncol = n.rep)

  # Creates an n.rep-spaced vector to store the Pearson score
  # after each iteration.
  P <- rep(0, n.rep)
  
  # Defines three steps to pick for changing t0 during the iterations.
  # The idea is when the iterations are 0-50% complete, a wider range
  # of steps (step_options) can be picked. Then, when it is more than
  # 50% done, a smaller range of steps (step_options_50) can be picked.
  # The range of steps is narrowed even more (step_options_75) when the
  # iteration is 75% done.
  step_options <- c(-8, -4, -2, -1, 1, 2, 4, 8)
  step_options_50 <- c(-4, -2, -1, 1, 2, 4)
  step_options_75 <- c(-2, -1, 1, 2)
  
  # Loops over n.rep times.
  for (i in 1:n.rep) {
    # The first part of this loop is to generate the variables that will
    # be used, together with evaluating the Pearson score
    # 
    sim_inf_dur <- sample(1:80, n, prob = inf_to_d_normalised, replace = TRUE)
    sim_death_days <- t0 + sim_inf_dur
    sim_death_days <- pmax(sim_death_days, 1) # done to filter out values < 1
    sim_death_days <- pmin(sim_death_days, 310) # done to filter out values > 310
    total_sim_deaths_by_day <- tabulate(sim_death_days, nbins = 310)
    pearson_value <- pearson_eval(total_deaths_by_day, total_sim_deaths_by_day)
    
    # The second part of this loop is to shift sim_death_days by a random step
    # to try to improve the Pearson value. The shifting and evaluation is done
    # per-element, so another loop of n times will be done in this loop.
    # Later on in the code, sim_death_days will be subtracted by sim_inf_dur
    # to get the new t0 which would fit the data better.
    
    # Defines which step distribution to sample based on the
    # current iteration number.
    if (i <= round(n.rep * 0.5)) {
      step_to_sample <- step_options
    } else if (i <= round(n.rep * 0.75)) {
      step_to_sample <- step_options_50
    } else {
      step_to_sample <- step_options_75
    }
    
    # Samples the step n times.
    random_steps <- sample(step_to_sample, n, replace = TRUE)
    
    # Starts the loop to shift and evaluate sim_death_days.
    # Each element in sim_death_days is accessed and updated in a random order
    # to avoid ______
    # TODO: complete desc above
    for (index in sample(1:n)) {
      # Saves the old day of death value.
      old_death_day <- sim_death_days[index]
      
      # Generates the new day of death value, with filtering to make sure
      # it is within the (1, 310) range.
      new_death_day <- sim_death_days[index] + random_steps[index]
      new_death_day <- max(new_death_day, 1)
      new_death_day <- min(new_death_day, 310)
      
      # Copies the total_sim_deaths_by_day tabulation 
      # to be modified according to the new day of death value.
      new_total_sdbd <- total_sim_deaths_by_day
      new_total_sdbd[new_death_day] <- new_total_sdbd[new_death_day] + 1
      new_total_sdbd[old_death_day] <- new_total_sdbd[old_death_day] - 1
      
      # Evaluates the Pearson value with the new death day values.
      current_pv <- pearson_eval(total_deaths_by_day, new_total_sdbd)
      
      # If the Pearson value decreases, 
      # update the new variables made within this loop.
      if (current_pv < pearson_value) {
        pearson_value <- current_pv
        sim_death_days[index] <- new_death_day
        total_sim_deaths_by_day <- new_total_sdbd
      }
    }
    
    # The third part of this loop is to update all variables
    # ___
    # TODO: complete above comment
    
    t0 <- sim_death_days - sim_inf_dur
    t0_freq <- tabulate(t0, nbins = 310)
    inft[, i] <- t0_freq
    
    # Puts the final Pearson score for the i-th iteration to P.
    P[i] <- pearson_value
    
    # Generates 3 line charts for this iteration, with the x axis
    # being the possible days of data (in this case 1:310),
    # and the y axis being the number of people.
    # The first line is the estimated line of COVID-19 infection,
    # the second line is the actual data of COVID-19 deaths, and
    # the third line is the estimated data of COVID-19 deaths.
    matplot(1:310, 
            cbind(t0_freq, total_deaths_by_day, total_sim_deaths_by_day),
            )
    
    # TODO: "prettify" plot above to include grids, titles, fixed axis scale,
    # axis labels, etc.
  }
  
  return(list(P, inft, t0))
}

t <- data$julian
deaths <- data$nhs
t0 <- deconv(t = t, deaths = deaths, n.rep = 100)[[3]]
t0_pearson <- deconv(t = t, deaths = data$nhs, n.rep = 100, t0=t0, bs=TRUE)[[3]]

# TODO: create final graphs
# TODO: check bootstrap implementation, is it correct already?