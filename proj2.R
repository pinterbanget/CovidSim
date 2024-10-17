## Authors:
# Ryanson Jonathan (s2570340@ed.ac.uk)
# Joseph Gill (s1910643@ed.ac.uk)
# Fransiskus Budi Kurnia Agung (s2670828@ed.ac.uk)

## Contributions:
# Ryan (xx%) - pearson_eval function, initial deconv function, commenting
# Joseph (xx%) - 
# Frans (xx%) - 

##################################
##################################

## Overview:
# This R code contains functions to implement a simple
# simulation method to infer fatal COVID-19 incidence rates
# using data from COVID-19 deaths in English hospitals.
# The idea is since there is sufficient information on distribution
# of time from infection to death from COVID-19, we can estimate the
# days that patients first contracted COVID-19. The result from this
# simulation can be used to determine the effectiveness of control methods.

# The simulation works by first getting NHS data of the number of patients
# that died each day and on what day they died. Then, using infection-to-death
# time distribution, we can estimate the days that patients most likely got
# COVID-19 by subtracting the death days by a random time from the distribution.
# This fatal infection day estimation is then added by a resampled time from the
# infection-to-death distribution (making them effectively death days)
# and comparing them with the actual data of death days to judge the
# goodness of the estimation (done using a modified Pearson formula).
# This process is done multiple times to get an estimation of infection days
# that fits most with the distribution data.

# When the data converges, the simulation is run again, this time to get
# the sense of uncertainty from the estimation. This is done by sampling
# using a Poisson distrubution where the mean is the deaths each day from
# the NHS data. Since the data is very large, a Poisson distrubution will
# be a good approximation.
#
# TODO: still pretty unclear IMO, thoughts?

##################################
##################################
# Sets the working directories for the coders.
# setwd("/Users/rj/Documents/Codes/StatProg/covidsim") # Ryan's path
# setwd("/Users/josephgill/covidsim") # Joseph's path
setwd("/Users/fransiskusbudi/uoe/stat_prog/covidsim") # Frans' path


pearson_eval <- function(real_deaths, sim_deaths) {
  # Evaluates the fitness of the model compared to the actual data
  # using a modified Pearson value formula: the difference between
  # each element, di and di_s, is squared, and then divided by di_s
  # (or 1 if di_s is lower than 1), where di indicates a point from the
  # actual data, and di_s indicates a point from the simulated data.
  # The sum of all the points differences is the Pearson value.
  #
  # Parameters:
  #   real_deaths (vec) : actual data of deaths by day.
  #   sim_deaths (vec)  : simulated data of deaths by day.
  #
  # Returns:
  #   p (float)         : the Pearson value calculated by the modified
  #                       Pearson formula.

  # di and di_s are declared as vectors,
  # so vectorisation computation can be done.
  di <- real_deaths
  di_s <- sim_deaths
  pearson_val_vec <- (di - di_s)**2 / pmax(di_s, 1)

  # The sum of each element in the vectors is the Pearson value.
  pearson_val <- sum(pearson_val_vec)

  return(pearson_val)
}

deconv <- function(t, deaths, n.rep = 100, bs = FALSE, t0 = NULL) {
  # Computes the estimated COVID-19 infection date for patients,
  # given the death date for said patients. This is done by ____
  # TODO: complete this comment. especially about bootstrapping
  #
  # Parameters:
  #     t (vec)         : a vector of the days of the year.
  #     deaths (vec)    : a vector of number of deaths happening on
  #                       the days supplied by t.
  #     n.rep (int)     : the number of iterations to be done
  #                       (default: 100).
  #     bs (bool)       : a flag indicating bootstrapping application
  #                       (default: FALSE).
  #     t0 (vec)        : a vector of estimations for the days of infections
  #                       (default: NULL).
  #
  # Returns:
  #     P (vec)         : an n.rep-sized vector containing the history of the
  #                       Pearson values for each iteration.
  #     inft (mat)      : a matrix of shape (310, n.rep) containing deaths
  #                       by t0 for each iteration.
  #     t0 (vec)        : a vector of estimations for the days of infections.
  #     total_sdbd (vec): a vector of simulated days of deaths from the model.
  #

  # Calculates the total number of deaths occurred from the data.
  n <- sum(deaths)

  # Creates an n-spaced vector, with each element of the vector representing
  # a patient, to store the day of the year said patient died.
  death_days <- rep(t, deaths)

  # Creates a 310-spaced vector to store the total deaths occurred on each day.
  # For this simulation, a limit is set where the deaths can only occur
  # between day 1 (inclusive) and day 310 (inclusive) of the year.
  total_deaths_by_day <- tabulate(death_days, nbins = 310)

  # Generates a probability of possible infection-to-death time for COVID-19
  # patients. The log normal distribution is used here, with mean = 3.152
  # and sd = 0.451, in line with the ISARIC study of COVID-19 distribution.
  inf_to_d_dist <- dlnorm(1:80, meanlog = 3.152, sdlog = 0.451)

  # This distribution is then normalised.
  inf_to_d_dist <- inf_to_d_dist / sum(inf_to_d_dist)

  # If not provided, generates t0, an n-spaced vector to estimate the
  # day of COVID-19 infection for each patient. Elements in t0 are gained
  # from subtracting death_days by a randomised number from 1 to 80,
  # picked using inf_to_d_dist as the weights for each number.
  if (is.null(t0)) {
    t0 <- death_days - sample(1:80, n, prob = inf_to_d_dist,
                              replace = TRUE)
  }

  # Creates a matrix of shape (310, n.rep) to store the estimated t0
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
    # The first part of this loop is to generate the variables that will be
    # used to update t0, together with evaluating the initial Pearson score.

    # If bootstrapping is enabled, the number of deaths by day is changed to
    # samples from a Poisson distribution with the mean given by the real data.
    # If not, the number of deaths by day used for the model is the real data.
    if (bs) {
      # Generates the number of deaths by day using Poisson distribution.
      death_pois <- rpois(length(deaths), lambda = deaths)

      # Generates an n-spaced vector of sampled death days for each patient.
      death_days_pois <- rep(t, death_pois)

      # Tabulates the vector to get the total number of (Poisson) deaths by day.
      total_deaths_by_day <- tabulate(death_days_pois, nbins = 310)
    }

    # Samples random infection-to-death times for each patient.
    sim_inf_dur <- sample(1:80, n, prob = inf_to_d_dist, replace = TRUE)

    # Adds the times to t0 to obtain simulated days of death.
    sim_death_days <- t0 + sim_inf_dur

    # Tabulates the result to get the total simulated deaths by day.
    total_sim_deaths_by_day <- tabulate(sim_death_days, nbins = 310)

    # Computes the initial Pearson value for this iteration.
    pearson_value <- pearson_eval(total_deaths_by_day, total_sim_deaths_by_day)

    # The second part of this loop is to shift elements in t0 by random steps
    # to try to improve the Pearson value. Because the shifting and evaluation
    # is done per-element, another n-time loop will be done inside this loop.

    # Defines which step distribution to sample based on the
    # current iteration number.
    if (i <= round(n.rep * 0.5)) {
      step_to_sample <- step_options
    } else if (i <= round(n.rep * 0.75)) {
      step_to_sample <- step_options_50
    } else {
      step_to_sample <- step_options_75
    }

    # Creates an n-spaced vector containing samples of steps.
    random_steps <- sample(step_to_sample, n, replace = TRUE)

    # Starts the loop to shift and evaluate elements in t0.
    # Each element in t0 is accessed and updated in a random order
    # to avoid biases.
    for (index in sample(1:n)) {
      # Saves the old day of infection value.
      old_infected_day <- t0[index]

      # Gets a new day of infection value by adding a step from random_steps,
      # with filtering to make sure it is within the (1, 310) range.
      new_infected_day <- old_infected_day + random_steps[index]
      new_infected_day <- max(new_infected_day, 1)
      new_infected_day <- min(new_infected_day, 310)

      # To calculate the Pearson value, the day of death values are needed.
      # Gets the old day of death value, together with the new day of death.
      old_death_day <- sim_death_days[index]
      new_death_day <- new_infected_day + sim_inf_dur[index]

      # Copies the total simulated deaths by day tabulation
      # and modifies according to the new day of death values.
      new_total_sdbd <- total_sim_deaths_by_day
      new_total_sdbd[new_death_day] <- new_total_sdbd[new_death_day] + 1
      new_total_sdbd[old_death_day] <- new_total_sdbd[old_death_day] - 1

      # Evaluates the Pearson value with the new day of death values.
      new_pv <- pearson_eval(total_deaths_by_day, new_total_sdbd)

      # If the Pearson value decreases,
      # updates the new variables made within this loop.
      if (new_pv < pearson_value) {
        pearson_value <- new_pv
        t0[index] <- new_infected_day
        sim_death_days[index] <- new_death_day
        total_sim_deaths_by_day <- new_total_sdbd
      }
    }

    # After t0 has converged, t0 is then tabulated,
    # to get the number of new fatal infections by day.
    t0_freq <- tabulate(t0, nbins = 310)

    # This iteration's number of new infections by day is stored in inft.
    inft[, i] <- t0_freq

    # Puts the final Pearson score for the i-th iteration to P.
    P[i] <- pearson_value

    # Generates three line charts for this iteration, with the x axis
    # being the possible days of data (in this case 1:310),
    # and the y axis being the number of people.

    # A label is produced to indicate if bootstrapping was done or not.
    if (bs) {
      bs_label <- "(Bootstrapped)"
    } else {
      bs_label <- ""
    }

    # The first generated ine is the estimated data of COVID-19 infection.
    plot(1:310, t0_freq, col = "black", type = "l", lwd = 1,
         xlab = "Amount of days since 1st January 2020",
         ylab = "Num. of people", ylim = c(0, 1800),
         main = paste("COVID-19 Infections and Deaths Simulation ",
                      bs_label, "\n(iter #", i, ")", sep = ""))

    # The second line is the actual data of COVID-19 deaths.
    lines(1:310, total_deaths_by_day, col = "blue",
          type = "l", lwd = 1, lty = "solid")

    # The third line is the estimated data of COVID-19 deaths.
    lines(1:310, total_sim_deaths_by_day, col = "red",
          type = "l", lwd = 1, lty = "dotdash")

    # A line is shown to indicate the first day of lockdown (day 84).
    abline(v = 84, col = "red", lwd = 2)

    # A legend is put on the top right of the graph.
    legend(x = "topright", inset = 0.05,
           legend = c("Est. Incidence Trajectory", "Est. Deaths",
                      "Actual Deaths", "UK Lockdown Start (Day 84)"),
           col = c("black", "red", "blue", "red"),
           lty = c("solid", "dotdash", "solid", "solid"),
           lwd = c(1, 1, 1, 2), cex = 0.8)
  }

  # Saves the simulation of death days from the final iteration.
  # This is necessary to know the fit of the simulation model
  # vs. the real data.
  total_sdbd <- total_sim_deaths_by_day

  return(list(P, inft, t0, total_sdbd))
}


# Loads the data.
data <- read.table("engcov.txt", header = TRUE)

# Takes the first 150 rows of the data.
data <- data[1:150, ]

# Gets t, the days of the year, from the data.
t <- data$julian

# Gets deaths, the number of deaths happening by day, from the data.
deaths <- data$nhs
death_days <- rep(t, deaths)

# Calls the deconv function to run the initial simulation for estimated
# infection days, just by feeding t and deaths to the function.
initial_sim <- deconv(t, deaths)

# To get an idea of uncertainty within the estimated infection days,
# calls the deconv function again, this time with bootstrapping,
# and with the previously-converged t0 being passed.
bootstrapped_sim <- deconv(t, deaths, bs = TRUE, t0 = initial_sim[[3]])

# Accesses the t0 iteration history for both the initial
# and the bootstrapped simulations.
inft_t0 <- initial_sim[[2]]
inft_bs <- bootstrapped_sim[[2]]

# Gets the 2.5% and 97.5% quantile values for t0 across different
# iterations per day (with the bootstrapped data).
# min_inft_bs <- apply(inft_bs, 1, min)
# max_inft_bs <- apply(inft_bs, 1, max)

min_inft_bs <- apply(inft_bs, 1,
                     function(x) quantile(x, probs = c(.025, .975))["2.5%"])
max_inft_bs <- apply(inft_bs, 1,
                     function(x) quantile(x, probs = c(.025, .975))["97.5%"])

# TODO: remark from ryan: can we just do this once and extract both values?

# Accesses the simulated days of death to capture the model's fitness.
sim_deaths_tabulated <- initial_sim[[4]]

# Tabulates the actual days of death data to compare with the simulated result.
actual_deaths_tabulated <- tabulate(death_days, nbins = 310)

## Generate the final plot with the amount of days since the start of the year
# as the x-axis, and the number of people as the y-axis, that includes:
# - estimated fatal incidence trajectory (black line),
plot(1:310, inft_t0[, 100], col = "black", type = "l", lwd = 1,
     xlab = "Amount of days since 1st January 2020", ylab = "Num. of people",
     main = paste("Final COVID-19 Fatal Infections and Deaths Simulation"),
     ylim = c(0, 1800))

# - 95% uncertainty for the incidence trajectory (shaded),
lines(1:310, min_inft_bs, col = "gray", type = "l", lwd = 1, lty = "dashed")
lines(1:310, max_inft_bs, col = "gray", type = "l", lwd = 1, lty = "dashed")

# - estimated deaths (red dashed line),
lines(1:310, sim_deaths_tabulated, col = "red",
      type = "l", lwd = 1, lty = "dotdash")

# - actual deaths (blue line),
lines(1:310, actual_deaths_tabulated, col = "blue",
      type = "l", lwd = 1, lty = "solid")

# - a line indicating the first day of UK lockdown (red line).
abline(v = 84, col = "red", lwd = 2)

# Generates the legend for the final plot.
legend(x = "topright", inset = 0.05,
       legend = c("Est. Incidence Trajectory", "Est. Deaths", "Actual Deaths",
                  "UK Lockdown Start (Day 84)", "95% Uncertainty Area"),
       col = c("black", "red", "blue", "red", "gray"),
       lty = c("solid", "dotdash", "solid", "solid", "dotted"),
       lwd = c(1, 1, 1, 2, 1), cex = 0.8)

# TODO: what does this thing below do? @Joseph @Frans
x <- 1:length(min_inft_bs)
polygon(c(x, rev(x)), c(min_inft_bs, rev(max_inft_bs)), col = 'gray', density = 50, border = NA)
