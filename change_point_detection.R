# Import packages
library(tidyverse)
library(lubridate)
library(changepoint)
library(strucchange)
library(ggplot2)
library(dplyr)
library(ocp)

# Load the prepared NOAA dataset
df <- read_csv("./data/noaa_water_levels_dataset.csv")

# Quick sanity check to confirm the signal looks reasonable.
plot(df$anomaly_smooth, type = "l")

# Choose the monitoring signal and matching time vector.
signal <- df$anomaly_smooth
time   <- df$datetime

# Remove missing values and keep time aligned with signal.
keep <- which(!is.na(signal))
signal <- signal[keep]
time   <- time[keep]

##########################################
### Frequentist Change Point Detection ###

cp_fit <- cpt.mean(signal,
                   method = "PELT",
                   penalty = "BIC")
cps <- cpts(cp_fit)

# Plot the signal with detected change points.
plot(time, signal, type = "l", lwd = 1,
     xlab = "Time", ylab = "Smoothed anomaly (m)")

# Add vertical dashed lines at change points
abline(v = time[cps], col = "red", lty = 2, lwd = 2)

# Add segment means as horizontal lines
breaks <- c(1, cps + 1, length(signal) + 1)
for (i in 1:(length(breaks)-1)) {
  s <- breaks[i]
  e <- breaks[i+1] - 1
  mu <- mean(signal[s:e], na.rm = TRUE)
  segments(x0 = time[s], y0 = mu, x1 = time[e], y1 = mu,
           col = "red", lwd = 2)
}

#######################################
### Bayesian Change Point Detection ###

# Convert to numeric and re drop any missing values (keep alignment).
y <- as.numeric(signal)
t <- time
keep2 <- which(!is.na(y) & !is.na(t))
y <- y[keep2]
t <- t[keep2]

# Standardize so the Bayesian model works on a stable scale.
y <- as.numeric(scale(y))

# Define constant hazard with expected regime length ~600
haz <- function(x, lambda) const_hazard(x, lambda = 600)

# Run Bayesian online CPD.
# - cpthreshold: how confident the model must be to flag a change point
# - minsep: prevents change points from being too close together
# - minRlength: prevents very short regimes from being treated as meaningful
ocp_result <- onlineCPD(
  datapts = y,
  getR = TRUE,
  missPts = "none",
  hazard_func = haz,
  cpthreshold = 0.85,
  minsep = 100,
  minRlength = 40,
  printupdates = FALSE
)

summary(ocp_result)

plot(
  ocp_result,
  xlab = "Time Index (6-minute intervals)",
  ylab = "Standardized Smoothed Anomaly"
)