# This file is for testing changepoint functions with toy examples

# Call functions script
source("./changePoint/changePointFunctions.R")

# One change point test
x <- 1:15
y1 <- 2 + 3*x[1:8] + rnorm(8, mean = 0, sd = 2)
y2 <- 60 + -5*x[9:15] + rnorm(7, mean = 0, sd = 4)
y <- c(y1,y2)
plot(y)



results <- segment(x, y, 1)
lines(results$fittedValues)
# Two Change Point test
x <- 1:20
y1 <- 2 + 3*x[1:8] + rnorm(8, mean = 0, sd = 2)
y2 <- 60 + -5*x[9:15] + rnorm(7, mean = 0, sd = 8)
y3 <- -120 + 8*x[16:20] + rnorm(5, mean = 0, sd = 5)
y <- c(y1,y2,y3)
plot(y)

solvedParams <- twoChangePoints(x, y, 8, 15)
fittedValues <- fitTwoChangePoints(solvedParams, x, 8, 15)
lines(fittedValues)

















