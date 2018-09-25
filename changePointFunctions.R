### File for functions used in creating changepoint estimation


# This function calculates the parameters of a segmented linear relationship with one changepoint.
# Returns solvedParams with the interecept, slope, and change in slope
oneChangePoint <- function(xVector, yVector, givenChangePoint){
  
  # Creating moment data needed to compute parameter estimates
  xBeforeChangePoint <- xVector[xVector <= givenChangePoint]
  xAfterChangePoint <- xVector[xVector > givenChangePoint]
  xBeforeChangePoints <- 1:length(xBeforeChangePoint)
  xAfterChangePoint <- 1:length(xAfterChangePoint)
  
  yBeforeChangePoint <- yVector[1:length(xBeforeChangePoint)]
  yAfterChangePoint <- yVector[(length(xBeforeChangePoint) + 1) : length(xVector)]
  
  avgXBefore <- mean(xBeforeChangePoint)
  avgXAfter <- mean(xAfterChangePoint)
  avgYBefore <- mean(yBeforeChangePoint)
  avgYAfter <- mean(yAfterChangePoint)
  
  
  avgXbeforeSq <- mean(xBeforeChangePoint ^ 2)
  avgXafterSq <- mean(xAfterChangePoint ^ 2)
  
  xyBefore <- mean(xBeforeChangePoint * yBeforeChangePoint)
  xyAfter <- mean(xAfterChangePoint * yAfterChangePoint)
  
  cp <- givenChangePoint
  
  # creating right and left side of of system of equations to solve for parameters
  lhsMatrix <- matrix(c(
            2, (avgXBefore + cp + avgXAfter), avgXAfter,
            (avgXBefore + cp + avgXAfter), (avgXbeforeSq + cp^2 + 2*cp*avgXAfter + avgXafterSq), (cp*avgXAfter + avgXafterSq),
            avgXAfter, (cp*avgXAfter + avgXafterSq), avgXafterSq
            
  ), 3, 3)
  
  rhs <- c( (avgYAfter + avgYBefore), 
            (xyBefore + cp*avgYAfter + xyAfter),
            (xyAfter)
    
  )
  

  solvedParams <- solve(lhsMatrix, rhs)
  return(solvedParams)
  
  
}

# This function is used to calculate segmented linear relationships with two changepoints
# Returns solvedParams with the intercept, slope, first change in slope, and second change in slope
twoChangePoints <- function(xVector, yVector, c1, c2){
  
  # Creating Moment data needed to create parameter estimates
  xBeforeChangePoint <- xVector[xVector <= c1]
  xMidChangePoint <- xVector[xVector > c1 & xVector <= c2]
  xAfterChangePoint <- xVector[xVector > c2]
  
  xBeforeChangePoints <- 1:length(xBeforeChangePoint)
  xMidChangePoint <- 1:length(xMidChangePoint)  
  xAfterChangePoint <- 1:length(xAfterChangePoint)
  
  
  yBeforeChangePoint <- yVector[1:length(xBeforeChangePoint)]
  yMidChangePoint <- yVector[(length(xBeforeChangePoint) + 1) : (length(xMidChangePoint) + length(xBeforeChangePoint))]
  yAfterChangePoint <- yVector[((length(xMidChangePoint) + length(xBeforeChangePoint)) + 1) : (length(yMidChangePoint) + length(xAfterChangePoint) + length(xBeforeChangePoint))]
  
  
  
  avgXBefore <- mean(xBeforeChangePoint)
  avgXMid <- mean(xMidChangePoint)
  avgXAfter <- mean(xAfterChangePoint)
  
  avgYBefore <- mean(yBeforeChangePoint)
  avgYAfter <- mean(yAfterChangePoint)
  avgYMid <- mean(yMidChangePoint)
  
  
  avgXbeforeSq <- mean(xBeforeChangePoint ^ 2)
  avgXafterSq <- mean(xAfterChangePoint ^ 2)
  avgXmidSq <- mean(xMidChangePoint ^ 2)
  
  xyBefore <- mean(xBeforeChangePoint * yBeforeChangePoint)
  xyAfter <- mean(xAfterChangePoint * yAfterChangePoint)
  xyMid <- mean(xMidChangePoint * yMidChangePoint)
  
  # computing the difference between changepoints as the value of second change point
  c2 <- c2 - c1
  
  # Creating the right and left hand side of the system of equations
  lhsMatrix <- matrix(c(
    3, (avgXBefore + c1 + avgXMid + c1 + c2 + avgXAfter), (avgXMid + c2 + avgXAfter), (avgXAfter),
    (avgXBefore + c1 + avgXMid + c1 + c2 + avgXAfter), (avgXbeforeSq + c1^2 + 2*c1*avgXMid + avgXmidSq + c1^2 + 2*c1*c2 + 2*c1*avgXAfter + c2^2 + 2*avgXAfter*c2 + avgXafterSq), (c1*avgXMid + avgXmidSq + c1*c2 + c1*avgXAfter + c2^2 + 2*c2*avgXAfter + avgXafterSq), (c1*avgXAfter + c2*avgXAfter + avgXafterSq),
    (avgXMid + c2 + avgXAfter), (c1*avgXMid + avgXmidSq + c1*c2 + c2^2 + 2*c2*avgXAfter + c1*avgXAfter + avgXafterSq), (avgXmidSq + c2^2 + 2*c2*avgXAfter + avgXafterSq), (c2*avgXAfter + avgXafterSq),
    (avgXAfter), (c1*avgXAfter + c2*avgXAfter + avgXafterSq), (c2*avgXAfter + avgXafterSq), (avgXafterSq)
    
  ), 4, 4)
  
  rhs <- c( (avgYMid + avgYBefore + avgYAfter), 
            (xyBefore + c1*avgYMid + xyMid + c1*avgYAfter + c2*avgYAfter + xyAfter),
            (xyMid + xyAfter + c2*avgYAfter),
            (xyAfter)
            
  )
  
  # solving system and returning parameters
  solvedParams <- solve(lhsMatrix, rhs)
  return(solvedParams)
  
  
}

# This function fits the parameters of a one change point model to the data, and returns a vector of the fitted values
fitOneChangePoint <- function(paramsVector, x, c1){
  
  
  c1 <- length(x[x<=c1])
  
  
  firstHalf <- paramsVector[1] + paramsVector[2]*x[1:(c1)]
  secondHalf <- firstHalf[c1] + (paramsVector[2] + paramsVector[3]) * x[1:(length(x) - c1)]
  fittedValues <- c(firstHalf, secondHalf)
  
  return(fittedValues)
  
}

# This function fits the parameters of a two change point model to the data, and returns a vector the fitted values
fitTwoChangePoints <- function(paramsVector, x, c1, c2){
  
  c1 <- length(x[x <= c1])
  c2 <- length(x[x <= c2])
  first <- paramsVector[1] + paramsVector[2]*x[1:(c1)]
  second <- first[c1] + (paramsVector[2] + paramsVector[3]) * x[1:(c2-c1)]
  third <- second[(c2 - c1)] + (paramsVector[2] + paramsVector[3] + paramsVector[4]) * x[1:(length(x) - length(first) - length(second))]
  fittedValues <- c(first, second, third)
  
  return(fittedValues)
  
}

# This function returns the residual sum of squares
rssFunc <- function(fittedY, actualY){
  
  rss <- sum((actualY - fittedY)^2)
  return(rss)
  
}


# This function finds the most optimal changepoint for either 1 or 2 changepoints

segment <- function(x, y, numChangePoint){
  
  if (numChangePoint == 1){
    rss <- numeric(length((2:(length(x) - 1))))
    count <- 1
    for (c in (2:(length(x)-1))){
      
      solvedParams <- oneChangePoint(x, y, c)
      fittedValues <- fitOneChangePoint(solvedParams, x, c)
      rss[count] <- rssFunc(fittedValues, y)
      count <- count + 1
      
    }
    
    changePoint <- which(rss == min(rss)) + 1 
    params <- oneChangePoint(x, y, changePoint)
    params <- c(params, changePoint)
    newValues <- oneCpOptimizer(params, x, y)
    params <- newValues$par
    
    return(params)
  }
  
  if (numChangePoint == 2){
    
    pairs <- combn((2:(length(x)-1)), m = 2) 
    rss <- numeric(length(pairs[1,]))
    
    for (i in 1:length(pairs[1,])){
    
      params <- twoChangePoints(x, y, pairs[1,i], pairs[2,i])
      fittedValues <- fitTwoChangePoints(params, x, pairs[1,i], pairs[2,i])
      rss[i] <- rssFunc(fittedValues, y)
    }
    
    indexMin <- which(rss == min(rss))
    c1 <- pairs[1,indexMin]
    c2 <- pairs[2,indexMin]
    
    solvedParams <- twoChangePoints(x, y, c1, c2)
    params <- c(solvedParams, c1, c2)
    newValues <- twoCpOptimizer(params, x, y)
    params <- newValues$par
    return(params)
    
  }
}


# Trying to use the optim function to make some tweaks to our found change point

oneCpOptimizer <- function(params, x, y){
  
  # alpha is intercept, beta is initial slope, betaDelta is second slope, c is changepoint
  functionToOptimize <- function(params){
    alpha <- params[1]
    beta <- params[2]
    betaDelta <- params[3]
    c <- params[4]
    
    x1 <- x[x <= c] 
    x1 <- 1:length(x1)
    x2 <- x[x > c]
    x2 <- 1:length(x2)
    y1 <- y[1:length(x1)]
    y2 <- y[(length(x1) + 1):length(y)]
    
    # function we are optimizing - sum of squares
    (sum((y1 - alpha - beta*x1)^2))*(1/length(x1)) + sum((y2 - alpha - beta*c - (beta + betaDelta)*x2)^2)*(1/length(x2))
    
    
  }
  
  newParams <- optim(params, functionToOptimize, method = "BFGS")
  return(newParams)
  
}

twoCpOptimizer <- function(params, x, y){
  
  # alpha is intercept, beta is initial slope, betaDelta is second slope, c is changepoint
  functionToOptimize <- function(params){
    alpha <- params[1]
    beta <- params[2]
    betaDeltaOne <- params[3]
    betaDeltaTwo <- params[4]
    c1 <- params[5]
    c2 <- params[6]
    
    x1 <- x[x <= c1] 
    x1 <- 1:length(x1)
    x2 <- x[x > c1 & x <= c2]
    x2 <- 1:length(x2)
    x3 <- x[x > c2]
    x3 <- 1:length(x3)
    y1 <- y[1:length(x1)]
    y2 <- y[(length(x1) + 1):(length(x1) + length(x2))]
    y3 <- y[(length(x1)+length(x2) + 1) : length(y)]
    
    # function we are optimizing - sum of squares
    (sum((y1 - alpha - beta*x1)^2))*(1/length(x1)) + sum((y2 - alpha - beta*c1 - (beta + betaDeltaOne)*x2)^2)*(1/length(x2)) + (sum((y3 - alpha - beta*c1 - (beta + betaDeltaOne)*(c2 - c1) - (beta + betaDeltaOne + betaDeltaTwo)*x3)^2)) * (1/length(x3))
    
    
  }
  
  newParams <- optim(params, functionToOptimize, method = "BFGS")
  return(newParams)
  
}









