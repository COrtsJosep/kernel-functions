### EVALUATION AT POINT ###
kernel_point <- function(ref, x, y, width, kernel,NW_lin_or_q){
  weights <- rep(0,length(y))
  for (i in 1:length(x)){
    u <- abs( (ref-x[i])/width )
    if (kernel == "gau"){
      weights[i] <- (( 1/sqrt(2*pi) ) * exp( -(1/2)*u^2 ))
    }
    else {
      if (u>1){ weights[i]<-0} 
      else{ 
        if (kernel == "ep"){ 
          weights[i]<-(3/4)*((1-u)^2) 
        } else if (kernel == "tri"){ 
          weights[i]<-(70/81)*((1-abs(u)^3)^3) 
        } else if (kernel != "gau"){
          return("Please select and adequate kernel: 'ep' for Epanechnikov, 'tri' for Tricube or 'gau' for Gaussian")
        }
      }
    }
  }
  if (NW_lin_or_q=="NW"){ return(sum(weights*y)/sum(weights)) }
  else if (NW_lin_or_q=="lin"){
    W <- diag(weights)
    B <- matrix( c(rep(1, length(x)), x), length(x), 2)
    return(t(c(1, ref)) %*% solve( t(B) %*% W %*% B ) %*% t(B) %*% W %*% y)
  }
  else if (NW_lin_or_q=="q"){
    W <- diag(weights)
    B <- matrix( c(rep(1, length(x)), x, x^2), length(x), 3)
    return(t(c(1, ref, ref^2)) %*% solve( t(B) %*% W %*% B ) %*% t(B) %*% W %*% y)
  }
  else {return("Please select and adequate fitted model: 'NW' for a weighted average, 'lin' for a least squares linear model or 'q' for a least squares quadratic model")
  }
}

### REGRESSION ###
kernel_reg <- function(x, y, width, kernel, NW_lin_or_q){
  if (kernel!="ep"&kernel!="tri"&kernel!="gau"){
    return("Please select and adequate kernel: 'ep' for Epanechnikov, 'tri' for Tricube or 'gau' for Gaussian")
  }
  if (NW_lin_or_q!="NW"&NW_lin_or_q!="lin"&NW_lin_or_q!="q"){
    return("Please select and adequate fitted model: 'NW' for a weighted average, 'lin' for a least squares linear model or 'q' for a least squares quadratic model")
  }
  if (width =="optimal"){
    width <- LOOCV(x,y,kernel, NW_lin_or_q, precision = 200)$optimal
  }
  sapply(seq(min(x), max(x), length =length(x)), FUN = kernel_point, x, y, width, kernel, NW_lin_or_q)
}

### LOOCV ERRORS AND OPTIMISATION ###
LOOCV <- function(x, y, kernel, NW_lin_or_q, precision, lambda_min = max_distance(x)*scalar_min, lambda_max = 5*sigma/nthroot(n,5)){
  if (kernel!="ep"&kernel!="tri"&kernel!="gau"){
    return("Please select and adequate kernel: 'ep' for Epanechnikov, 'tri' for Tricube or 'gau' for Gaussian")
  }
  if (NW_lin_or_q!="NW"&NW_lin_or_q!="lin"&NW_lin_or_q!="q"){
    return("Please select and adequate fitted model: 'NW' for a weighted average, 'lin' for a least squares linear model or 'q' for a least squares quadratic model")
  }
  library(pracma)
  sigma <- sd(x)
  n <- length(x)
  max_distance <- function(x){
    distances <- c()
    sorted <- sort(x)
    for (i in 1:n-1){
      distances <- append(distances, abs(sorted[i]-sorted[i+1]))
    }
    return(max(distances))
  }
  if (NW_lin_or_q == "NW"){
    scalar_min <- 1.05
  } else if (NW_lin_or_q == "lin"){
    scalar_min <- 2.05
  } else { scalar_min <- 3.5 }
  if (kernel == "gau"){ scalar_min <- scalar_min/3}
  kernel_t <- function(ref, x, y, width, kernel,NW_lin_or_q, obs0){
    weights <- rep(0,n)
    for (i in 1:n){
      u <- abs( (ref-x[i])/width )
      if (kernel == "gau"){
        weights[i] <- (( 1/sqrt(2*pi) ) * exp( -(1/2)*u^2 ))
      }
      else {
        if (u>1){ weights[i]<-0} 
        else{ 
          if (kernel == "ep"){ 
            weights[i]<-(3/4)*((1-u)^2) 
          } else if (kernel == "tri"){ 
            weights[i]<-(70/81)*((1-abs(u)^3)^3) 
          } else if (kernel != "gau"){
            return("Please select and adequate kernel: 'ep' for Epanechnikov, 'tri' for Tricube or 'gau' for Gaussian")
          }
        }
      }
    }
    if (NW_lin_or_q=="NW"){ W <- weights/sum(weights) }
    else if (NW_lin_or_q=="lin"){
      W <- diag(weights)
      B <- matrix( c(rep(1, n), x), n, 2)
      W <- t(c(1, ref)) %*% solve( t(B) %*% W %*% B ) %*% t(B) %*% W
    }
    else {
      W <- diag(weights)
      B <- matrix( c(rep(1, n), x, x^2), n, 3)
      W <- t(c(1, ref, ref^2)) %*% solve( t(B) %*% W %*% B ) %*% t(B) %*% W
    }
    return( list( Skk = W[obs0] , predicted = W %*% y) )
  }
  CVE <- rep(0, precision)
  lambdas <- seq( lambda_min, lambda_max, length = precision )
  for (i in 1:precision ){
    lmb <- lambdas[i]
    MSEs <- rep(0, n)
    for (j in 1:n){
      MSEs[j] <- ( (y[j]-kernel_t(x[j], x, y, lmb, kernel, NW_lin_or_q, j)$predicted)/(1-kernel_t(x[j], x, y, lmb, kernel, NW_lin_or_q, j)$Skk) )^2
    }
    CVE[i] <- mean(MSEs)
  }
  return(list(optimal = lambdas[which( CVE==min(CVE) )], set_lambdas = lambdas, errors = CVE, rule_thumb = sigma/nthroot(n,5)))
}

### BOOTSTRAP CI USING K-FOLD CV ###
CI <- function(x, y, width, kernel, NW_lin_or_q, interval, times){
  n <- length(x)
  max <- max(x)
  min <- min(x)
  
  leave_out_reg <- function(x, y, width, kernel, NW_lin_or_q){
    indices <- sample(1:n, n*0.8)
    if (kernel!="ep"&kernel!="tri"&kernel!="gau"){
      return("Please select and adequate kernel: 'ep' for Epanechnikov, 'tri' for Tricube or 'gau' for Gaussian")
    }
    if (NW_lin_or_q!="NW"&NW_lin_or_q!="lin"&NW_lin_or_q!="q"){
      return("Please select and adequate fitted model: 'NW' for a weighted average, 'lin' for a least squares linear model or 'q' for a least squares quadratic model")
    }
    if (width =="optimal"){
      width <- LOOCV(x,y,kernel, NW_lin_or_q, precision = 200)$optimal
    }
    sapply(seq(min, max, length = n), FUN = kernel_point, x[indices], y[indices], width, kernel, NW_lin_or_q)
  }
  
  results <- replicate(times, leave_out_reg(x, y, width, kernel, NW_lin_or_q))
  
  intervals <- c()
  for (i in 1:n){
    intervals <- rbind(intervals, as.vector(quantile(results[i,], probs = c(interval/2, 1-interval/2))) )
  }
  return(list(upper = intervals[,2], lower = intervals[,1]))
}
