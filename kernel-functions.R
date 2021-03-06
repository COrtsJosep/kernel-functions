################ PART ONE: DEFINITIONS ################
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
          weights[i]<-(3/4)*(1-u^2) 
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
            weights[i]<-(3/4)*(1-u^2) 
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





################ PART TWO: REPRODUCTION OF THE PLOTS ################
### FIRST PLOT: COMPARISON AMONG KERNEL WEIGHTS ###
NW_weights <- function(ref, x, y, width, kernel){
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
          weights[i]<-(3/4)*(1-u^2) 
        } else if (kernel == "tri"){ 
          weights[i]<-(70/81)*((1-abs(u)^3)^3) 
        } else if (kernel != "gau"){
          return("Please select and adequate kernel: 'ep' for Epanechnikov, 'tri' for Tricube or 'gau' for Gaussian")
        }
      }
    }
  }
  return(weights/sum(weights))
}
set.seed(11)
n <- 120
x <- seq(-2.7,2.7, length = n)
y <- 1/(1+(exp(-x)))
plot(seq(min(x), max(x), length = length(x)), NW_weights(0,x,y, 1, kernel = "ep"), col = "black", main = "Comparison Among Kernel Weights", xlab = "X", ylab = "Weights", ylim = c(0,0.045))
points(seq(min(x), max(x), length = length(x)), NW_weights(0,x,y, 1, kernel = "tri"), col = "red")
points(seq(min(x), max(x), length = length(x)), NW_weights(0,x,y, 1, kernel = "gau"), col = "blue")
legend(1.3,0.0435, legend=c("Epanechnikov", "Tricube", "Gaussian"), col=c("black", "red", "blue"), pch = 1, cex=0.8)

### SECOND PLOT: COMPARISON OF THE PERFORMANCES AMONG KERNELS ###
set.seed(11)
n <- 40
x <- seq(-5,5, length = n)
y <- 1/(1+(exp(-x))) + rnorm(n, 0, 0.1)
plot(x,y, main = "Comparison of the Performances Among Kernels",xlab = "X", ylab = "Y")
lines(x, 1/(1+exp(-x)),lwd = 2)
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y, LOOCV(x,y, "ep", "NW", 100)$optimal, "ep", "NW"), lty = 2)
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y, LOOCV(x,y, "tri", "NW", 100)$optimal, "tri", "NW"), lty = 2, col = "red")
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y, LOOCV(x,y, "gau", "NW", 100)$optimal, "gau", "NW"), lty = 2, col = "blue")
legend(2,0.36, legend=c("True Curve", "Epanechnikov", "Tricube", "Gaussian"), col=c("black","black", "red", "blue"), lty = c(1,2,2,2), lwd = c(2,1,1,1), cex = 0.8 )

### THIRD PLOT: COMPARISON BETWEEN NW AND LOCAL LINEAR REGRESSION ###
set.seed(11)
n <- 60
x <- seq(-5,5, length = n*2/3)
x <- sort(c(x, rnorm(n/3, 0, 1.7)))
y <- cos(x)+rnorm(n,0, 0.3)
plot(x,y, main = "Comparison Between NW and Local Linear Regression",xlab = "X", ylab = "Y")
lines(x, cos(x))
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"gau", "NW", 100)$optimal  , "gau", "NW"), col = "red")
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"gau", "lin", 100)$optimal  , "gau", "lin"), col = "blue")
legend(1.8,1.2, legend=c("True Curve", "Nadaraya-Watson", "Local Linear Reg."), col=c("black", "red", "blue"), lty=1, cex=0.8)

### FORTH PLOT: COMPARISON BETWEEN LOCAL LINEAR AND QUADRATIC REGRESSIONS ###
set.seed(8)
n <- 50
x <- seq(-4,4, length = n)
y <- x+2*x^2+rnorm(n,0, 4)
plot(x,y, ylim = c(-10,35), main = "Comparison Between Local Linear and Quadratic Regressions",xlab = "X", ylab = "Y")
lines(x, x+2*x^2)
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"tri", "lin", 100)$optimal  , "tri", "lin"), col = "blue")
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"tri", "q", 100)$optimal  , "tri", "q"), col = "red")
legend(1.45,2, legend=c("True Curve", "Local Linear Reg.", "Local Quadratic Reg."), col=c("black", "blue", "red"), lty=1, cex=0.8)

### FIFTH PLOT: ERROR CURVES ###
set.seed(10)
n <- 50
x <- seq(0,6, length = n)
x <- c(x, rnorm(2*n/3, 3, 1))
y <- sin(x)+rnorm(n+2*n/3, 0, 0.4)

result <- LOOCV(x,y,"tri", "lin", 500,lambda_max =  2)
plot(result$set_lambdas, result$errors, type = "l", ylim = c(0.12, 0.245), xlab = "Widths", ylab = "Associated Error", main = "Error Curves")
points( result$optimal, min(result$errors), pch = 4, col = "blue" , cex = 2, lwd = 2)
lines(x = c(result$optimal,result$optimal), y = c(0.12, 0.125), col = "blue", lwd = 2.5)

newcurve <- function(n){
  x <- seq(0,6, length = n)
  x <- c(x, rnorm(2*n/3, 3, 1))
  y <- sin(x)+rnorm(n+2*n/3, 0, 0.4)
  result <- LOOCV(x,y,"tri", "lin", 500,lambda_max =  2)
  lines(result$set_lambdas, result$errors)
  points( result$optimal, min(result$errors), pch = 4, col = "blue" , cex = 2, lwd = 2)
  lines(x = c(result$optimal,result$optimal), y = c(0.12, 0.125), col = "blue", lwd = 2.5)
}
replicate(4, newcurve(60))





################ PART THREE: EXAMPLES OF USE ################
n <- 60
x <- seq(-5,5, length = n)
y <- sin(x)+rnorm(n, 0, 0.4)
plot(x,y)
lines(x, sin(x))
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"ep", "NW", 100)$optimal  , "ep", "NW"))
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"ep", "lin", 100)$optimal  , "ep", "lin"), col = "blue")
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"ep", "q", 100)$optimal  , "ep", "q"), col = "red")

data(iris)
x <- iris$Sepal.Length
y <- iris$Petal.Length
plot(x,y)
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"gau", "lin", 100)$optimal  , "gau", "lin"), col = "black")
CIE <- CI(x,y,LOOCV(x,y,"gau", "lin", 100)$optimal, "gau", "lin", interval = 0.01, times = 800)
lines( seq(min(x), max(x), length = length(x)), CIE$upper, col = "blue", lty = 2)
lines( seq(min(x), max(x), length = length(x)), CIE$lower, col = "blue", lty = 2)

set.seed(9)
n <- 50
x <- seq(-4,6, length = n*3/4)
x <- sort(c(x, rnorm(n/4, 2, 1.5)))
y <- -x+2*x^2-1.5*x^3+rnorm(n, 0, 20)
plot(x,y)
lines(x, -x+2*x^2-1.5*x^3)
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"tri", "NW", 100)$optimal  , "tri", "NW"), col = "red")
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"tri", "lin", 100)$optimal  , "tri", "lin"), col = "blue")
lines(seq(min(x), max(x), length = length(x)), kernel_reg(x,y,LOOCV(x,y,"tri", "q", 100)$optimal  , "tri", "q"), col = "green")
