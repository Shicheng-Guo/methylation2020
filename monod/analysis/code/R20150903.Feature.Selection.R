
run.name <-"feature-bc"
install.packages("caret")
install.packages("randomForest")
install.packages("ipred")
install.packages("gbm")

## Load early to get the warnings out of the way:
library("caret")
library("randomForest")
library("ipred")
library("gbm")

set.seed(1)

## Set up artificial test data for our analysis
n.var <- 20
n.obs <- 200
x <- data.frame(V = matrix(rnorm(n.var*n.obs), n.obs, n.var))
n.dep <- floor(n.var/5)
cat( "Number of dependent variables is", n.dep, "\n")
m <- diag(n.dep:1)

## These are our four test targets
y.1 <- factor( ifelse( x$V.1 >= 0, 'A', 'B' ) )
y.2 <- ifelse( rowSums(as.matrix(x[, 1:n.dep]) %*% m) >= 0, "A", "B" )
y.2 <- factor(y.2)
y.3 <- factor(rowSums(x[, 1:n.dep] >= 0))
y.4 <- factor(rowSums(x[, 1:n.dep] >= 0) %% 2)
  
control <- rfeControl(functions = rfFuncs, method = "boot", verbose = FALSE,returnResamp = "final", number = 50)

if ( require("multicore", quietly = TRUE, warn.conflicts = FALSE) ) {
  control$workers <- multicore:::detectCores()
  control$computeFunction <- mclapply
  control$computeArgs <- list(mc.preschedule = FALSE, mc.set.seed = FALSE)
}


sizes <- 1:6
## Use randomForest for prediction
profile.1 <- rfe(x, y.1, sizes = sizes, rfeControl = control)
cat( "rf : Profile 1 predictors:", predictors(profile.1), fill = TRUE )
profile.2 <- rfe(x, y.2, sizes = sizes, rfeControl = control)
cat( "rf : Profile 2 predictors:", predictors(profile.2), fill = TRUE )
profile.3 <- rfe(x, y.3, sizes = sizes, rfeControl = control)
cat( "rf : Profile 3 predictors:", predictors(profile.3), fill = TRUE )
profile.4 <- rfe(x, y.4, sizes = sizes, rfeControl = control)
cat( "rf : Profile 4 predictors:", predictors(profile.4), fill = TRUE )


