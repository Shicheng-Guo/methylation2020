regessinf <- function(x, pooled, sortvar, level.comb=x$level.comb){

# paremter passing 
  if (!inherits(x, "meta"))
    stop("Argument 'x' must be an object of class \"meta\"")  
  if (missing(pooled)){
    if (length(x$comb.fixed)==0 & length(x$comb.random)==0)
      pooled <- "fixed"
    if (length(x$comb.fixed)>0 & length(x$comb.random)==0)
      if (x$comb.fixed) pooled <- "fixed"
    else pooled <- "NoNe"
    if (length(x$comb.fixed)==0 & length(x$comb.random)>0)
      if (x$comb.random) pooled <- "random"
    else pooled <- "NoNe"
    if (length(x$comb.fixed)>0 & length(x$comb.random)>0){
      if (x$comb.fixed)
        pooled <- "fixed"
      if (!x$comb.fixed & x$comb.random)
        pooled <- "random"
      if (!x$comb.fixed & !x$comb.random)
        pooled <- "NoNe"
    }
  }
  ##
  if (pooled=="NoNe")
    stop("Parameters \"comb.fixed\" and \"comb.random\" in object '",
         deparse(substitute(x)),
         "' are either 'FALSE' or 'NULL'. ",
         "Please use argument \"pooled=fixed\" or \"pooled=random\" ",
         "to select meta-analytical model.")
  imeth <- charmatch(tolower(pooled), c("fixed", "random"), nomatch = NA)
  ##
  if (is.na(imeth)) 
    stop("'pooled' should be \"fixed\" or \"random\"")
  ##
  pooled <- c("fixed", "random")[imeth]
  if (length(level.comb)==0){
    warning("level.comb set to 0.95")
    level.comb <- 0.95
  }
  
  
# body
  
k.all <- length(x$TE)
for (i in 1:k.all){
  sel <- -i
  ##
  n.e <- x$n.e
  n.c <- x$n.c
  n   <- x$n
  ##
  event.e <- x$event.e
  event.c <- x$event.c
  event   <- x$event
  
  if (inherits(x, "metabin"))
    m <- metabin(event.e[sel], n.e[sel], event.c[sel], n.c[sel],
                 method=x$method, sm=x$sm,
                 incr=x$incr, allincr=x$allincr, addincr=x$addincr,
                 allstudies=x$allstudies, MH.exact=x$MH.exact,
                 RR.cochrane=x$RR.cochrane,
                 level=level.comb, level.comb=level.comb,
                 hakn=x$hakn,
                 method.tau=x$method.tau,
                 tau.preset=x$tau.preset, TE.tau=x$TE.tau,
                 warn=FALSE)
  ##
}




}