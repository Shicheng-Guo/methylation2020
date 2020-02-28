function (x, k, method = NULL, m = 1, weights = 1, control = list()) 
{
  mc <- match.call()
  if (!all(row_norms(x) > 0)) 
    stop("Zero rows are not allowed.")
  args <- list(x = x, k = k, control = control)
  if (!missing(m) && !identical(m, 1)) 
    args <- c(args, list(m = m))
  if (!missing(weights) && !all(weights == 1)) 
    args <- c(args, list(weights = rep(weights, length.out = nrow(x))))
  skmeans_methods <- c(genetic = "genetic", pclust = "pclust", 
                       CLUTO = "CLUTO", gmeans = "gmeans", kmndirs = "kmndirs", 
                       LIH = "local_improvement_heuristic", LIHC = "local_improvement_heuristic_with_chains")
  if (!is.function(method)) {
    method <- if (is.null(method)) {
      if (is.null(args$m)) 
        "genetic"
      else "pclust"
    }
    else if (is.character(method)) {
      pos <- pmatch(tolower(method), tolower(names(skmeans_methods)))
      if (is.na(pos)) 
        stop(gettextf("Invalid skmeans method '%s'.", 
                      method))
      method <- skmeans_methods[pos]
    }
    else {
      stop("Invalid skmeans method.")
    }
    method <- get(sprintf(".skmeans_%s", method))
  }
  na <- names(args)
  nf <- names(formals(method))
  if (any(ind <- is.na(match(na, nf)))) 
    stop(gettextf("Given skmeans method lacks formals %s", 
                  paste(sQuote(na[ind]), collapse = " and ")))
  if (("m" %in% nf) && !("m" %in% na)) 
    args <- c(args, list(m = m))
  if (("weights" %in% nf) && !("weights" %in% na)) 
    args <- c(args, list(weights = rep(weights, length.out = nrow(x))))
  y <- do.call(method, args)
  y$call <- mc
  y
}
  