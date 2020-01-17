# To-do list for package:

# - clean up NAMESPACE and DESCRIPTION files; some info is unique to each file, and Tamara ran into issues using the package b/c she was able to install it without also installing minpack.lm and zoo packages. See: http://r-pkgs.had.co.nz/namespace.html#imports 

# - make it possible to truncate temperature ranges in predict.nbcurve when applying it iteratively to multiple TPCs.

# - provide more informative error messages and method selection for short time series (<=3 time points?)

# - allow for weighted regression based on reported standard error values

# - check out how warnings() are handled in mle2; can I bundle together error messages from iterated get.decurve.tpc() in a similar fashion?

# - consider wrapping suppressWarnings() around the uniroot call in get.decurve.tpc and other locations, to avoid redundant reporting of "Error in uniroot(f = objective, interval = c(-400, cf$topt)) : f() values at end points not of opposite sign"

#   - maybe target specific error message?

#   - would also be nice to turn off 'convergence failure: cod=1 (iteration limit 'maxit' reached) in grid.mle2...

