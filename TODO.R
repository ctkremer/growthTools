# To-do list for package:

# - dependencies: try using @importFrom rather than @import to specify specific functions to import from other packages. This might help avoid problems like bbmle::slice and dplyr::slice, where function names are redundance across packages, but potentially with quite different functionality.

# - clean up example data sets - they contain more columns than are actually used/useful in examples and vignettes. also, consider renaming example_TPC_data

# - make it possible to truncate temperature ranges in predict.nbcurve when applying it iteratively to multiple TPCs.

# - provide more informative error messages and method selection for short time series (<=3 time points?)

# - allow for weighted regression based on reported standard error values

# - check out how warnings() are handled in mle2; can I bundle together error messages from iterated get.decurve.tpc() in a similar fashion?

# - consider wrapping suppressWarnings() around the uniroot call in get.decurve.tpc and other locations, to avoid redundant reporting of "Error in uniroot(f = objective, interval = c(-400, cf$topt)) : f() values at end points not of opposite sign"

#   - maybe target specific error message?

#   - would also be nice to turn off 'convergence failure: cod=1 (iteration limit 'maxit' reached) in grid.mle2...

# - create unit testing for the package to ensure stability of calculations...

