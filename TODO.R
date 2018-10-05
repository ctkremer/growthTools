# To-do list for package:

# - allow for weighted regression based on reported standard error values
# - combine nbcurve and nbcurve2 into a single function that depends on either copt or opt (see alternative parameterizations of distributions, e.g.)
# - check out how warnings() are handled in mle2; can I bundle together error messages from iterated get.decurve.tpc() in a similar fashion?
# - consider wrapping suppressWarnings() around the uniroot call in get.decurve.tpc and other locations, to avoid redundant reporting of "Error in uniroot(f = objective, interval = c(-400, cf$topt)) : f() values at end points not of opposite sign"
#   - maybe target specific error message?
#   - would also be nice to turn off 'convergence failure: cod=1 (iteration limit 'maxit' reached) in grid.mle2...



# NOTES:

#library(devtools)

# at setup:
#devtools::use_readme_rmd()
#devtools::use_build_ignore("NEWS.md")
#devtools::use_build_ignore("TODO.R")

# during active development
#devtools::load_all()  # re-load all files in /R

# devtools::build_vignettes()
