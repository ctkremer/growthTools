

# Need to install devtools, if not already there
library(devtools)

#### growthTools depends on mleTools, which you can install via:

devtools::install_github("ctkremer/mleTools",upgrade='ask')


#### Install development version of growthTools:

# remove any old versions of the package before installing
# remove.packages('growthTools')

# Install without building vignette:
devtools::install_github("ctkremer/growthTools@develop",build_vignettes = F,upgrade='ask')

# Note: if you want to build the vignette, change build_vignettes to T, may take ~10 min to build. Alternatively, view pre-built version.

# Load growthTools
library(growthTools)

# IF you built the vignette, view it via:
vignette("growthTools_vignette",package="growthTools")

