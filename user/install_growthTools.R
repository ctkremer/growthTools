
# Need to install devtools, if not already there
library(devtools)

# remove any old versions of the package before installing
remove.packages('growthTools')

# Install the growthTools package!
devtools::install_github("ctkremer/growthTools",auth_token = '6af60cb683ef4cbab1d6b25e0ec2ed6b925c831b',build_vignettes = TRUE)

# Load package!
library(growthTools)

# To view a document outlining how this package works and what it contains, try:
vignette("growthTools_vignette",package="growthTools")