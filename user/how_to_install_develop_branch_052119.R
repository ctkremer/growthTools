

# Need to install devtools, if not already there
library(devtools)

#issues to build vignettes, need a development version of pillar package
devtools::install_github("r-lib/pillar")

#### Install development version of growthTools:

# remove any old versions of the package before installing
# remove.packages('growthTools')

devtools::install_github("ctkremer/growthTools@develop",auth_token = '6af60cb683ef4cbab1d6b25e0ec2ed6b925c831b',build_vignettes = T,upgrade_dependencies=F,build_opts = c("--no-resave-data", "--no-manual"))

library(growthTools)

vignette("growthTools_vignette",package="growthTools")



