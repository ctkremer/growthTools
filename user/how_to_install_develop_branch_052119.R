

# Need to install devtools, if not already there
library(devtools)

#issues to build vignettes, need a development version of pillar package
devtools::install_github("r-lib/pillar")

#### Install development version of growthTools:

# remove any old versions of the package before installing
# remove.packages('growthTools')

devtools::install_github("ctkremer/growthTools@develop",auth_token = 'f8fe6d72858a84490060478a38802cc873b172e9',build_vignettes = T,upgrade_dependencies=F,build_opts = c("--no-resave-data", "--no-manual"))

library(growthTools)

vignette("growthTools_vignette",package="growthTools")




devtools::install_github("ctkremer/growthTools@master",auth_token = 'f8fe6d72858a84490060478a38802cc873b172e9',build_vignettes = T,upgrade_dependencies=F,build_opts = c("--no-resave-data", "--no-manual"))


