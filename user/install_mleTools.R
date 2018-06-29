
# Need to install devtools, if not already there
library(devtools)

# remove any old versions of the package before installing
remove.packages('mleTools')

# Install the growthTools package!
devtools::install_github("ctkremer/mleTools",auth_token = '6af60cb683ef4cbab1d6b25e0ec2ed6b925c831b')

# Load package!
library(mleTools)

?grid.mle2
