
# Notes on useful commands for building/installing package

# Loading package development tools:
#library(devtools)

# at setup:
#devtools::use_readme_rmd()
#devtools::use_build_ignore("NEWS.md")
#devtools::use_build_ignore("TODO.R")

# during active development
#devtools::load_all()  # re-load all files in /R

# devtools::build_vignettes()


# To merge master and develop branch:
# - in terminal, run: git merge develop
# - resolve conflicts as necesary (indicated by orange coloring in Git pane of R Studio)
# - commit files, push to GitHub

# To increment the version number in the description file:
# - see use_version() function

# Track updates in versions within the NEWS.md file, check that it's readable via news()