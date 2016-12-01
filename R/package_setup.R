# This is to install all necessary R packages
# On the Alcedo server I am using the homebrew version of R
# /home/jeremias/.linuxbrew/Cellar/r/3.3.2/bin/R

packages <- c("lme4", "lattice", "tidyverse", "cowplot")
install.packages(packages)
