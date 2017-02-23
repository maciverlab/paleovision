#making sure all packages are installed

list.of.packages <- c("ape",
                      "phytools",
                      "geiger",
                      "denstrip",
                      "ouch",
                      "bayou",
                      "MASS",
                      "nlme",
                      "paleotree",
                      "strap")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages, repos='https://cloud.r-project.org')