source("R/packages.R")  # Load all the packages you need.
source("R/functions.R") # Load all the functions into your environment.

source("R/plot-alt.R")
source("R/plot-alt-util.R")
source("R/plot-alt-extra.R")
source("R/plotMaker.R")

source("R/plan.R")      # Build your workflow plan data frame.

make(my_plan_immediately)
