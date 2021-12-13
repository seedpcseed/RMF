#' RMetaflow
#' Module: setup
#' Initiated March 2021
#' Patrick C. Seed
#' Stanley Manne Children's Research Institute
#' Ann and Robert H. Lurie Children's Hospital
#' Northwestern University
#' MIT license
#'

# Check and install needed libraries --------------------------------------
pkgs.to.install <- attachment::att_from_rscripts("R")
install.packages(pkgs.to.install, dependencies = TRUE)

# call libraries needed at start ------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))
library(tidyverse)
library(yaml)

# set project directory root ----------------------------------------------

source("R/utilities.R")
setwd(getProjDir())

# source needed functions -------------------------------------------------

getCommandSteps()

# make directory structure ------------------------------------------------

if( !commands$startup$completed ){
  directories = c("data", "qc", "classify", "assembly", "metadata", "log")
  sapply(directories, dir.create)

  commands$startup$completed = TRUE; updateYaml()

}

# set some global parameters

source("R/variables.R")

# get steps and set start point -------------------------------------------

STEP <<- 2
# printGlue("#---- STEP set as {STEP}", log = TRUE)
getCurrentStep()

# run through uncompleted steps -------------------------------------------

while(STEP <= length(commands)){
  current.step = commands[STEP]
  printGlue("#---- Current step set to {STEP}", log = TRUE)
  printGlue("#---- Checking if directory needed for {dirNormalize(names(current.step))}", log = TRUE)
  dirIf(dirNormalize(names(current.step)))
  stepIf(current.step[[1]]$source)
  ifelse( current.step$type == "slurm", slurm <<-TRUE, slurm<<-FALSE)
  printGlue("#---- Calling {current.step$call}. SLURM = {slurm}", log = TRUE)
  source.complete = FALSE
  source(current.step$call)
  if( !source.complete ){
    Sys.sleep(1)
  }
  STEP <- STEP + 1
  printGlue("#---- Updated STEP to {STEP}", log = TRUE)
}

printGlue("All steps of RMetaflow were completed.

          If you think there are more steps to be finished,
          please review your metaflow-steps.yml file in the
          folder 'yaml'", log = TRUE)
q('no')
