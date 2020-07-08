library(glpkAPI)
library(sybil)
library(sybilSBML)
library(dplyr)
library(ggplot2)
library(egg)
library(readODS)
SYBIL_SETTINGS("SOLVER", "glpkAPI")

# STEP 5
setwd("~/github/StepDiet/")

#READ MODEL, FIND EXCHANGE REACTIONS, THEIR VALUES, AND THEIR POSITION
#PLEASE ADD THE ACTUAL DIRECTORY OF YOUR MODEL
core <- readSBMLmod("/home/corona/github/RoadMap/e_coli_core/e_coli_core.xml")
reactions_vcore <- grep(pattern = "EX_",
                           core@react_id, value = T, fixed = T)
reactions_pcore <- grep(pattern = "EX_",
                           core@react_id, value = F, fixed = T)


#LOAD THE DIET (STEP 6), SELECT ITEMS TO EXCLUDE NA VALUES
diet <- read_ods("S1_File/LB_Steps_autolysed_yeast_nofatty acids.ods",
                  sheet = "STEP 6")
diet_ex <- as.character(diet[, 19])[1:54]
diet_val <- as.numeric(diet[, 20])[1:54]

#FIND THE DIET'S EXCHANGER REACTIONS IN THE MODEL AND THEIR VALUES
diet_pos_core <- which(diet_ex %in% reactions_vcore, arr.ind = T)
diet_ex_core <- diet_ex[diet_pos_core]
diet_val_core <- -1 * diet_val[diet_pos_core]

#REMOVING ORIGINAL DIET
core_new <- changeBounds(model = core,
                            lb = rep(x = 0, length(reactions_vcore)),
                            ub = rep(x = 1000, length(reactions_vcore)),
                            react = reactions_vcore)

#APPLYING STEP 6 DIET
core_new <- changeBounds(model = core_new,
                            lb = diet_val_core,
                            ub = rep(x = 1000, length(diet_val_core)),
                            react = diet_ex_core)

#OPTIMIZATION WITH THE ORIGINAL SETTING
optimizeProb(core_new)

#check acetate production
 acetat <- grep(pattern = "EX_ac(e)",x = core_new@react_id, fixed = T)
 acetat_opt <- optimizeProb(core_new, algorithm = "fba", retOptSol = FALSE)
 acetat_opt$fluxes[acetat] # flux of acetate
 
