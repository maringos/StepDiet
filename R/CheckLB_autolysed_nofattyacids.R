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
iJO1366 <- readSBMLmod("/home/corona/github/RoadMap/iJO1366/iJO1366.xml")
reactions_viJO1366 <- grep(pattern = "EX_",
                           iJO1366@react_id, value = T, fixed = T)
reactions_piJO1366 <- grep(pattern = "EX_",
                           iJO1366@react_id, value = F, fixed = T)


#LOAD THE DIET (STEP 5), SELECT ITEMS TO EXCLUDE NA VALUES
diet <- read_ods("S1_File/LB_Steps_autolysed_yeast_nofatty acids.ods",
                  sheet = "STEP 5")
diet_ex <- as.character(diet[, 19])[1:53]
diet_val <- as.numeric(diet[, 20])[1:53]

#FIND THE DIET'S EXCHANGER REACTIONS IN THE MODEL AND THEIR VALUES
diet_pos_iJO1366 <- which(diet_ex %in% reactions_viJO1366, arr.ind = T)
diet_ex_iJO1366 <- diet_ex[diet_pos_iJO1366]
diet_val_iJO1366 <- -1 * diet_val[diet_pos_iJO1366]

#OPTIMIZATION WITH THE ORIGINAL SETTING
optimizeProb(iJO1366)

#REMOVING ORIGINAL DIET
iJO1366_new <- changeBounds(model = iJO1366,
                            lb = rep(x = 0, length(reactions_viJO1366)),
                            ub = rep(x = 1000, length(reactions_viJO1366)),
                            react = reactions_viJO1366)

#APPLYING STEP 5 DIET
iJO1366_new <- changeBounds(model = iJO1366_new,
                            lb = diet_val_iJO1366,
                            ub = rep(x = 1000, length(diet_val_iJO1366)),
                            react = diet_ex_iJO1366)

#EXCLUDE FLUXES FROM/TO WT BIOMASS REACTION
iJO1366_new <- changeBounds(model = iJO1366_new, lb = 0,
                            react = iJO1366_new@react_id[[7]], ub = 0)

#OPTIMIZATION
optimizeProb(iJO1366_new)
# RESULTS MUST SHOW "value of objective function (fba): 0.000000"
# RESULTS MUST SHOW "solution is optimal" IN iJO1366_new

# STEP 6
#LOAD THE DIET (STEP 6), SELECT ITEMS TO EXCLUDE NA VALUES
diet_STEP6 <- read_ods("S1_File/LB_Steps_autolysed_yeast_nofatty acids.ods",
                       sheet = "STEP 6")
diet_ex_STEP6 <- as.character(diet_STEP6[, 22])[1:54]
diet_val_STEP6 <- as.numeric(diet_STEP6[, 23])[1:54]

#FIND EXCHANGE REACTIONS WITH MINIMUM NEGATIVE REDUCED COSTS / FUNCTION 1
#KEEP ADDING THE FOUND COMPOUNDS UNTIL THERE IS GROWTH
GM_compounds_added_until_growth <- function(model, cutoff_for_growth) {
  vector.exchanges <- vector() # initial vector to save the exchange reactions
  mpd <- model # model
  vec.reactions <-  grep(pattern = "EX_", x = mpd@react_id, fixed = T,
                         value = T) # report all exchange reactions
  vec.reactions.positions <-  grep(pattern = "EX_",
                                   x = mpd@react_id, fixed = T,
                      value = F)  # report positions of all exchangre reactions
  # optimization results - as seen in example of sybil documentation
  x2 <- optimizeProb(mpd, poCmd = list("getRedCosts"))
  vd <- vector()  # interim vector to save the found exchangw reactions
  while (x2@lp_obj <= cutoff_for_growth) { # 1e-6 cut-off value
    x3 <- postProc(x2) # as seen in example of sybil documentation
    x4 <- x3@pa[[1]] # as seen in example of sybil documentation
    # make a dataframe to save the results - formatting
    df <- matrix(nrow = length(vec.reactions), ncol = 2)
    df[, 1] <- vec.reactions
    df1 <- as.data.frame(df)
    # only the reduced costs of the exchange reactions
    df1$V2 <- as.numeric(x4[vec.reactions.positions])
    df2 <- df1[abs((df1[, 2])) > (1e-6), ] # default cutoff
    df2 <- as.data.frame(df2)
    # alphabetically ordered
    dfmin <- sort(df2$V1[which(df2$V2 == min(df2$V2))])[[1]]
    #optimization - preparing the results for the next round of the loop
    #add -1 (assumption)
    mpd <- changeBounds(model = mpd, lb = -1, ub = 1000,
                        react = as.character(dfmin))
    x2 <- optimizeProb(mpd, poCmd = list("getRedCosts"))
    vd <- unique(c(vd, as.character(dfmin)))
  }
  cat("These compounds have been added until growth achieved: \n")
  print(vd)
  #initial vector to save the exchange reactions
  ve <- unique(c(vector.exchanges, vd))
  return(ve)
  }

#FIND EXCHANGE REACTIONS WITH MINIMUM NEGATIVE REDUCED COSTS / FUNCTION 2
#RETURNS A DATAFRAME WITH THE REDUCED COSTS AFTER ONE OPTIMIZATION
GM_reduced_costs_after_1_optimization <- function(model) {
  mpd <- model # model
  # report all exchangre reactions
  vec.reactions <-  grep(pattern = "EX_", x = mpd@react_id,
                         fixed = T, value = T)
  # report positions of all exchange reactions
  vec.reactions.positions <-  grep(pattern = "EX_", x = mpd@react_id, fixed = T,
                                   value = F)
  # optimization results - as seen in example of sybil documentation
  x2 <- optimizeProb(mpd, poCmd = list("getRedCosts"))
    x3 <- postProc(x2) # as seen in example of sybil documentation
    x4 <- x3@pa[[1]] # as seen in example of sybil documentation
    # make a dataframe to save the results - formatting
    df <- matrix(nrow = length(vec.reactions), ncol = 2)
    df[, 1] <- vec.reactions
    df1 <- as.data.frame(df)
    # only the reduced costs of the exchange reactions
    df1$V2 <- as.numeric(x4[vec.reactions.positions])
    df2 <- df1[abs((df1[, 2])) > (1e-6), ] # default cutoff
    df2 <- as.data.frame(df2)
    colnames(df2)[[1]] <- "Exchange Reaction"
    colnames(df2)[[2]] <- "Reduced Costs"
    # ordering
    print(df2[order(df2[, 2], decreasing = T), ])
    return(df2[order(df2[, 2], decreasing = T), ])
    }

# TROUBLESHOOTING FOR iJO1366_new MODEL
iJO1366_new_df_c <- GM_reduced_costs_after_1_optimization(model = iJO1366_new)
iJO1366_new_df_g <- GM_compounds_added_until_growth(model = iJO1366_new,
                                                    cutoff_for_growth = 1e-6)
# APPLYING STEP 6 DIET
iJO1366_new <- changeBounds(model = iJO1366_new,
                               lb = -diet_val_STEP6,
                               ub = rep(x = 1000, length(diet_ex_STEP6)),
                               react = diet_ex_STEP6)
# OPTIMIZATION
optimizeProb(iJO1366_new) # BINGO
iJO1366_new_optimized <- optimizeProb(iJO1366_new) # BINGO


# SUPPLEMENTATION EXPERIMENT
plot.supplement <- function(model, reactions_position,
                            reaction_value, add_percentage, add_value,
                            cutoff_percentage, exclude=NA) {

# dataframe that contains the original lb value,
# lb + add, add is positive
# lb - lb*add, add is positive, %
df.reactions <- data.frame(reactions = reaction_value,
                           original_lb = model@lowbnd[reactions_position],
                    lbx10 = (model@lowbnd[reactions_position] +
                    (add_percentage * model@lowbnd[reactions_position] / 100)),
                    lbp10 = (model@lowbnd[reactions_position] - add_value),
                    original_result = NA, x10_result = NA, p10_result = NA)
# exclude exchange reaction with lower bound 0
df.reactions <- filter(df.reactions, df.reactions$original_lb != 0)

#Across the reactions
for (i in 1:length(df.reactions$reactions)) {
  #Across the default values, + add, * add
  for (k in 2:4) {
  model_new <- changeBounds(model = model,
                            lb = df.reactions[[i, k]],
                            ub = 1000,
                            react = as.character(df.reactions$reactions[[i]]))
 res <- optimizeProb(model_new)
 df.reactions[[i, k + 3]] <- res@lp_obj
    }
  }
# vector to save cases with absolute relative difference < cutoff
posi.v <- vector()
for (i in 1:length(df.reactions$reactions)) {
  if (abs((df.reactions[[i, 6]] - df.reactions[[i, 5]]) / df.reactions[[i, 5]])
      < ((cutoff_percentage) / 100)) {
    posi.v <- c(posi.v, i)
    }
  }

for (i in 1:length(df.reactions$reactions)) {
  if (abs((df.reactions[[i, 7]] - df.reactions[[i, 5]]) / df.reactions[[i, 5]])
      < ((cutoff_percentage) / 100)) {
    posi.v <- c(posi.v, i)
    }
}
#remove such cases
df.reactionscorrected <- df.reactions[-posi.v, ]
# the % relative changes for +
df.reactionscorrected$percentage <- ((df.reactionscorrected$x10_result -
                                       df.reactionscorrected$original_result)
                                   / df.reactionscorrected$original_result) * 100
# the % relative changes for *
df.reactionscorrected$ADD <- ((df.reactionscorrected$p10_result -
                                df.reactionscorrected$original_result)
                             / df.reactionscorrected$original_result) * 100
#see below about the function subfunction
endsubfunction <- subfunction(model = model,
                              reactions_position = reactions_position,
                              reaction_value = reaction_value,
                              add_percentage = add_percentage,
                              cutoff_percentage = cutoff_percentage,
                              add_value = add_value, exclude = exclude,
                              df.reactionscorrected = df.reactionscorrected)
return(endsubfunction)
}

subfunction <- function(model,
                        reactions_position,
                        reaction_value,
                        add_percentage,
                        cutoff_percentage,
                        add_value, exclude,
                        df.reactionscorrected) {
  # data formating
  # new dataframe to collect * cases
  df.x10 <- data.frame(reaction = df.reactionscorrected$reactions,
                       rel_difference = df.reactionscorrected$percentage,
                       group = paste("x+", add_percentage, sep = ""))
  # if exclude is not NA
  df.x10 <- filter(df.x10, (df.x10$reaction %in% exclude) == F)

  # new dataframe to collect + cases
  df.p10 <- data.frame(reaction = df.reactionscorrected$reactions,
                       rel_difference = df.reactionscorrected$ADD,
                       group = paste("p+", add_value, sep = ""))
  # if exclude is not NA
  df.p10 <- filter(df.p10, (df.x10$reaction %in% exclude) == F)

  #replace exchange reaction with compounds
  #remove "exchange tag"
  for (i in 1:length(df.x10$reaction)) {
    df.x10[i, 4] <- strsplit(x = model@react_name[which(model@react_id ==
                                          as.character(df.x10$reaction[[i]]))],
                             split = " ")[[1]][1]
  }

  for (i in 1:length(df.p10$reaction)) {
    df.p10[i, 4] <- strsplit(x = model@react_name[which(model@react_id ==
                                          as.character(df.p10$reaction[[i]]))],
                             split = " ")[[1]][1]
  }
  #renaming
  colnames(df.x10)[4] <- "Compounds"
  df.x10$plot <- paste("Relative Supplementation (+ ",
                       add_percentage, "%)", sep = "")
  colnames(df.p10)[4] <- "Compounds"
  df.p10$plot <- paste("Absolute Supplementation (+ ",
                       add_value, "mM)", sep = "")

  #creating a common list
  retp <- vector("list", 3)
  retp[[1]] <- df.x10
  retp[[2]] <- df.p10
  retp[[3]] <- df.reactionscorrected
  return(retp)
}


#APPLICATION
plot <- plot.supplement(model = iJO1366_new,
                        reactions_position = reactions_piJO1366,
                        reaction_value = reactions_viJO1366,
                        add_percentage = 30,
                        cutoff_percentage = 0.0001, #corresponds to 1e-6 (not %)
                        add_value = 1)

#FORMATING
plot[[1]]$Category <- c("Nucleobases", "Amino acids", "Amino acids",
                        "Amino acids", "Amino acids", "Nucleobases",
                        "Amino acids", "Amino acids",
                        "Amino acids", "Amino acids", "Nucleobases",
                        "Carbohydrates", "Inorganic nutrients", "Amino acids",
                        "Amino acids", "Amino acids", "Nucleobases")
plot[[2]]$Category <- c("Nucleobases", "Amino acids", "Amino acids",
                        "Amino acids", "Amino acids", "Nucleobases",
                        "Amino acids", "Amino acids",
                        "Amino acids", "Amino acids", "Nucleobases",
                        "Carbohydrates", "Inorganic nutrients", "Amino acids",
                        "Amino acids", "Amino acids", "Nucleobases")
#ADDING CARBOHYDRATES
plot[[2]][18, ] <- plot[[2]][12, ]
plot[[2]][19, ] <- plot[[2]][12, ]
plot[[1]][18, ] <- plot[[1]][12, ]
plot[[1]][19, ] <- plot[[1]][12, ]
#GLUCOSE
plot[[2]][18, 1] <- "EX_glc__D(e)"
plot[[2]][18, 4] <- "D-Glucose"
plot[[1]][18, 1] <- "EX_glc__D(e)"
plot[[1]][18, 4] <- "D-Glucose"
iJO1366_new1 <- changeBounds(model = iJO1366_new,
                             lb = -1,
                             ub = 1000,
                             react = "EX_glc__D(e)")
#OPTIMIZATION
gl <- (((optimizeProb(iJO1366_new1))@lp_obj - (optimizeProb(iJO1366_new))@lp_obj) /
         (optimizeProb(iJO1366_new))@lp_obj) * 100
plot[[2]][18, 2] <- gl
plot[[1]][18, 2] <- 0 # 0 as original value 0 (e.g 0 + 30%*0 = 0)

#SPECIAL CASE 20 mM of glucose
iJO1366_new1 <- changeBounds(model = iJO1366_new,
                             lb = -20,
                             ub = 1000,
                             react = "EX_glc__D(e)")
# ADDING 20 mM OF GLUCOSE TO THE SIMULATIONS -> RESULT IN YIELD
(((optimizeProb(iJO1366_new1))@lp_obj - (optimizeProb(iJO1366_new))@lp_obj) /
    (optimizeProb(iJO1366_new))@lp_obj) * 100

#ARABINOSE
plot[[2]][19, 4] <- "L-Arabinose"
plot[[2]][19, 1] <- "EX_arab__L(e)"
plot[[1]][19, 1] <- "EX_arab__L(e)"
plot[[1]][19, 4] <- "L-Arabinose"
iJO1366_new1 <- changeBounds(model = iJO1366_new,
                             lb = -1,
                             ub = 1000,
                             react = "EX_arab__L(e)")
gl <- (((optimizeProb(iJO1366_new1))@lp_obj - (optimizeProb(iJO1366_new))@lp_obj) /
         (optimizeProb(iJO1366_new))@lp_obj) * 100
plot[[2]][19, 2] <- gl
plot[[1]][19, 2] <-  0 # 0 as original value 0 (e.g 0 + 30%*0 = 0)

#ORDERING
plot[[1]]$Compounds <- factor(plot[[1]]$Compounds,
                              levels = plot[[1]]$Compounds
                              [order(plot[[1]]$Category, decreasing = T)])
plot[[2]]$Compounds <- factor(plot[[2]]$Compounds, levels = plot[[2]]$Compounds
                              [order(plot[[2]]$Category, decreasing = T)])

#PLOTTING WITHOUT OXYGEN
c.plot <- rbind(plot[[1]], plot[[2]])
c.plot1 <- c.plot
c.plot1 <- c.plot1[-which(c.plot1$reaction %in% c("EX_o2(e)")), ]

plottt <- ggplot(data = c.plot1, aes(x = Compounds, y = rel_difference,
                                     fill = Category))  + coord_flip() +
  facet_wrap(~plot)  +
  geom_bar(stat = "identity") +  theme(axis.text.x = element_text(angle = 0,
                                      vjust = 1, hjust = 1))  +
  labs(y = "Relative difference in predicted yield (%)") +
  labs(x = "Compounds")
plot(plottt)

#PLOTTING OF OXYGEN
c.plot2 <- rbind(plot[[1]], plot[[2]])
c.plot2 <- c.plot2[which(c.plot2$reaction %in% c("EX_o2(e)")), ]
colnames(c.plot2)[5] <- "Supplementation"
c.plot2[1, 5] <- "+ 30 %"
c.plot2[2, 5] <- "+ 1 mM"

plottt2 <- ggplot(data = c.plot2, aes(x = Compounds, y = rel_difference,
                                    fill = Supplementation)) + coord_flip() +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))  +
  labs(y = "Relative difference in predicted yield (%)") +
  labs(x = "Compounds") + scale_fill_manual(values = c("#E69F00", "#999999"))
plot(plottt2)

ggarrange(plottt, plottt2, labels = c("a", "b"), ncol = 1, nrow = 2,
          heights = c(2, 0.25))


# CALCULATIONS
# COMPOUNDS
plot[[3]]$reactions[c(1, 11, 6, 17)]
# MEAN OF THEIR VALUES
round(mean(plot[[3]]$original_lb[c(1, 11)]), digits = 2)
round(mean(plot[[3]]$original_lb[c(6, 17)]), digits = 2)
