#Measuring Political Polarization: The Cluster-Polarization Coefficient
#Isaac D. Mehlhaff
#February 2021
###################################################################################
#setup
###################################################################################

#install and load necessary packages

packages <- c("haven", "labelled", "pbapply", "parallel", "mixtools", "Metrics",
              "gdata", "psych", "cluster", "data.table", "GGally", "MASS",
              "animation", "tidyverse", "devtools")

for (i in packages) {
  ifelse(i %in% installed.packages(), print("Package already installed"),
         install.packages(i))
}

library(haven)
library(labelled)
library(pbapply)
library(parallel)
library(mixtools)
library(Metrics)
library(gdata)
library(psych)
library(cluster)
library(data.table)
library(GGally)
library(MASS)
library(animation)
library(tidyverse)
library(devtools)

install_github("imehlhaff/CPC@main")
library(CPC)

#define necessary settings

setwd(") #set to whatever directory contains data files
theme_set(theme_bw(base_size = 22))
cores <- detectCores()

#import data

nominate <- as.data.frame(read_csv("HSall_members.csv"))
inequality <- as.data.frame(read_csv("inequality.csv"))

###################################################################################
#functions
###################################################################################

#write difference function (univariate)

diff.uni <- function(x, y, k) {
  input <- na.omit(matrix(x, ncol = y))
  
  if(length(unique(input)) < k){
    return(NA)
  }
  
  else {
    output_kmeans <- kmeans(x = input, centers = k, nstart = 30)
    diff_1 <- abs(as.numeric(output_kmeans$centers[1, 1] -
                               output_kmeans$centers[2, 1]))
    diff_2 <- ifelse(k %in% c(3, 4),
                     abs(as.numeric(output_kmeans$centers[1, 1] -
                                      output_kmeans$centers[3, 1])), NA)
    diff_3 <- ifelse(k %in% c(3, 4),
                     abs(as.numeric(output_kmeans$centers[2, 1] -
                                      output_kmeans$centers[3, 1])), NA)
    diff_4 <- ifelse(k == 4,
                     abs(as.numeric(output_kmeans$centers[1, 1] -
                                      output_kmeans$centers[4, 1])), NA)
    diff_5 <- ifelse(k == 4,
                     abs(as.numeric(output_kmeans$centers[2, 1] -
                                      output_kmeans$centers[4, 1])), NA)
    diff_6 <- ifelse(k == 4,
                     abs(as.numeric(output_kmeans$centers[3, 1] -
                                      output_kmeans$centers[4, 1])), NA)
    diff <- mean(c(diff_1, diff_2, diff_3, diff_4, diff_5, diff_6), na.rm = TRUE)
    return(diff)
  }
}    

#write difference function (bivariate)

diff.bi <- function(x, y, k) {
  input <- na.omit(matrix(x, ncol = y))
  
  if(length(unique(input)) < k){
    return(NA)
  }
  
  else {
    output_kmeans <- kmeans(x = input, centers = k, nstart = 30)
    diff_1 <- as.numeric(sqrt((output_kmeans$centers[1, 1] -
                                 output_kmeans$centers[1, 2])^2 +
                                (output_kmeans$centers[2, 1] -
                                   output_kmeans$centers[2, 2])^2))
    diff_2 <- ifelse(k %in% c(3, 4),
                     as.numeric(sqrt((output_kmeans$centers[1, 1] -
                                        output_kmeans$centers[1, 2])^2 +
                                       (output_kmeans$centers[3, 1] -
                                          output_kmeans$centers[3, 2])^2)), NA)
    diff_3 <- ifelse(k %in% c(3, 4),
                     as.numeric(sqrt((output_kmeans$centers[2, 1] -
                                        output_kmeans$centers[2, 2])^2 +
                                       (output_kmeans$centers[3, 1] -
                                          output_kmeans$centers[3, 2])^2)), NA)
    diff_4 <- ifelse(k == 4,
                     as.numeric(sqrt((output_kmeans$centers[1, 1] -
                                        output_kmeans$centers[1, 2])^2 +
                                       (output_kmeans$centers[4, 1] -
                                          output_kmeans$centers[4, 2])^2)), NA)
    diff_5 <- ifelse(k == 4,
                     as.numeric(sqrt((output_kmeans$centers[2, 1] -
                                        output_kmeans$centers[2, 2])^2 +
                                       (output_kmeans$centers[4, 1] -
                                          output_kmeans$centers[4, 2])^2)), NA)
    diff_6 <- ifelse(k == 4,
                     as.numeric(sqrt((output_kmeans$centers[3, 1] -
                                        output_kmeans$centers[3, 2])^2 +
                                       (output_kmeans$centers[4, 1] -
                                          output_kmeans$centers[4, 2])^2)), NA)
    diff <- mean(c(diff_1, diff_2, diff_3, diff_4, diff_5, diff_6), na.rm = TRUE)
    return(diff)
  }
}    

#write function to calculate Euclidean distance in two dimensions

Euclidean <- function(x) {
  x <- as.data.frame(na.omit(x))
  colnames(x) <- c("V1", "V2")
  x <- mutate(x, x_mean = rep(base::mean(x$V1), nrow(x)),
              y_mean = rep(base::mean(x$V2), nrow(x)),
              distance = sqrt((x$V1 - x_mean)^2 + (x$V2 - y_mean)^2))
  return(sum(x$distance)/nrow(x))
}

#write function to create kmeans visualization (adapted from kmeans.ani function in animation package - makes plots more readable)

kmeans.ani.alt <- function(x = cbind(X1 = runif(50), X2 = runif(50)), centers = 3,
                           hints = c("Move centers!", "Find cluster?"), pch = 1:3,
                           pch_centers = 1:3, col = 1:3) {
  x = as.matrix(x)
  ocluster = sample(centers, nrow(x), replace = TRUE)
  if (length(centers) == 1) 
    centers = x[sample(nrow(x), centers), ]
  else centers = as.matrix(centers)
  numcent = nrow(centers)
  dst = matrix(nrow = nrow(x), ncol = numcent)
  j = 1
  pch = rep(pch, length = numcent)
  col = rep(col, length = numcent)
  for (j in 1:ani.options("nmax")) {
    dev.hold()
    plot(x, pch = pch[ocluster], col = col[ocluster], panel.first = grid(),
         xaxt = "n", yaxt = "n", ann = FALSE, cex = 2)
    mtext(hints[1], 4)
    points(centers, pch = pch_centers[1:numcent], cex = 4, lwd = 2, 
           col = col[1:numcent])
    ani.pause()
    for (i in 1:numcent) {
      dst[, i] = sqrt(apply((t(t(x) - unlist(centers[i,])))^2, 1, sum))
    }
    ncluster = apply(dst, 1, which.min)
    plot(x, type = "n", xaxt = "n", yaxt = "n", ann = FALSE, cex = 2)
    mtext(hints[2], 4)
    grid()
    ocenters = centers
    for (i in 1:numcent) {
      xx = subset(x, ncluster == i)
      polygon(xx[chull(xx), ], density = 10, col = col[i], lty = 2)
      points(xx, pch = pch[i], col = col[i], cex = 2)
      centers[i, ] = apply(xx, 2, mean)
    }
    points(ocenters, cex = 4, col = col[1:numcent], pch = pch_centers[1:numcent],
           lwd = 2)
    ani.pause()
    if (all(centers == ocenters)) 
      break
    ocluster = ncluster
  }
  invisible(list(cluster = ncluster, centers = centers))
}

###################################################################################
#data cleaning
###################################################################################

#remove observations from uncompleted 117th Congress and select necessary variables

nominate <- filter(nominate, congress != 117) %>%
  dplyr::select(congress, chamber, nominate_dim1, nominate_dim2, nokken_poole_dim1,
                nokken_poole_dim2)

#break apart into House and Senate; eliminates President observations

nom_house <- filter(nominate, chamber == "House")
nom_senate <- filter(nominate, chamber == "Senate")

#convert all inequality measures to decimals

inequality <- mutate(inequality,
                     `Gini coefficient` = `Gini coefficient`*0.01,
                     `Top 1% income` = `Top 1% income`*0.01,
                     `Top 1% wealth` = `Top 1% wealth`*0.01)

###################################################################################
#calculation of polarization estimates
###################################################################################

#calculate polarization estimates for each chamber using three measures and scale measures for presentation on same plot

house_polarization <- group_by(nom_house, congress) %>%
  dplyr::summarize(nom_diff = diff.uni(nominate_dim1, y = 1, k = 2),
                   nom_var = var(nominate_dim1, na.rm = TRUE),
                   nom_kurt = -(kurtosi(nominate_dim1, na.rm = TRUE)),
                   nom_CPC = CPC(nominate_dim1, 2, "kmeans", adjust = TRUE,
                                 nstart = 30),
                   np_diff = diff.uni(nokken_poole_dim1, y = 1, k = 2),
                   np_var = var(nokken_poole_dim1, na.rm = TRUE),
                   np_kurt = -(kurtosi(nokken_poole_dim1, na.rm = TRUE)),
                   np_CPC = CPC(nokken_poole_dim1, 2, "kmeans", adjust = TRUE,
                                nstart = 30)) %>%
  mutate(nom_diff_scale = scales::rescale(nom_diff, to = c(0, 1)),
         nom_var_scale = scales::rescale(nom_var, to = c(0, 1)),
         nom_kurt_scale = scales::rescale(nom_kurt, to = c(0, 1)),
         nom_CPC_scale = scales::rescale(nom_CPC, to = c(0, 1)),
         np_diff_scale = scales::rescale(np_diff, to = c(0, 1)),
         np_var_scale = scales::rescale(np_var, to = c(0, 1)),
         np_kurt_scale = scales::rescale(np_kurt, to = c(0, 1)),
         np_CPC_scale = scales::rescale(np_CPC, to = c(0, 1)),
         chamber = "House of Representatives")

senate_polarization <- group_by(nom_senate, congress) %>%
  dplyr::summarize(nom_diff = diff.uni(nominate_dim1, y = 1, k = 2),
                   nom_var = var(nominate_dim1, na.rm = TRUE),
                   nom_kurt = -(kurtosi(nominate_dim1, na.rm = TRUE)),
                   nom_CPC = CPC(nominate_dim1, 2, "kmeans", adjust = TRUE,
                                 nstart = 30),
                   np_diff = diff.uni(nokken_poole_dim1, y = 1, k = 2),
                   np_var = var(nokken_poole_dim1, na.rm = TRUE),
                   np_kurt = -(kurtosi(nokken_poole_dim1, na.rm = TRUE)),
                   np_CPC = CPC(nokken_poole_dim1, 2, "kmeans", adjust = TRUE,
                                nstart = 30)) %>%
  mutate(nom_diff_scale = scales::rescale(nom_diff, to = c(0, 1)),
         nom_var_scale = scales::rescale(nom_var, to = c(0, 1)),
         nom_kurt_scale = scales::rescale(nom_kurt, to = c(0, 1)),
         nom_CPC_scale = scales::rescale(nom_CPC, to = c(0, 1)),
         np_diff_scale = scales::rescale(np_diff, to = c(0, 1)),
         np_var_scale = scales::rescale(np_var, to = c(0, 1)),
         np_kurt_scale = scales::rescale(np_kurt, to = c(0, 1)),
         np_CPC_scale = scales::rescale(np_CPC, to = c(0, 1)),
         chamber = "Senate")

#combine House and Senate polarization estimates into single data frame

nom_polarization <- rbind(house_polarization, senate_polarization)

###################################################################################
#plots of DW-NOMINATE polarization estimates
###################################################################################

#divide polarization data frame by estimate type and mean/se, collapse measures into single column, and add column denoting estimate type

nom <- dplyr::select(nom_polarization, congress, chamber, starts_with("nom")) %>%
  pivot_longer(cols = c("nom_diff", "nom_var", "nom_kurt", "nom_CPC",
                        "nom_diff_scale", "nom_var_scale", "nom_kurt_scale",
                        "nom_CPC_scale"), names_to = "measure",
               values_to = "value") %>%
  mutate(type = "NOMINATE",
         measure = substr(measure, 5, nchar(measure)),
         scaled = ifelse(measure %in% c("diff_scale", "var_scale", "kurt_scale",
                                        "CPC_scale"), 1, 0),
         measure = ifelse(measure %in% c("diff", "diff_scale"), "difference",
                          ifelse(measure %in% c("var", "var_scale"), "variance",
                                 ifelse(measure %in% c("kurt", "kurt_scale"),
                                        "kurtosis",
                                        ifelse(measure %in% c("CPC", "CPC_scale"),
                                               "CPC", NA)))))

np <- dplyr::select(nom_polarization, congress, chamber, starts_with("np")) %>%
  pivot_longer(cols = c("np_diff", "np_var", "np_kurt", "np_CPC", "np_diff_scale",
                        "np_var_scale", "np_kurt_scale", "np_CPC_scale"),
               names_to = "measure", values_to = "value") %>%
  mutate(type = "Nokken-Poole",
         measure = substr(measure, 4, nchar(measure)),
         scaled = ifelse(measure %in% c("diff_scale", "var_scale", "kurt_scale",
                                        "CPC_scale"), 1, 0),
         measure = ifelse(measure %in% c("diff", "diff_scale"), "difference",
                          ifelse(measure %in% c("var", "var_scale"), "variance",
                                 ifelse(measure %in% c("kurt", "kurt_scale"),
                                        "kurtosis",
                                        ifelse(measure %in% c("CPC", "CPC_scale"),
                                               "CPC", NA)))))

#merge data frames back together and relevel faceting variables

nom_polarization <- rbind(nom, np)

nom_polarization$type <- factor(nom_polarization$type,
                                levels = c("NOMINATE", "Nokken-Poole"))
nom_polarization$measure <- factor(nom_polarization$measure,
                                   levels = c("difference", "variance", "kurtosis",
                                              "CPC"))

#add labels to data frame

nom_polarization$label <- NA
nom_polarization$y <- 0.97 #0.86 for bidimensional

nom_polarization$label[nom_polarization$chamber == "House of Representatives" &
                         nom_polarization$congress == 50] <- "Gilded\nAge"
nom_polarization$label[nom_polarization$chamber == "House of Representatives" &
                         nom_polarization$congress == 80] <- "Textbook\nCongress"
nom_polarization$label[nom_polarization$chamber == "House of Representatives" &
                         nom_polarization$congress == 109] <- "Modern\nPolarization"

#plot unscaled CPC (Figure 5)

ggplot(subset(nom_polarization, scaled == 0 & measure == "CPC"),
       aes(x = congress, y = value, linetype = type)) +
  geom_rect(aes(xmin = 45, xmax = 56, ymin = -Inf, ymax = Inf), fill = "gray",
            alpha = 0.01) +
  geom_rect(aes(xmin = 73, xmax = 88, ymin = -Inf, ymax = Inf), fill = "gray",
            alpha = 0.01) +
  geom_rect(aes(xmin = 103, xmax = 116, ymin = -Inf, ymax = Inf), fill = "gray",
            alpha = 0.01) +
  geom_text(aes(y = y, label = label), size = 3) +
  geom_line() +
  facet_wrap(~ chamber, ncol = 1) +
  labs(x = "Congress", y = "CPC Estimate", linetype = "Ideology Estimate")

#plot scaled measures (Figure S13)

ggplot(subset(nom_polarization, scaled == 1),
       aes(x = congress, y = value, linetype = type)) +
  geom_line() +
  facet_grid(measure ~ chamber) +
  labs(x = "Congress", y = "Polarization Estimates", linetype = "Ideology Estimate")

###################################################################################
#income inequality analysis
###################################################################################

#separate polarization data by chamber, merge in inequality data, and retain only regular NOMINATE estimates

house_short <- dplyr::filter(house_polarization, congress %in% c(63:114)) %>%
  dplyr::select(1:5) %>%
  dplyr::rename(difference = nom_diff, variance = nom_var, kurtosis = nom_kurt,
                CPC = nom_CPC) %>%
  inner_join(., inequality) %>%
  mutate(chamber = "House of Representatives")

senate_short <- dplyr::filter(senate_polarization, congress %in% c(63:114)) %>%
  dplyr::select(1:5) %>%
  dplyr::rename(difference = nom_diff, variance = nom_var, kurtosis = nom_kurt,
                CPC = nom_CPC) %>%
  inner_join(., inequality) %>%
  mutate(chamber = "Senate")

#combine chambers and shorten chamber and inequality names

polar_short <- rbind(house_short, senate_short)

polar_short <- mutate(polar_short,
                      chamber = ifelse(chamber == "House of Representatives",
                                       "House",
                                       ifelse(chamber == "Senate", "Senate", NA)))

#plot correlation between inequality and each polarization measure for each chamber (Figure 6)

theme_set(theme_bw(base_size = 21))

ggpairs(polar_short, columns = 5:8, mapping = aes(color = chamber),
        diag = list(continuous = wrap("densityDiag", alpha = 0.2)),
        upper = list(continuous = wrap("cor", stars = FALSE, size = 6)),
        lower = list(continuous = wrap("smooth", alpha = 0.3, size = 1))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#rescale all variables and pivot to longer format for plotting

house_short <- mutate(house_short,
                      `Gini coefficient` = scales::rescale(`Gini coefficient`,
                                                           to = c(0, 1)),
                      `Top 1% income` = scales::rescale(`Top 1% income`,
                                                        to = c(0, 1)),
                      `Top 1% wealth` = scales::rescale(`Top 1% wealth`,
                                                        to = c(0, 1)),
                      difference = scales::rescale(difference, to = c(0, 1)),
                      variance = scales::rescale(variance, to = c(0, 1)),
                      kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                      CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "measure_polar",
               values_to = "value_polar") %>%
  pivot_longer(cols = 2:4, names_to = "measure_ineq", values_to = "value_ineq")

senate_short <- mutate(senate_short,
                       `Gini coefficient` = scales::rescale(`Gini coefficient`,
                                                            to = c(0, 1)),
                       `Top 1% income` = scales::rescale(`Top 1% income`,
                                                         to = c(0, 1)),
                       `Top 1% wealth` = scales::rescale(`Top 1% wealth`,
                                                         to = c(0, 1)),
                       difference = scales::rescale(difference, to = c(0, 1)),
                       variance = scales::rescale(variance, to = c(0, 1)),
                       kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                       CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "measure_polar",
               values_to = "value_polar") %>%
  pivot_longer(cols = 2:4, names_to = "measure_ineq", values_to = "value_ineq")

#combine rescaled polarization data

polar_short <- rbind(house_short, senate_short)

#relevel polarization variable to get plots in correct order

polar_short$measure_polar <- factor(polar_short$measure_polar,
                                    levels = c("difference", "variance", "kurtosis",
                                               "CPC"))
polar_short$measure_ineq <- factor(polar_short$measure_ineq,
                                   levels = c("Gini coefficient", "Top 1% income",
                                              "Top 1% wealth"))

#plot polarization and Gini coefficient (Figure S14)

ggplot(subset(polar_short, measure_polar %in% c("difference", "variance",
                                                "kurtosis", "CPC") &
                measure_ineq == "Gini coefficient" & !is.na(value_ineq)),
       aes(x = congress)) +
  geom_line(aes(y = value_polar, linetype = "Estimated polarization")) +
  geom_line(aes(y = value_ineq, linetype = "Gini coefficient")) +
  facet_grid(measure_polar ~ chamber) +
  labs(x = "Congress", y = "Scaled Estimates of Polarization and Inequality",
       linetype = "Variable")

#plot polarization and top 1% income (Figure S15)

ggplot(subset(polar_short, measure_polar %in% c("difference", "variance",
                                                "kurtosis", "CPC") &
                measure_ineq == "Top 1% income" & !is.na(value_ineq)),
       aes(x = congress)) +
  geom_line(aes(y = value_polar, linetype = "Estimated polarization")) +
  geom_line(aes(y = value_ineq, linetype = "Top 1% income")) +
  facet_grid(measure_polar ~ chamber) +
  labs(x = "Congress", y = "Scaled Estimates of Polarization and Inequality",
       linetype = "Variable")

#plot polarization and top 1% wealth (Figure S16)

ggplot(subset(polar_short, measure_polar %in% c("difference", "variance",
                                                "kurtosis", "CPC") &
                measure_ineq == "Top 1% wealth" & !is.na(value_ineq)),
       aes(x = congress)) +
  geom_line(aes(y = value_polar, linetype = "Estimated polarization")) +
  geom_line(aes(y = value_ineq, linetype = "Top 1% wealth")) +
  facet_grid(measure_polar ~ chamber) +
  labs(x = "Congress", y = "Scaled Estimates of Polarization and Inequality",
       linetype = "Variable")

#resent plot base size

theme_set(theme_bw(base_size = 22))

###################################################################################
#k-means demonstration
###################################################################################

set.seed(999)

#generate bimodal data

data <- rmvnormmix(200, c(0.5, 0.5), 
                   matrix(data = c(0.25, 0.2500001, 0.75, 0.7500001), nrow = 2,
                          ncol = 2, byrow = TRUE), 
                   matrix(data = 0.125, nrow = 2, ncol = 2, byrow = TRUE))

#convert to data frame, rename columns, and filter out observations outside bounds

data <- as.data.frame(matrix(data = data, ncol = 2, byrow = FALSE))
colnames(data) <- c("x", "y")
data <- filter(data, x <= 1, x >= 0, y <= 1, y >= 0)

#plot raw data

plot(data$x, data$y, xaxt = "n", yaxt = "n", ann = FALSE, panel.first = grid(),
     cex = 2)

#initiate k-means animation

kmeans.ani.alt(data, centers = 2, hints = c(), pch = c(1, 2),
               pch_centers = c(16, 17), col = c("red", "blue"))

###################################################################################
#CPC vs. adjusted CPC simulation (increasing dimensions)
###################################################################################

#set seed

set.seed(123)

#generate 20 sets of simulated Guassian mixture distributions using parameters giving middling level of polarization, adding one dimension each time

CPC_data2 <- list()
CPC_data3 <- list()
CPC_data4 <- list()
CPC_compare2 <- c()
CPC_compare3 <- c()
CPC_compare4 <- c()
CPC_adj_compare2 <- c()
CPC_adj_compare3 <- c()
CPC_adj_compare4 <- c()

for (i in 1:20) {
  data2 <- rmvnormmix(1000, c(0.5, 0.5),
                      matrix(data = c(-3, 3), nrow = 2, ncol = i),
                      matrix(data = 1, nrow = 2, ncol = i))
  CPC_data2[[i]] <- data2
  
  data3 <- rmvnormmix(1000, rep(0.33, 3),
                      matrix(data = c(-3, 0, 3), nrow = 3, ncol = i),
                      matrix(data = 1, nrow = 3, ncol = i))
  CPC_data3[[i]] <- data3
  
  data4 <- rmvnormmix(1000, rep(0.25, 4),
                      matrix(data = c(-3, -1.5, 1.5, 3), nrow = 4, ncol = i),
                      matrix(data = 1, nrow = 4, ncol = i))
  CPC_data4[[i]] <- data4
}

#calculate CPC and adjusted CPC for each simulated distribution

for (i in 1:20) {
  output2 <- CPC(CPC_data2[[i]], 2, "kmeans", nstart = 30)
  output_adj2 <- CPC(CPC_data2[[i]], 2, "kmeans", adjust = TRUE, nstart = 30)
  
  output3 <- CPC(CPC_data3[[i]], 3, "kmeans", nstart = 30)
  output_adj3 <- CPC(CPC_data3[[i]], 3, "kmeans", adjust = TRUE, nstart = 30)
  
  output4 <- CPC(CPC_data4[[i]], 4, "kmeans", nstart = 30)
  output_adj4 <- CPC(CPC_data4[[i]], 4, "kmeans", adjust = TRUE, nstart = 30)
  
  CPC_compare2 <- c(CPC_compare2, output2)
  CPC_adj_compare2 <- c(CPC_adj_compare2, output_adj2)
  
  CPC_compare3 <- c(CPC_compare3, output3)
  CPC_adj_compare3 <- c(CPC_adj_compare3, output_adj3)
  
  CPC_compare4 <- c(CPC_compare4, output4)
  CPC_adj_compare4 <- c(CPC_adj_compare4, output_adj4)
}

#combine CPC and adjusted CPC calculations

CPC_compare <- matrix(c(CPC_compare2, CPC_compare3, CPC_compare4), ncol = 1)
CPC_adj_compare <- matrix(c(CPC_adj_compare2, CPC_adj_compare3, CPC_adj_compare4),
                          ncol = 1)

#combine CPC and adjusted CPC calculations into tidy data frame, coerce cluster numbers to factor

CPC_plot <- as.data.frame(cbind(rep(seq(1, 20), 3), rep(seq(2, 4), each = 20),
                                CPC_compare, CPC_adj_compare))

colnames(CPC_plot) <- c("Number of Dimensions", "Number of Clusters", "CPC",
                        "adj. CPC")

CPC_plot <- pivot_longer(CPC_plot, cols = c("CPC", "adj. CPC"),
                         names_to = "Measure", values_to = "Polarization Estimate")

CPC_plot$`Number of Clusters` <- as.factor(CPC_plot$`Number of Clusters`)

#plot CPC vs. adjusted CPC to visualize correction (Figure S1)

ggplot(CPC_plot, aes(x = `Number of Dimensions`, y = `Polarization Estimate`,
                     linetype = Measure, color = `Number of Clusters`)) +
  geom_line() +
  scale_linetype_manual(values = c("longdash", "solid"))

###################################################################################
#CPC vs. adjusted CPC simulation (increasing clusters)
###################################################################################

#set seed

set.seed(77)

#generate 19 sets of simulated Guassian mixture distributions, adding one cluster each time

CPC_data1 <- list()
CPC_data2 <- list()
CPC_compare1 <- c()
CPC_compare2 <- c()
CPC_adj_compare1 <- c()
CPC_adj_compare2 <- c()

for (i in 2:20) {
  data1 <- rnormmix(1000, rep(1/i, i), c(seq(1, i)/(i*2)), 0.01)
  CPC_data1[[i]] <- data1
  
  data2 <- rmvnormmix(1000, rep(1/i, i),
                      matrix(data = c(seq(1, i)/(i*2)), nrow = i, ncol = 2),
                      matrix(data = 0.01, nrow = i, ncol = 2))
  CPC_data2[[i]] <- data2
}

#calculate CPC and adjusted CPC for each simulated distribution

for (i in 2:20) {
  output1 <- CPC(CPC_data1[[i]], i, "kmeans", nstart = 30)
  output_adj1 <- CPC(CPC_data1[[i]], i, "kmeans", adjust = TRUE, nstart = 30)
  
  output2 <- CPC(CPC_data2[[i]], i, "kmeans", nstart = 30)
  output_adj2 <- CPC(CPC_data2[[i]], i, "kmeans", adjust = TRUE, nstart = 30)
  
  CPC_compare1 <- c(CPC_compare1, output1)
  CPC_adj_compare1 <- c(CPC_adj_compare1, output_adj1)
  
  CPC_compare2 <- c(CPC_compare2, output2)
  CPC_adj_compare2 <- c(CPC_adj_compare2, output_adj2)
}

#combine CPC and adjusted CPC calculations

CPC_compare <- matrix(c(CPC_compare1, CPC_compare2), ncol = 1)
CPC_adj_compare <- matrix(c(CPC_adj_compare1, CPC_adj_compare2), ncol = 1)

#combine CPC and adjusted CPC calculations into tidy data frame, coerce cluster numbers to factor

CPC_plot <- as.data.frame(cbind(rep(seq(2, 20), 2), rep(seq(1, 2), each = 19),
                                CPC_compare, CPC_adj_compare))

colnames(CPC_plot) <- c("Number of Clusters", "Number of Dimensions", "CPC",
                        "adj. CPC")

CPC_plot <- pivot_longer(CPC_plot, cols = c("CPC", "adj. CPC"),
                         names_to = "Measure", values_to = "Polarization Estimate")

CPC_plot$`Number of Dimensions` <- as.factor(CPC_plot$`Number of Dimensions`)

#plot CPC vs. adjusted CPC to visualize correction (Figure S2)

ggplot(CPC_plot, aes(x = `Number of Clusters`, y = `Polarization Estimate`,
                     linetype = Measure, color = `Number of Dimensions`)) +
  geom_line() +
  scale_linetype_manual(values = c("longdash", "solid"))

###################################################################################
#simulated univariate distributions for visual comparison
###################################################################################

#generate simulated bimodal distributions

set.seed(123)

variance_values <- c(0.5, 1, 1.5, 2)
means_values <- c(5, 4, 3, 2)
plots <- c()

for (i in variance_values) {
  for (j in means_values) {
    data <- rnormmix(1000, c(0.5, 0.5), c(-j, j), c(i, i))
    plots <- cbind(plots, data)
    base::print(j)
  }
  base::print(i)
}

#add plot labels and convert to tidy data frame

plots <- as.data.frame(plots)
plots_univariate <- c("5, 0.5", "4, 0.5", "3, 0.5", "2, 0.5", "5, 1", "4, 1", 
                      "3, 1", "2, 1", "5, 1.5", "4, 1.5", "3, 1.5", "2, 1.5",
                      "5, 2", "4, 2", "3, 2", "2, 2")
colnames(plots) <- plots_univariate
plots <- pivot_longer(plots, cols = everything()) %>%
  dplyr::arrange(name)

#add plot labels and specify levels for ordering plots

means_labels <- c(rep("means = (-2, 2)", 4000), rep("means = (-3, 3)", 4000),
                  rep("means = (-4, 4)", 4000), rep("means = (-5, 5)", 4000))
variance_labels <- rep(c(rep("sd = 0.5", 1000), rep("sd = 1", 1000),
                         rep("sd = 1.5", 1000), rep("sd = 2", 1000)), 4)

plots <- cbind(plots, means_labels, variance_labels)

plots$means_labels <- factor(plots$means_labels, levels = c("means = (-2, 2)", 
                                                            "means = (-3, 3)",
                                                            "means = (-4, 4)",
                                                            "means = (-5, 5)"))
plots$variance_labels <- factor(plots$variance_labels, levels = c("sd = 2",
                                                                  "sd = 1.5",
                                                                  "sd = 1",
                                                                  "sd = 0.5"))

#plot univariate distributions with variance, kurtosis, and CPC (Figure 2)

ggplot(plots, aes(x = value)) +
  geom_density() +
  facet_grid(variance_labels ~ means_labels) +
  labs(x = "Value", y = "Density")

###################################################################################
#Monte Carlo simulations for univariate distributions (2 clusters)
###################################################################################

#set seed

set.seed(123)

#generate Gaussian mixture distributions

MC_uni_ranvar <- list()
MC_uni_ranmean <- list()
MC_uni_variances <- c()
MC_uni_means <- c()

for (i in means_values) {
  for (j in 1:1000) {
    variance <- runif(1, min = 0.5, max = 2)
    data <- rnormmix(1000, c(0.5, 0.5), c(-i, i), c(variance, variance))
    MC_uni_ranvar[[length(MC_uni_ranvar) + 1]] <- data
    MC_uni_variances[[length(MC_uni_variances) + 1]] <- variance
  }
  base::print(i)
}

for (i in variance_values) {
  for (j in 1:1000) {
    mean <- runif(1, min = 2, max = 5)
    data <- rnormmix(1000, c(0.5, 0.5), c(-mean, mean), c(i, i))
    MC_uni_ranmean[[length(MC_uni_ranmean) + 1]] <- data
    MC_uni_means[[length(MC_uni_means) + 1]] <- mean
  }
  base::print(i)
}

#calculate various polarization scores for simulated distributions

MC_uni_ranvar_diff <- pblapply(MC_uni_ranvar,
                               function(x) diff.uni(x = x, y = 1, k = 2),
                               cl = cores)
MC_uni_ranvar_var <- pblapply(MC_uni_ranvar, var, cl = cores)
MC_uni_ranvar_kurt <- pblapply(MC_uni_ranvar, kurtosi, cl = cores)
MC_uni_ranvar_CPC <- pblapply(MC_uni_ranvar,
                              function(x) CPC(x, 2, "kmeans", adjust = TRUE,
                                              nstart = 30), cl = cores)

MC_uni_ranmean_diff <- pblapply(MC_uni_ranmean,
                                function(x) diff.uni(x = x, y = 1, k = 2),
                                cl = cores)
MC_uni_ranmean_var <- pblapply(MC_uni_ranmean, var, cl = cores)
MC_uni_ranmean_kurt <- pblapply(MC_uni_ranmean, kurtosi, cl = cores)
MC_uni_ranmean_CPC <- pblapply(MC_uni_ranmean,
                               function(x) CPC(x, 2, "kmeans", adjust = TRUE,
                                               nstart = 30), cl = cores)

#combine random variance polarization results into data frames, holding means constant within data frame for accurate rank calculations

MC_uni_ranvar1 <- as.data.frame(matrix(data = c(unlist(MC_uni_variances[1:1000]),
                                                unlist(MC_uni_ranvar_diff[1:1000]),
                                                unlist(MC_uni_ranvar_var[1:1000]),
                                                unlist(MC_uni_ranvar_kurt[1:1000]),
                                                unlist(MC_uni_ranvar_CPC[1:1000])),
                                       nrow = 1000))
colnames(MC_uni_ranvar1) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC_uni_ranvar2 <- as.data.frame(matrix(data = c(unlist(MC_uni_variances[1001:2000]),
                                                unlist(MC_uni_ranvar_diff[1001:2000]),
                                                unlist(MC_uni_ranvar_var[1001:2000]),
                                                unlist(MC_uni_ranvar_kurt[1001:2000]),
                                                unlist(MC_uni_ranvar_CPC[1001:2000])),
                                       nrow = 1000))
colnames(MC_uni_ranvar2) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC_uni_ranvar3 <- as.data.frame(matrix(data = c(unlist(MC_uni_variances[2001:3000]),
                                                unlist(MC_uni_ranvar_diff[2001:3000]),
                                                unlist(MC_uni_ranvar_var[2001:3000]),
                                                unlist(MC_uni_ranvar_kurt[2001:3000]),
                                                unlist(MC_uni_ranvar_CPC[2001:3000])),
                                       nrow = 1000))
colnames(MC_uni_ranvar3) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC_uni_ranvar4 <- as.data.frame(matrix(data = c(unlist(MC_uni_variances[3001:4000]),
                                                unlist(MC_uni_ranvar_diff[3001:4000]),
                                                unlist(MC_uni_ranvar_var[3001:4000]),
                                                unlist(MC_uni_ranvar_kurt[3001:4000]),
                                                unlist(MC_uni_ranvar_CPC[3001:4000])),
                                       nrow = 1000))
colnames(MC_uni_ranvar4) <- c("ranvar", "difference" ,"variance", "kurtosis", "CPC")

#combine random mean polarization results into data frames, holding variance constant within data frame for accurate rank calculations

MC_uni_ranmean1 <- as.data.frame(matrix(data = c(unlist(MC_uni_means[1:1000]),
                                                 unlist(MC_uni_ranmean_diff[1:1000]),
                                                 unlist(MC_uni_ranmean_var[1:1000]),
                                                 unlist(MC_uni_ranmean_kurt[1:1000]),
                                                 unlist(MC_uni_ranmean_CPC[1:1000])),
                                        nrow = 1000))
colnames(MC_uni_ranmean1) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

MC_uni_ranmean2 <- as.data.frame(matrix(data = c(unlist(MC_uni_means[1001:2000]),
                                                 unlist(MC_uni_ranmean_diff[1001:2000]),
                                                 unlist(MC_uni_ranmean_var[1001:2000]),
                                                 unlist(MC_uni_ranmean_kurt[1001:2000]),
                                                 unlist(MC_uni_ranmean_CPC[1001:2000])),
                                        nrow = 1000))
colnames(MC_uni_ranmean2) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

MC_uni_ranmean3 <- as.data.frame(matrix(data = c(unlist(MC_uni_means[2001:3000]),
                                                 unlist(MC_uni_ranmean_diff[2001:3000]),
                                                 unlist(MC_uni_ranmean_var[2001:3000]),
                                                 unlist(MC_uni_ranmean_kurt[2001:3000]),
                                                 unlist(MC_uni_ranmean_CPC[2001:3000])),
                                        nrow = 1000))
colnames(MC_uni_ranmean3) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

MC_uni_ranmean4 <- as.data.frame(matrix(data = c(unlist(MC_uni_means[3001:4000]),
                                                 unlist(MC_uni_ranmean_diff[3001:4000]),
                                                 unlist(MC_uni_ranmean_var[3001:4000]),
                                                 unlist(MC_uni_ranmean_kurt[3001:4000]),
                                                 unlist(MC_uni_ranmean_CPC[3001:4000])),
                                        nrow = 1000))
colnames(MC_uni_ranmean4) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

#merge data sets together for plotting

MC_uni_ranvar_all <- rbind(MC_uni_ranvar1, MC_uni_ranvar2, MC_uni_ranvar3,
                           MC_uni_ranvar4)
MC_uni_ranmean_all <- rbind(MC_uni_ranmean1, MC_uni_ranmean2, MC_uni_ranmean3,
                            MC_uni_ranmean4)

#add means and standard deviation labels and scale measures

ranvar_labels <- c(rep(c("means = (-5, 5)", "means = (-4, 4)", "means = (-3, 3)",
                         "means = (-2, 2)"), each = 1000))
ranmean_labels <- c(rep(c("sd = 0.5", "sd = 1", "sd = 1.5", "sd = 2"), each = 1000))

MC_uni_ranvar_all <- cbind(MC_uni_ranvar_all, label = ranvar_labels)
MC_uni_ranmean_all <- cbind(MC_uni_ranmean_all, label = ranmean_labels)

MC_uni_ranvar_all <- mutate(MC_uni_ranvar_all,
                            difference = scales::rescale(difference, to = c(0, 1)),
                            variance = scales::rescale(variance, to = c(0, 1)),
                            kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                            CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

MC_uni_ranmean_all <- mutate(MC_uni_ranmean_all,
                             difference = scales::rescale(difference, to = c(0, 1)),
                             variance = scales::rescale(variance, to = c(0, 1)),
                             kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                             CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

#plot simulation results (Figure S3)

ggplot(MC_uni_ranvar_all, aes(x = ranvar, y = value, linetype = Measure,
                              color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Standard Deviations", y = "Estimated Level of Polarization")

ggplot(MC_uni_ranmean_all, aes(x = ranmean, y = value, linetype = Measure,
                               color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Means", y = "Estimated Level of Polarization")

#calculate polarization ranks for random variance distributions

MC_uni_ranvar1 <- mutate(MC_uni_ranvar1, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-5, 5)", 1000))
MC_uni_ranvar2 <- mutate(MC_uni_ranvar2, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-4, 4)", 1000))
MC_uni_ranvar3 <- mutate(MC_uni_ranvar3, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-3, 3)", 1000))
MC_uni_ranvar4 <- mutate(MC_uni_ranvar4, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-2, 2)", 1000))

#calculate polarization ranks for random mean distributions

MC_uni_ranmean1 <- mutate(MC_uni_ranmean1, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 0.5", 1000))
MC_uni_ranmean2 <- mutate(MC_uni_ranmean2, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 1", 1000))
MC_uni_ranmean3 <- mutate(MC_uni_ranmean3, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 1.5", 1000))
MC_uni_ranmean4 <- mutate(MC_uni_ranmean4, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 2", 1000))

#combine polarization rank data into single data set for plotting

MC_uni_ranvar_plot <- rbind(MC_uni_ranvar1, MC_uni_ranvar2, MC_uni_ranvar3,
                            MC_uni_ranvar4)
MC_uni_ranmean_plot <- rbind(MC_uni_ranmean1, MC_uni_ranmean2, MC_uni_ranmean3,
                             MC_uni_ranmean4)

#tidy data and create labels for plotting

MC_uni_ranvar_plot <- pivot_longer(MC_uni_ranvar_plot, cols = 2:5,
                                   names_to = "measure", values_to = "value")
MC_uni_ranmean_plot <- pivot_longer(MC_uni_ranmean_plot, cols = 2:5,
                                    names_to = "measure", values_to = "value")

MC_uni_ranvar_plot$mean_label <- factor(MC_uni_ranvar_plot$mean_label,
                                        levels = c("means = (-2, 2)",
                                                   "means = (-3, 3)",
                                                   "means = (-4, 4)",
                                                   "means = (-5, 5)"))
MC_uni_ranmean_plot$variance_label <- factor(MC_uni_ranmean_plot$variance_label,
                                             levels = c("sd = 2", "sd = 1.5",
                                                        "sd = 1", "sd = 0.5"))

MC_uni_ranvar_plot$measure <- factor(MC_uni_ranvar_plot$measure,
                                     levels = c("difference", "variance",
                                                "kurtosis", "CPC"))
MC_uni_ranmean_plot$measure <- factor(MC_uni_ranmean_plot$measure,
                                      levels = c("difference", "variance",
                                                 "kurtosis", "CPC"))

#plot results of univariate Monte Carlo simulations (Figure 3)

ggplot(MC_uni_ranvar_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ mean_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(MC_uni_ranmean_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ variance_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###################################################################################
#Monte Carlo simulations for bivariate distributions (2 clusters)
###################################################################################

#generate Gaussian mixture distributions

MC_bi_ranvar <- list()
MC_bi_ranmean <- list()
MC_bi_variances <- c()
MC_bi_means <- c()

for (i in means_values) {
  for (j in 1:1000) {
    variance <- runif(1, min = 0.5, max = 2)
    data <- rmvnormmix(1000, c(0.5, 0.5), 
                       matrix(data = c(i + 0.00001, -i - 0.00001, -i, i), nrow = 2,
                              ncol = 2, byrow = TRUE), 
                       matrix(data = variance, nrow = 2, ncol = 2, byrow = TRUE))
    MC_bi_ranvar[[length(MC_bi_ranvar) + 1]] <- data
    MC_bi_variances[[length(MC_bi_variances) + 1]] <- variance
  }
  base::print(i)
}

for (i in variance_values) {
  for (j in 1:1000) {
    mean <- runif(1, min = 2, max = 5)
    data <- rmvnormmix(1000, c(0.5, 0.5), 
                       matrix(data = c(mean + 0.00001, -mean - 0.00001, -mean,
                                       mean), nrow = 2, ncol = 2, byrow = TRUE), 
                       matrix(data = i, nrow = 2, ncol = 2, byrow = TRUE))
    MC_bi_ranmean[[length(MC_bi_ranmean)+1]] <- data
    MC_bi_means[[length(MC_bi_means)+1]] <- mean
  }
  base::print(i)
}

#calculate various polarization scores for simulated distributions

MC_bi_ranvar_diff <- pblapply(MC_bi_ranvar,
                              function(x) diff.bi(x = x, y = 2, k = 2), cl = cores)
MC_bi_ranvar_var <- pblapply(MC_bi_ranvar, Euclidean, cl = cores)
MC_bi_ranvar_mardia <- pblapply(MC_bi_ranvar,
                                function(x) mardia(x = x, plot = FALSE), cl = cores)
MC_bi_ranvar_CPC <- pblapply(MC_bi_ranvar,
                             function(x) CPC(x, 2, "kmeans", adjust = TRUE,
                                             nstart = 30), cl = cores)

MC_bi_ranmean_diff <- pblapply(MC_bi_ranmean,
                               function(x) diff.bi(x = x, y = 2, k = 2), cl = cores)
MC_bi_ranmean_var <- pblapply(MC_bi_ranmean, Euclidean, cl = cores)
MC_bi_ranmean_mardia <- pblapply(MC_bi_ranmean,
                                 function(x) mardia(x = x, plot = FALSE),
                                 cl = cores)
MC_bi_ranmean_CPC <- pblapply(MC_bi_ranmean,
                              function(x) CPC(x, 2, "kmeans", adjust = TRUE,
                                              nstart = 30), cl = cores)

#extract kurtosis element from mardia object

MC_bi_ranvar_kurt <- list()
MC_bi_ranmean_kurt <- list()

for (i in 1:4000) {
  MC_bi_ranvar_kurt[[length(MC_bi_ranvar_kurt) + 1]] <- MC_bi_ranvar_mardia[[i]]$b2p
  MC_bi_ranmean_kurt[[length(MC_bi_ranmean_kurt) + 1]] <- MC_bi_ranmean_mardia[[i]]$b2p
}

#combine random variance polarization results into data frames, holding means constant within data frame for accurate rank calculations

MC_bi_ranvar1 <- as.data.frame(matrix(data = c(unlist(MC_bi_variances[1:1000]),
                                               unlist(MC_bi_ranvar_diff[1:1000]),
                                               unlist(MC_bi_ranvar_var[1:1000]),
                                               unlist(MC_bi_ranvar_kurt[1:1000]),
                                               unlist(MC_bi_ranvar_CPC[1:1000])),
                                      nrow = 1000))
colnames(MC_bi_ranvar1) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC_bi_ranvar2 <- as.data.frame(matrix(data = c(unlist(MC_bi_variances[1001:2000]),
                                               unlist(MC_bi_ranvar_diff[1001:2000]),
                                               unlist(MC_bi_ranvar_var[1001:2000]),
                                               unlist(MC_bi_ranvar_kurt[1001:2000]),
                                               unlist(MC_bi_ranvar_CPC[1001:2000])),
                                      nrow = 1000))
colnames(MC_bi_ranvar2) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC_bi_ranvar3 <- as.data.frame(matrix(data = c(unlist(MC_bi_variances[2001:3000]),
                                               unlist(MC_bi_ranvar_diff[2001:3000]),
                                               unlist(MC_bi_ranvar_var[2001:3000]),
                                               unlist(MC_bi_ranvar_kurt[2001:3000]),
                                               unlist(MC_bi_ranvar_CPC[2001:3000])),
                                      nrow = 1000))
colnames(MC_bi_ranvar3) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC_bi_ranvar4 <- as.data.frame(matrix(data = c(unlist(MC_bi_variances[3001:4000]),
                                               unlist(MC_bi_ranvar_diff[3001:4000]),
                                               unlist(MC_bi_ranvar_var[3001:4000]),
                                               unlist(MC_bi_ranvar_kurt[3001:4000]),
                                               unlist(MC_bi_ranvar_CPC[3001:4000])),
                                      nrow = 1000))
colnames(MC_bi_ranvar4) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

#combine random mean polarization results into data frames, holding variance constant within data frame for accurate rank calculations

MC_bi_ranmean1 <- as.data.frame(matrix(data = c(unlist(MC_bi_means[1:1000]),
                                                unlist(MC_bi_ranmean_diff[1:1000]),
                                                unlist(MC_bi_ranmean_var[1:1000]),
                                                unlist(MC_bi_ranmean_kurt[1:1000]),
                                                unlist(MC_bi_ranmean_CPC[1:1000])),
                                       nrow = 1000))
colnames(MC_bi_ranmean1) <- c("ranmean", "difference", "variance", "kurtosis",
                              "CPC")

MC_bi_ranmean2 <- as.data.frame(matrix(data = c(unlist(MC_bi_means[1001:2000]),
                                                unlist(MC_bi_ranmean_diff[1001:2000]),
                                                unlist(MC_bi_ranmean_var[1001:2000]),
                                                unlist(MC_bi_ranmean_kurt[1001:2000]),
                                                unlist(MC_bi_ranmean_CPC[1001:2000])),
                                       nrow = 1000))
colnames(MC_bi_ranmean2) <- c("ranmean", "difference", "variance", "kurtosis",
                              "CPC")

MC_bi_ranmean3 <- as.data.frame(matrix(data = c(unlist(MC_bi_means[2001:3000]),
                                                unlist(MC_bi_ranmean_diff[2001:3000]),
                                                unlist(MC_bi_ranmean_var[2001:3000]),
                                                unlist(MC_bi_ranmean_kurt[2001:3000]),
                                                unlist(MC_bi_ranmean_CPC[2001:3000])),
                                       nrow = 1000))
colnames(MC_bi_ranmean3) <- c("ranmean", "difference", "variance", "kurtosis",
                              "CPC")

MC_bi_ranmean4 <- as.data.frame(matrix(data = c(unlist(MC_bi_means[3001:4000]),
                                                unlist(MC_bi_ranmean_diff[3001:4000]),
                                                unlist(MC_bi_ranmean_var[3001:4000]),
                                                unlist(MC_bi_ranmean_kurt[3001:4000]),
                                                unlist(MC_bi_ranmean_CPC[3001:4000])),
                                       nrow = 1000))
colnames(MC_bi_ranmean4) <- c("ranmean", "difference", "variance", "kurtosis",
                              "CPC")

#merge data sets together for plotting

MC_bi_ranvar_all <- rbind(MC_bi_ranvar1, MC_bi_ranvar2, MC_bi_ranvar3,
                          MC_bi_ranvar4)
MC_bi_ranmean_all <- rbind(MC_bi_ranmean1, MC_bi_ranmean2, MC_bi_ranmean3,
                           MC_bi_ranmean4)

#add means and standard deviation labels and scale measures

ranvar_labels <- c(rep(c("means = (-5, 5)", "means = (-4, 4)", "means = (-3, 3)",
                         "means = (-2, 2)"), each = 1000))
ranmean_labels <- c(rep(c("sd = 0.5", "sd = 1", "sd = 1.5", "sd = 2"), each = 1000))

MC_bi_ranvar_all <- cbind(MC_bi_ranvar_all, label = ranvar_labels)
MC_bi_ranmean_all <- cbind(MC_bi_ranmean_all, label = ranmean_labels)

MC_bi_ranvar_all <- mutate(MC_bi_ranvar_all,
                           difference = scales::rescale(difference, to = c(0, 1)),
                           variance = scales::rescale(variance, to = c(0, 1)),
                           kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                           CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

MC_bi_ranmean_all <- mutate(MC_bi_ranmean_all,
                            difference = scales::rescale(difference, to = c(0, 1)),
                            variance = scales::rescale(variance, to = c(0, 1)),
                            kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                            CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

#plot simulation results (Figure S4)

ggplot(MC_bi_ranvar_all, aes(x = ranvar, y = value, linetype = Measure,
                             color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Standard Deviations", y = "Estimated Level of Polarization")

ggplot(MC_bi_ranmean_all, aes(x = ranmean, y = value, linetype = Measure,
                              color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Means", y = "Estimated Level of Polarization")

#calculate polarization ranks for random variance distributions

MC_bi_ranvar1 <- mutate(MC_bi_ranvar1, rank_true = rank(-ranvar),
                        difference = rank(difference),
                        variance = rank(variance),
                        kurtosis = rank(-kurtosis),
                        CPC = rank(CPC),
                        mean_label = rep("means = (-5, 5)", 1000))
MC_bi_ranvar2 <- mutate(MC_bi_ranvar2, rank_true = rank(-ranvar),
                        difference = rank(difference),
                        variance = rank(variance),
                        kurtosis = rank(-kurtosis),
                        CPC = rank(CPC),
                        mean_label = rep("means = (-4, 4)", 1000))
MC_bi_ranvar3 <- mutate(MC_bi_ranvar3, rank_true = rank(-ranvar),
                        difference = rank(difference),
                        variance = rank(variance),
                        kurtosis = rank(-kurtosis),
                        CPC = rank(CPC),
                        mean_label = rep("means = (-3, 3)", 1000))
MC_bi_ranvar4 <- mutate(MC_bi_ranvar4, rank_true = rank(-ranvar),
                        difference = rank(difference),
                        variance = rank(variance),
                        kurtosis = rank(-kurtosis),
                        CPC = rank(CPC),
                        mean_label = rep("means = (-2, 2)", 1000))

#calculate polarization ranks for random mean distributions

MC_bi_ranmean1 <- mutate(MC_bi_ranmean1, rank_true = rank(ranmean),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         variance_label = rep("sd = 0.5", 1000))
MC_bi_ranmean2 <- mutate(MC_bi_ranmean2, rank_true = rank(ranmean),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         variance_label = rep("sd = 1", 1000))
MC_bi_ranmean3 <- mutate(MC_bi_ranmean3, rank_true = rank(ranmean),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         variance_label = rep("sd = 1.5", 1000))
MC_bi_ranmean4 <- mutate(MC_bi_ranmean4, rank_true = rank(ranmean),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         variance_label = rep("sd = 2", 1000))

#combine polarization rank data into single data set for plotting

MC_bi_ranvar_plot <- rbind(MC_bi_ranvar1, MC_bi_ranvar2, MC_bi_ranvar3,
                           MC_bi_ranvar4)
MC_bi_ranmean_plot <- rbind(MC_bi_ranmean1, MC_bi_ranmean2, MC_bi_ranmean3,
                            MC_bi_ranmean4)

#tidy data and create labels for plotting

MC_bi_ranvar_plot <- pivot_longer(MC_bi_ranvar_plot, cols = 2:5,
                                  names_to = "measure", values_to = "value")
MC_bi_ranmean_plot <- pivot_longer(MC_bi_ranmean_plot, cols = 2:5,
                                   names_to = "measure", values_to = "value")

MC_bi_ranvar_plot$mean_label <- factor(MC_bi_ranvar_plot$mean_label,
                                       levels = c("means = (-2, 2)",
                                                  "means = (-3, 3)",
                                                  "means = (-4, 4)",
                                                  "means = (-5, 5)"))
MC_bi_ranmean_plot$variance_label <- factor(MC_bi_ranmean_plot$variance_label,
                                            levels = c("sd = 2", "sd = 1.5",
                                                       "sd = 1", "sd = 0.5"))
MC_bi_ranvar_plot$measure <- factor(MC_bi_ranvar_plot$measure,
                                    levels = c("difference", "variance", "kurtosis",
                                               "CPC"))
MC_bi_ranmean_plot$measure <- factor(MC_bi_ranmean_plot$measure,
                                     levels = c("difference", "variance",
                                                "kurtosis", "CPC"))

#plot results of bivariate Monte Carlo simulations (Figure 4)

ggplot(MC_bi_ranvar_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ mean_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(MC_bi_ranmean_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ variance_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###################################################################################
#Monte Carlo simulations for univariate distributions (3 clusters)
###################################################################################

#set seed

set.seed(123)

#generate Gaussian mixture distributions

MC3_uni_ranvar <- list()
MC3_uni_ranmean <- list()
MC3_uni_variances <- c()
MC3_uni_means <- c()

for (i in means_values) {
  for (j in 1:1000) {
    variance <- runif(1, min = 0.5, max = 2)
    data <- rnormmix(1000, c(rep(0.33, 3)), c(-i, 0, i), c(rep(variance, 3)))
    MC3_uni_ranvar[[length(MC3_uni_ranvar) + 1]] <- data
    MC3_uni_variances[[length(MC3_uni_variances) + 1]] <- variance
  }
  base::print(i)
}

for (i in variance_values) {
  for (j in 1:1000) {
    mean <- runif(1, min = 2, max = 5)
    data <- rnormmix(1000, c(rep(0.33, 3)), c(-mean, 0, mean), c(rep(i, 3)))
    MC3_uni_ranmean[[length(MC3_uni_ranmean) + 1]] <- data
    MC3_uni_means[[length(MC3_uni_means) + 1]] <- mean
  }
  base::print(i)
}

#calculate various polarization scores for simulated distributions

MC3_uni_ranvar_diff <- pblapply(MC3_uni_ranvar,
                                function(x) diff.uni(x = x, y = 1, k = 3),
                                cl = cores)
MC3_uni_ranvar_var <- pblapply(MC3_uni_ranvar, var, cl = cores)
MC3_uni_ranvar_kurt <- pblapply(MC3_uni_ranvar, kurtosi, cl = cores)
MC3_uni_ranvar_CPC <- pblapply(MC3_uni_ranvar,
                               function(x) CPC(x, 3, "kmeans", adjust = TRUE,
                                               nstart = 30), cl = cores)

MC3_uni_ranmean_diff <- pblapply(MC3_uni_ranmean,
                                 function(x) diff.uni(x = x, y = 1, k = 3),
                                 cl = cores)
MC3_uni_ranmean_var <- pblapply(MC3_uni_ranmean, var, cl = cores)
MC3_uni_ranmean_kurt <- pblapply(MC3_uni_ranmean, kurtosi, cl = cores)
MC3_uni_ranmean_CPC <- pblapply(MC3_uni_ranmean,
                                function(x) CPC(x, 3, "kmeans", adjust = TRUE,
                                                nstart = 30), cl = cores)

#combine random variance polarization results into data frames, holding means constant within data frame for accurate rank calculations

MC3_uni_ranvar1 <- as.data.frame(matrix(data = c(unlist(MC3_uni_variances[1:1000]),
                                                 unlist(MC3_uni_ranvar_diff[1:1000]),
                                                 unlist(MC3_uni_ranvar_var[1:1000]),
                                                 unlist(MC3_uni_ranvar_kurt[1:1000]),
                                                 unlist(MC3_uni_ranvar_CPC[1:1000])),
                                        nrow = 1000))
colnames(MC3_uni_ranvar1) <- c("ranvar", "difference", "variance", "kurtosis",
                               "CPC")

MC3_uni_ranvar2 <- as.data.frame(matrix(data = c(unlist(MC3_uni_variances[1001:2000]),
                                                 unlist(MC3_uni_ranvar_diff[1001:2000]),
                                                 unlist(MC3_uni_ranvar_var[1001:2000]),
                                                 unlist(MC3_uni_ranvar_kurt[1001:2000]),
                                                 unlist(MC3_uni_ranvar_CPC[1001:2000])),
                                        nrow = 1000))
colnames(MC3_uni_ranvar2) <- c("ranvar", "difference", "variance", "kurtosis",
                               "CPC")

MC3_uni_ranvar3 <- as.data.frame(matrix(data = c(unlist(MC3_uni_variances[2001:3000]),
                                                 unlist(MC3_uni_ranvar_diff[2001:3000]),
                                                 unlist(MC3_uni_ranvar_var[2001:3000]),
                                                 unlist(MC3_uni_ranvar_kurt[2001:3000]),
                                                 unlist(MC3_uni_ranvar_CPC[2001:3000])),
                                        nrow = 1000))
colnames(MC3_uni_ranvar3) <- c("ranvar", "difference", "variance", "kurtosis",
                               "CPC")

MC3_uni_ranvar4 <- as.data.frame(matrix(data = c(unlist(MC3_uni_variances[3001:4000]),
                                                 unlist(MC3_uni_ranvar_diff[3001:4000]),
                                                 unlist(MC3_uni_ranvar_var[3001:4000]),
                                                 unlist(MC3_uni_ranvar_kurt[3001:4000]),
                                                 unlist(MC3_uni_ranvar_CPC[3001:4000])),
                                        nrow = 1000))
colnames(MC3_uni_ranvar4) <- c("ranvar", "difference", "variance", "kurtosis",
                               "CPC")

#combine random mean polarization results into data frames, holding variance constant within data frame for accurate rank calculations

MC3_uni_ranmean1 <- as.data.frame(matrix(data = c(unlist(MC3_uni_means[1:1000]),
                                                  unlist(MC3_uni_ranmean_diff[1:1000]),
                                                  unlist(MC3_uni_ranmean_var[1:1000]),
                                                  unlist(MC3_uni_ranmean_kurt[1:1000]),
                                                  unlist(MC3_uni_ranmean_CPC[1:1000])),
                                         nrow = 1000))
colnames(MC3_uni_ranmean1) <- c("ranmean", "difference", "variance", "kurtosis",
                                "CPC")

MC3_uni_ranmean2 <- as.data.frame(matrix(data = c(unlist(MC3_uni_means[1001:2000]),
                                                  unlist(MC3_uni_ranmean_diff[1001:2000]),
                                                  unlist(MC3_uni_ranmean_var[1001:2000]),
                                                  unlist(MC3_uni_ranmean_kurt[1001:2000]),
                                                  unlist(MC3_uni_ranmean_CPC[1001:2000])),
                                         nrow = 1000))
colnames(MC3_uni_ranmean2) <- c("ranmean", "difference", "variance", "kurtosis",
                                "CPC")

MC3_uni_ranmean3 <- as.data.frame(matrix(data = c(unlist(MC3_uni_means[2001:3000]),
                                                  unlist(MC3_uni_ranmean_diff[2001:3000]),
                                                  unlist(MC3_uni_ranmean_var[2001:3000]),
                                                  unlist(MC3_uni_ranmean_kurt[2001:3000]),
                                                  unlist(MC3_uni_ranmean_CPC[2001:3000])),
                                         nrow = 1000))
colnames(MC3_uni_ranmean3) <- c("ranmean", "difference", "variance", "kurtosis",
                                "CPC")

MC3_uni_ranmean4 <- as.data.frame(matrix(data = c(unlist(MC3_uni_means[3001:4000]),
                                                  unlist(MC3_uni_ranmean_diff[3001:4000]),
                                                  unlist(MC3_uni_ranmean_var[3001:4000]),
                                                  unlist(MC3_uni_ranmean_kurt[3001:4000]),
                                                  unlist(MC3_uni_ranmean_CPC[3001:4000])),
                                         nrow = 1000))
colnames(MC3_uni_ranmean4) <- c("ranmean", "difference", "variance", "kurtosis",
                                "CPC")

#merge data sets together for plotting

MC3_uni_ranvar_all <- rbind(MC3_uni_ranvar1, MC3_uni_ranvar2, MC3_uni_ranvar3,
                            MC3_uni_ranvar4)
MC3_uni_ranmean_all <- rbind(MC3_uni_ranmean1, MC3_uni_ranmean2, MC3_uni_ranmean3,
                             MC3_uni_ranmean4)

#add means and standard deviation labels and scale measures

ranvar_labels <- c(rep(c("means = (-5, 5)", "means = (-4, 4)", "means = (-3, 3)",
                         "means = (-2, 2)"), each = 1000))
ranmean_labels <- c(rep(c("sd = 0.5", "sd = 1", "sd = 1.5", "sd = 2"), each = 1000))

MC3_uni_ranvar_all <- cbind(MC3_uni_ranvar_all, label = ranvar_labels)
MC3_uni_ranmean_all <- cbind(MC3_uni_ranmean_all, label = ranmean_labels)

MC3_uni_ranvar_all <- mutate(MC3_uni_ranvar_all,
                             difference = scales::rescale(difference, to = c(0, 1)),
                             variance = scales::rescale(variance, to = c(0, 1)),
                             kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                             CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

MC3_uni_ranmean_all <- mutate(MC3_uni_ranmean_all,
                              difference = scales::rescale(difference,
                                                           to = c(0, 1)),
                              variance = scales::rescale(variance, to = c(0, 1)),
                              kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                              CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

#plot simulation results (Figure S5)

ggplot(MC3_uni_ranvar_all, aes(x = ranvar, y = value, linetype = Measure,
                               color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Standard Deviations", y = "Estimated Level of Polarization")

ggplot(MC3_uni_ranmean_all, aes(x = ranmean, y = value, linetype = Measure,
                                color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Means", y = "Estimated Level of Polarization")

#calculate polarization ranks for random variance distributions

MC3_uni_ranvar1 <- mutate(MC3_uni_ranvar1, rank_true = rank(-ranvar),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          mean_label = rep("means = (-5, 5)", 1000))
MC3_uni_ranvar2 <- mutate(MC3_uni_ranvar2, rank_true = rank(-ranvar),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          mean_label = rep("means = (-4, 4)", 1000))
MC3_uni_ranvar3 <- mutate(MC3_uni_ranvar3, rank_true = rank(-ranvar),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          mean_label = rep("means = (-3, 3)", 1000))
MC3_uni_ranvar4 <- mutate(MC3_uni_ranvar4, rank_true = rank(-ranvar),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          mean_label = rep("means = (-2, 2)", 1000))

#calculate polarization ranks for random mean distributions

MC3_uni_ranmean1 <- mutate(MC3_uni_ranmean1, rank_true = rank(ranmean),
                           difference = rank(difference),
                           variance = rank(variance),
                           kurtosis = rank(-kurtosis),
                           CPC = rank(CPC),
                           variance_label = rep("sd = 0.5", 1000))
MC3_uni_ranmean2 <- mutate(MC3_uni_ranmean2, rank_true = rank(ranmean),
                           difference = rank(difference),
                           variance = rank(variance),
                           kurtosis = rank(-kurtosis),
                           CPC = rank(CPC),
                           variance_label = rep("sd = 1", 1000))
MC3_uni_ranmean3 <- mutate(MC3_uni_ranmean3, rank_true = rank(ranmean),
                           difference = rank(difference),
                           variance = rank(variance),
                           kurtosis = rank(-kurtosis),
                           CPC = rank(CPC),
                           variance_label = rep("sd = 1.5", 1000))
MC3_uni_ranmean4 <- mutate(MC3_uni_ranmean4, rank_true = rank(ranmean),
                           difference = rank(difference),
                           variance = rank(variance),
                           kurtosis = rank(-kurtosis),
                           CPC = rank(CPC),
                           variance_label = rep("sd = 2", 1000))

#combine polarization rank data into single data set for plotting

MC3_uni_ranvar_plot <- rbind(MC3_uni_ranvar1, MC3_uni_ranvar2, MC3_uni_ranvar3,
                             MC3_uni_ranvar4)
MC3_uni_ranmean_plot <- rbind(MC3_uni_ranmean1, MC3_uni_ranmean2, MC3_uni_ranmean3,
                              MC3_uni_ranmean4)

#tidy data and create labels for plotting

MC3_uni_ranvar_plot <- pivot_longer(MC3_uni_ranvar_plot, cols = 2:5,
                                    names_to = "measure", values_to = "value")
MC3_uni_ranmean_plot <- pivot_longer(MC3_uni_ranmean_plot, cols = 2:5,
                                     names_to = "measure", values_to = "value")

MC3_uni_ranvar_plot$mean_label <- factor(MC3_uni_ranvar_plot$mean_label,
                                         levels = c("means = (-2, 2)",
                                                    "means = (-3, 3)",
                                                    "means = (-4, 4)",
                                                    "means = (-5, 5)"))
MC3_uni_ranmean_plot$variance_label <- factor(MC3_uni_ranmean_plot$variance_label,
                                              levels = c("sd = 2", "sd = 1.5",
                                                         "sd = 1", "sd = 0.5"))

MC3_uni_ranvar_plot$measure <- factor(MC3_uni_ranvar_plot$measure,
                                      levels = c("difference", "variance",
                                                 "kurtosis", "CPC"))
MC3_uni_ranmean_plot$measure <- factor(MC3_uni_ranmean_plot$measure,
                                       levels = c("difference", "variance",
                                                  "kurtosis", "CPC"))

#plot results of univariate Monte Carlo simulations (Figure S7)

ggplot(MC3_uni_ranvar_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ mean_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(MC3_uni_ranmean_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ variance_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###################################################################################
#Monte Carlo simulations for bivariate distributions (3 clusters)
###################################################################################

#generate Gaussian mixture distributions

MC3_bi_ranvar <- list()
MC3_bi_ranmean <- list()
MC3_bi_variances <- c()
MC3_bi_means <- c()

for (i in means_values) {
  for (j in 1:1000) {
    variance <- runif(1, min = 0.5, max = 2)
    data <- rmvnormmix(1000, c(rep(0.33, 3)), 
                       matrix(data = c(i + 0.00001, 0, -i - 0.00001, -i, 0, i),
                              nrow = 3, ncol = 2, byrow = TRUE), 
                       matrix(data = variance, nrow = 3, ncol = 2, byrow = TRUE))
    MC3_bi_ranvar[[length(MC3_bi_ranvar) + 1]] <- data
    MC3_bi_variances[[length(MC3_bi_variances) + 1]] <- variance
  }
  base::print(i)
}

for (i in variance_values) {
  for (j in 1:1000) {
    mean <- runif(1, min = 2, max = 5)
    data <- rmvnormmix(1000, c(rep(0.33, 3)), 
                       matrix(data = c(mean + 0.00001, 0, -mean - 0.00001, -mean, 0,
                                       mean), nrow = 3, ncol = 2, byrow = TRUE), 
                       matrix(data = i, nrow = 3, ncol = 2, byrow = TRUE))
    MC3_bi_ranmean[[length(MC3_bi_ranmean)+1]] <- data
    MC3_bi_means[[length(MC3_bi_means)+1]] <- mean
  }
  base::print(i)
}

#calculate various polarization scores for simulated distributions

MC3_bi_ranvar_diff <- pblapply(MC3_bi_ranvar,
                               function(x) diff.bi(x = x, y = 2, k = 3), cl = cores)
MC3_bi_ranvar_var <- pblapply(MC3_bi_ranvar, Euclidean, cl = cores)
MC3_bi_ranvar_mardia <- pblapply(MC3_bi_ranvar,
                                 function(x) mardia(x = x, plot = FALSE),
                                 cl = cores)
MC3_bi_ranvar_CPC <- pblapply(MC3_bi_ranvar,
                              function(x) CPC(x, 3, "kmeans", adjust = TRUE,
                                              nstart = 30), cl = cores)

MC3_bi_ranmean_diff <- pblapply(MC3_bi_ranmean,
                                function(x) diff.bi(x = x, y = 2, k = 3),
                                cl = cores)
MC3_bi_ranmean_var <- pblapply(MC3_bi_ranmean, Euclidean, cl = cores)
MC3_bi_ranmean_mardia <- pblapply(MC3_bi_ranmean,
                                  function(x) mardia(x = x, plot = FALSE),
                                  cl = cores)
MC3_bi_ranmean_CPC <- pblapply(MC3_bi_ranmean,
                               function(x) CPC(x, 3, "kmeans", adjust = TRUE,
                                               nstart = 30), cl = cores)

#extract kurtosis element from mardia object

MC3_bi_ranvar_kurt <- list()
MC3_bi_ranmean_kurt <- list()

for (i in 1:4000) {
  MC3_bi_ranvar_kurt[[length(MC3_bi_ranvar_kurt) + 1]] <- MC3_bi_ranvar_mardia[[i]]$b2p
  MC3_bi_ranmean_kurt[[length(MC3_bi_ranmean_kurt) + 1]] <- MC3_bi_ranmean_mardia[[i]]$b2p
}

#combine random variance polarization results into data frames, holding means constant within data frame for accurate rank calculations

MC3_bi_ranvar1 <- as.data.frame(matrix(data = c(unlist(MC3_bi_variances[1:1000]),
                                                unlist(MC3_bi_ranvar_diff[1:1000]),
                                                unlist(MC3_bi_ranvar_var[1:1000]),
                                                unlist(MC3_bi_ranvar_kurt[1:1000]),
                                                unlist(MC3_bi_ranvar_CPC[1:1000])),
                                       nrow = 1000))
colnames(MC3_bi_ranvar1) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC3_bi_ranvar2 <- as.data.frame(matrix(data = c(unlist(MC3_bi_variances[1001:2000]),
                                                unlist(MC3_bi_ranvar_diff[1001:2000]),
                                                unlist(MC3_bi_ranvar_var[1001:2000]),
                                                unlist(MC3_bi_ranvar_kurt[1001:2000]),
                                                unlist(MC3_bi_ranvar_CPC[1001:2000])),
                                       nrow = 1000))
colnames(MC3_bi_ranvar2) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC3_bi_ranvar3 <- as.data.frame(matrix(data = c(unlist(MC3_bi_variances[2001:3000]),
                                                unlist(MC3_bi_ranvar_diff[2001:3000]),
                                                unlist(MC3_bi_ranvar_var[2001:3000]),
                                                unlist(MC3_bi_ranvar_kurt[2001:3000]),
                                                unlist(MC3_bi_ranvar_CPC[2001:3000])),
                                       nrow = 1000))
colnames(MC3_bi_ranvar3) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC3_bi_ranvar4 <- as.data.frame(matrix(data = c(unlist(MC3_bi_variances[3001:4000]),
                                                unlist(MC3_bi_ranvar_diff[3001:4000]),
                                                unlist(MC3_bi_ranvar_var[3001:4000]),
                                                unlist(MC3_bi_ranvar_kurt[3001:4000]),
                                                unlist(MC3_bi_ranvar_CPC[3001:4000])),
                                       nrow = 1000))
colnames(MC3_bi_ranvar4) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

#combine random mean polarization results into data frames, holding variance constant within data frame for accurate rank calculations

MC3_bi_ranmean1 <- as.data.frame(matrix(data = c(unlist(MC3_bi_means[1:1000]),
                                                 unlist(MC3_bi_ranmean_diff[1:1000]),
                                                 unlist(MC3_bi_ranmean_var[1:1000]),
                                                 unlist(MC3_bi_ranmean_kurt[1:1000]),
                                                 unlist(MC3_bi_ranmean_CPC[1:1000])),
                                        nrow = 1000))
colnames(MC3_bi_ranmean1) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

MC3_bi_ranmean2 <- as.data.frame(matrix(data = c(unlist(MC3_bi_means[1001:2000]),
                                                 unlist(MC3_bi_ranmean_diff[1001:2000]),
                                                 unlist(MC3_bi_ranmean_var[1001:2000]),
                                                 unlist(MC3_bi_ranmean_kurt[1001:2000]),
                                                 unlist(MC3_bi_ranmean_CPC[1001:2000])),
                                        nrow = 1000))
colnames(MC3_bi_ranmean2) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

MC3_bi_ranmean3 <- as.data.frame(matrix(data = c(unlist(MC3_bi_means[2001:3000]),
                                                 unlist(MC3_bi_ranmean_diff[2001:3000]),
                                                 unlist(MC3_bi_ranmean_var[2001:3000]),
                                                 unlist(MC3_bi_ranmean_kurt[2001:3000]),
                                                 unlist(MC3_bi_ranmean_CPC[2001:3000])),
                                        nrow = 1000))
colnames(MC3_bi_ranmean3) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

MC3_bi_ranmean4 <- as.data.frame(matrix(data = c(unlist(MC3_bi_means[3001:4000]),
                                                 unlist(MC3_bi_ranmean_diff[3001:4000]),
                                                 unlist(MC3_bi_ranmean_var[3001:4000]),
                                                 unlist(MC3_bi_ranmean_kurt[3001:4000]),
                                                 unlist(MC3_bi_ranmean_CPC[3001:4000])),
                                        nrow = 1000))
colnames(MC3_bi_ranmean4) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

#merge data sets together for plotting

MC3_bi_ranvar_all <- rbind(MC3_bi_ranvar1, MC3_bi_ranvar2, MC3_bi_ranvar3,
                           MC3_bi_ranvar4)
MC3_bi_ranmean_all <- rbind(MC3_bi_ranmean1, MC3_bi_ranmean2, MC3_bi_ranmean3,
                            MC3_bi_ranmean4)

#add means and standard deviation labels and scale measures

ranvar_labels <- c(rep(c("means = (-5, 5)", "means = (-4, 4)", "means = (-3, 3)",
                         "means = (-2, 2)"), each = 1000))
ranmean_labels <- c(rep(c("sd = 0.5", "sd = 1", "sd = 1.5", "sd = 2"), each = 1000))

MC3_bi_ranvar_all <- cbind(MC3_bi_ranvar_all, label = ranvar_labels)
MC3_bi_ranmean_all <- cbind(MC3_bi_ranmean_all, label = ranmean_labels)

MC3_bi_ranvar_all <- mutate(MC3_bi_ranvar_all,
                            difference = scales::rescale(difference, to = c(0, 1)),
                            variance = scales::rescale(variance, to = c(0, 1)),
                            kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                            CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

MC3_bi_ranmean_all <- mutate(MC3_bi_ranmean_all,
                             difference = scales::rescale(difference,
                                                          to = c(0, 1)),
                             variance = scales::rescale(variance, to = c(0, 1)),
                             kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                             CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

#plot simulation results (Figure S6)

ggplot(MC3_bi_ranvar_all, aes(x = ranvar, y = value, linetype = Measure,
                              color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Standard Deviations", y = "Estimated Level of Polarization")

ggplot(MC3_bi_ranmean_all, aes(x = ranmean, y = value, linetype = Measure,
                               color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Means", y = "Estimated Level of Polarization")

#calculate polarization ranks for random variance distributions

MC3_bi_ranvar1 <- mutate(MC3_bi_ranvar1, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-5, 5)", 1000))
MC3_bi_ranvar2 <- mutate(MC3_bi_ranvar2, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-4, 4)", 1000))
MC3_bi_ranvar3 <- mutate(MC3_bi_ranvar3, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-3, 3)", 1000))
MC3_bi_ranvar4 <- mutate(MC3_bi_ranvar4, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-2, 2)", 1000))

#calculate polarization ranks for random mean distributions

MC3_bi_ranmean1 <- mutate(MC3_bi_ranmean1, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 0.5", 1000))
MC3_bi_ranmean2 <- mutate(MC3_bi_ranmean2, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 1", 1000))
MC3_bi_ranmean3 <- mutate(MC3_bi_ranmean3, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 1.5", 1000))
MC3_bi_ranmean4 <- mutate(MC3_bi_ranmean4, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 2", 1000))

#combine polarization rank data into single data set for plotting

MC3_bi_ranvar_plot <- rbind(MC3_bi_ranvar1, MC3_bi_ranvar2, MC3_bi_ranvar3,
                            MC3_bi_ranvar4)
MC3_bi_ranmean_plot <- rbind(MC3_bi_ranmean1, MC3_bi_ranmean2, MC3_bi_ranmean3,
                             MC3_bi_ranmean4)

#tidy data and create labels for plotting

MC3_bi_ranvar_plot <- pivot_longer(MC3_bi_ranvar_plot, cols = 2:5,
                                   names_to = "measure", values_to = "value")
MC3_bi_ranmean_plot <- pivot_longer(MC3_bi_ranmean_plot, cols = 2:5,
                                    names_to = "measure", values_to = "value")

MC3_bi_ranvar_plot$mean_label <- factor(MC3_bi_ranvar_plot$mean_label,
                                        levels = c("means = (-2, 2)",
                                                   "means = (-3, 3)",
                                                   "means = (-4, 4)",
                                                   "means = (-5, 5)"))
MC3_bi_ranmean_plot$variance_label <- factor(MC3_bi_ranmean_plot$variance_label,
                                             levels = c("sd = 2", "sd = 1.5",
                                                        "sd = 1", "sd = 0.5"))
MC3_bi_ranvar_plot$measure <- factor(MC3_bi_ranvar_plot$measure,
                                     levels = c("difference", "variance",
                                                "kurtosis", "CPC"))
MC3_bi_ranmean_plot$measure <- factor(MC3_bi_ranmean_plot$measure,
                                      levels = c("difference", "variance",
                                                 "kurtosis", "CPC"))

#plot results of bivariate Monte Carlo simulations (Figure S8)

ggplot(MC3_bi_ranvar_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ mean_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(MC3_bi_ranmean_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ variance_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###################################################################################
#Monte Carlo simulations for univariate distributions (4 clusters)
###################################################################################

#set seed

set.seed(123)

#generate Gaussian mixture distributions

MC4_uni_ranvar <- list()
MC4_uni_ranmean <- list()
MC4_uni_variances <- c()
MC4_uni_means <- c()

for (i in means_values) {
  for (j in 1:1000) {
    variance <- runif(1, min = 0.5, max = 2)
    data <- rnormmix(1000, c(rep(0.25, 4)), c(-i*2, -i/2, i/2, i*2),
                     c(rep(variance, 4)))
    MC4_uni_ranvar[[length(MC4_uni_ranvar) + 1]] <- data
    MC4_uni_variances[[length(MC4_uni_variances) + 1]] <- variance
  }
  base::print(i)
}

for (i in variance_values) {
  for (j in 1:1000) {
    mean <- runif(1, min = 2, max = 5)
    data <- rnormmix(1000, c(rep(0.25, 4)), c(-mean*2, -mean/2, mean/2, mean*2),
                     c(rep(i, 4)))
    MC4_uni_ranmean[[length(MC4_uni_ranmean) + 1]] <- data
    MC4_uni_means[[length(MC4_uni_means) + 1]] <- mean
  }
  base::print(i)
}

#calculate various polarization scores for simulated distributions

MC4_uni_ranvar_diff <- pblapply(MC4_uni_ranvar,
                                function(x) diff.uni(x = x, y = 1, k = 4),
                                cl = cores)
MC4_uni_ranvar_var <- pblapply(MC4_uni_ranvar, var, cl = cores)
MC4_uni_ranvar_kurt <- pblapply(MC4_uni_ranvar, kurtosi, cl = cores)
MC4_uni_ranvar_CPC <- pblapply(MC4_uni_ranvar,
                               function(x) CPC(x, 4, "kmeans", adjust = TRUE,
                                               nstart = 30), cl = cores)

MC4_uni_ranmean_diff <- pblapply(MC4_uni_ranmean,
                                 function(x) diff.uni(x = x, y = 1, k = 4),
                                 cl = cores)
MC4_uni_ranmean_var <- pblapply(MC4_uni_ranmean, var, cl = cores)
MC4_uni_ranmean_kurt <- pblapply(MC4_uni_ranmean, kurtosi, cl = cores)
MC4_uni_ranmean_CPC <- pblapply(MC4_uni_ranmean,
                                function(x) CPC(x, 4, "kmeans", adjust = TRUE,
                                                nstart = 30), cl = cores)

#combine random variance polarization results into data frames, holding means constant within data frame for accurate rank calculations

MC4_uni_ranvar1 <- as.data.frame(matrix(data = c(unlist(MC4_uni_variances[1:1000]),
                                                 unlist(MC4_uni_ranvar_diff[1:1000]),
                                                 unlist(MC4_uni_ranvar_var[1:1000]),
                                                 unlist(MC4_uni_ranvar_kurt[1:1000]),
                                                 unlist(MC4_uni_ranvar_CPC[1:1000])),
                                        nrow = 1000))
colnames(MC4_uni_ranvar1) <- c("ranvar", "difference", "variance", "kurtosis",
                               "CPC")

MC4_uni_ranvar2 <- as.data.frame(matrix(data = c(unlist(MC4_uni_variances[1001:2000]),
                                                 unlist(MC4_uni_ranvar_diff[1001:2000]),
                                                 unlist(MC4_uni_ranvar_var[1001:2000]),
                                                 unlist(MC4_uni_ranvar_kurt[1001:2000]),
                                                 unlist(MC4_uni_ranvar_CPC[1001:2000])),
                                        nrow = 1000))
colnames(MC4_uni_ranvar2) <- c("ranvar", "difference", "variance", "kurtosis",
                               "CPC")

MC4_uni_ranvar3 <- as.data.frame(matrix(data = c(unlist(MC4_uni_variances[2001:3000]),
                                                 unlist(MC4_uni_ranvar_diff[2001:3000]),
                                                 unlist(MC4_uni_ranvar_var[2001:3000]),
                                                 unlist(MC4_uni_ranvar_kurt[2001:3000]),
                                                 unlist(MC4_uni_ranvar_CPC[2001:3000])),
                                        nrow = 1000))
colnames(MC4_uni_ranvar3) <- c("ranvar", "difference", "variance", "kurtosis",
                               "CPC")

MC4_uni_ranvar4 <- as.data.frame(matrix(data = c(unlist(MC4_uni_variances[3001:4000]),
                                                 unlist(MC4_uni_ranvar_diff[3001:4000]),
                                                 unlist(MC4_uni_ranvar_var[3001:4000]),
                                                 unlist(MC4_uni_ranvar_kurt[3001:4000]),
                                                 unlist(MC4_uni_ranvar_CPC[3001:4000])),
                                        nrow = 1000))
colnames(MC4_uni_ranvar4) <- c("ranvar", "difference", "variance", "kurtosis",
                               "CPC")

#combine random mean polarization results into data frames, holding variance constant within data frame for accurate rank calculations

MC4_uni_ranmean1 <- as.data.frame(matrix(data = c(unlist(MC4_uni_means[1:1000]),
                                                  unlist(MC4_uni_ranmean_diff[1:1000]),
                                                  unlist(MC4_uni_ranmean_var[1:1000]),
                                                  unlist(MC4_uni_ranmean_kurt[1:1000]),
                                                  unlist(MC4_uni_ranmean_CPC[1:1000])),
                                         nrow = 1000))
colnames(MC4_uni_ranmean1) <- c("ranmean", "difference", "variance", "kurtosis",
                                "CPC")

MC4_uni_ranmean2 <- as.data.frame(matrix(data = c(unlist(MC4_uni_means[1001:2000]),
                                                  unlist(MC4_uni_ranmean_diff[1001:2000]),
                                                  unlist(MC4_uni_ranmean_var[1001:2000]),
                                                  unlist(MC4_uni_ranmean_kurt[1001:2000]),
                                                  unlist(MC4_uni_ranmean_CPC[1001:2000])),
                                         nrow = 1000))
colnames(MC4_uni_ranmean2) <- c("ranmean", "difference", "variance", "kurtosis",
                                "CPC")

MC4_uni_ranmean3 <- as.data.frame(matrix(data = c(unlist(MC4_uni_means[2001:3000]),
                                                  unlist(MC4_uni_ranmean_diff[2001:3000]),
                                                  unlist(MC4_uni_ranmean_var[2001:3000]),
                                                  unlist(MC4_uni_ranmean_kurt[2001:3000]),
                                                  unlist(MC4_uni_ranmean_CPC[2001:3000])),
                                         nrow = 1000))
colnames(MC4_uni_ranmean3) <- c("ranmean", "difference", "variance", "kurtosis",
                                "CPC")

MC4_uni_ranmean4 <- as.data.frame(matrix(data = c(unlist(MC4_uni_means[3001:4000]),
                                                  unlist(MC4_uni_ranmean_diff[3001:4000]),
                                                  unlist(MC4_uni_ranmean_var[3001:4000]),
                                                  unlist(MC4_uni_ranmean_kurt[3001:4000]),
                                                  unlist(MC4_uni_ranmean_CPC[3001:4000])),
                                         nrow = 1000))
colnames(MC4_uni_ranmean4) <- c("ranmean", "difference", "variance", "kurtosis",
                                "CPC")

#merge data sets together for plotting

MC4_uni_ranvar_all <- rbind(MC4_uni_ranvar1, MC4_uni_ranvar2, MC4_uni_ranvar3,
                            MC4_uni_ranvar4)
MC4_uni_ranmean_all <- rbind(MC4_uni_ranmean1, MC4_uni_ranmean2, MC4_uni_ranmean3,
                             MC4_uni_ranmean4)

#add means and standard deviation labels and scale measures

ranvar_labels <- c(rep(c("means = (-5, 5)", "means = (-4, 4)", "means = (-3, 3)",
                         "means = (-2, 2)"), each = 1000))
ranmean_labels <- c(rep(c("sd = 0.5", "sd = 1", "sd = 1.5", "sd = 2"), each = 1000))

MC4_uni_ranvar_all <- cbind(MC4_uni_ranvar_all, label = ranvar_labels)
MC4_uni_ranmean_all <- cbind(MC4_uni_ranmean_all, label = ranmean_labels)

MC4_uni_ranvar_all <- mutate(MC4_uni_ranvar_all,
                             difference = scales::rescale(difference, to = c(0, 1)),
                             variance = scales::rescale(variance, to = c(0, 1)),
                             kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                             CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

MC4_uni_ranmean_all <- mutate(MC4_uni_ranmean_all,
                              difference = scales::rescale(difference,
                                                           to = c(0, 1)),
                              variance = scales::rescale(variance, to = c(0, 1)),
                              kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                              CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

#plot simulation results (Figure S9)

ggplot(MC4_uni_ranvar_all, aes(x = ranvar, y = value, linetype = Measure,
                               color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Standard Deviations", y = "Estimated Level of Polarization")

ggplot(MC4_uni_ranmean_all, aes(x = ranmean, y = value, linetype = Measure,
                                color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Means", y = "Estimated Level of Polarization")

#calculate polarization ranks for random variance distributions

MC4_uni_ranvar1 <- mutate(MC4_uni_ranvar1, rank_true = rank(-ranvar),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          mean_label = rep("means = (-5, 5)", 1000))
MC4_uni_ranvar2 <- mutate(MC4_uni_ranvar2, rank_true = rank(-ranvar),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          mean_label = rep("means = (-4, 4)", 1000))
MC4_uni_ranvar3 <- mutate(MC4_uni_ranvar3, rank_true = rank(-ranvar),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          mean_label = rep("means = (-3, 3)", 1000))
MC4_uni_ranvar4 <- mutate(MC4_uni_ranvar4, rank_true = rank(-ranvar),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          mean_label = rep("means = (-2, 2)", 1000))

#calculate polarization ranks for random mean distributions

MC4_uni_ranmean1 <- mutate(MC4_uni_ranmean1, rank_true = rank(ranmean),
                           difference = rank(difference),
                           variance = rank(variance),
                           kurtosis = rank(-kurtosis),
                           CPC = rank(CPC),
                           variance_label = rep("sd = 0.5", 1000))
MC4_uni_ranmean2 <- mutate(MC4_uni_ranmean2, rank_true = rank(ranmean),
                           difference = rank(difference),
                           variance = rank(variance),
                           kurtosis = rank(-kurtosis),
                           CPC = rank(CPC),
                           variance_label = rep("sd = 1", 1000))
MC4_uni_ranmean3 <- mutate(MC4_uni_ranmean3, rank_true = rank(ranmean),
                           difference = rank(difference),
                           variance = rank(variance),
                           kurtosis = rank(-kurtosis),
                           CPC = rank(CPC),
                           variance_label = rep("sd = 1.5", 1000))
MC4_uni_ranmean4 <- mutate(MC4_uni_ranmean4, rank_true = rank(ranmean),
                           difference = rank(difference),
                           variance = rank(variance),
                           kurtosis = rank(-kurtosis),
                           CPC = rank(CPC),
                           variance_label = rep("sd = 2", 1000))

#combine polarization rank data into single data set for plotting

MC4_uni_ranvar_plot <- rbind(MC4_uni_ranvar1, MC4_uni_ranvar2, MC4_uni_ranvar3,
                             MC4_uni_ranvar4)
MC4_uni_ranmean_plot <- rbind(MC4_uni_ranmean1, MC4_uni_ranmean2, MC4_uni_ranmean3,
                              MC4_uni_ranmean4)

#tidy data and create labels for plotting

MC4_uni_ranvar_plot <- pivot_longer(MC4_uni_ranvar_plot, cols = 2:5,
                                    names_to = "measure", values_to = "value")
MC4_uni_ranmean_plot <- pivot_longer(MC4_uni_ranmean_plot, cols = 2:5,
                                     names_to = "measure", values_to = "value")

MC4_uni_ranvar_plot$mean_label <- factor(MC4_uni_ranvar_plot$mean_label,
                                         levels = c("means = (-2, 2)",
                                                    "means = (-3, 3)",
                                                    "means = (-4, 4)",
                                                    "means = (-5, 5)"))
MC4_uni_ranmean_plot$variance_label <- factor(MC4_uni_ranmean_plot$variance_label,
                                              levels = c("sd = 2", "sd = 1.5",
                                                         "sd = 1", "sd = 0.5"))

MC4_uni_ranvar_plot$measure <- factor(MC4_uni_ranvar_plot$measure,
                                      levels = c("difference", "variance",
                                                 "kurtosis", "CPC"))
MC4_uni_ranmean_plot$measure <- factor(MC4_uni_ranmean_plot$measure,
                                       levels = c("difference", "variance",
                                                  "kurtosis", "CPC"))

#plot results of univariate Monte Carlo simulations (Figure S11)

ggplot(MC4_uni_ranvar_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ mean_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(MC4_uni_ranmean_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ variance_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###################################################################################
#Monte Carlo simulations for bivariate distributions (4 clusters)
###################################################################################

#generate Gaussian mixture distributions

MC4_bi_ranvar <- list()
MC4_bi_ranmean <- list()
MC4_bi_variances <- c()
MC4_bi_means <- c()

for (i in means_values) {
  for (j in 1:1000) {
    variance <- runif(1, min = 0.5, max = 2)
    data <- rmvnormmix(1000, c(rep(0.25, 4)), 
                       matrix(data = c(i*2 + 0.00001, i/2 + 0.00001, -i/2 - 0.00001,
                                       -i*2 - 0.00001, -i*2, -i/2, i/2, i*2),
                              nrow = 4, ncol = 2, byrow = TRUE), 
                       matrix(data = variance, nrow = 4, ncol = 2, byrow = TRUE))
    MC4_bi_ranvar[[length(MC4_bi_ranvar) + 1]] <- data
    MC4_bi_variances[[length(MC4_bi_variances) + 1]] <- variance
  }
  base::print(i)
}

for (i in variance_values) {
  for (j in 1:1000) {
    mean <- runif(1, min = 2, max = 5)
    data <- rmvnormmix(1000, c(rep(0.25, 4)), 
                       matrix(data = c(mean*2 + 0.00001, mean/2 + 0.00001,
                                       -mean/2 - 0.00001, -mean*2 - 0.00001, -mean*2,
                                       -mean/2, mean/2, mean*2), nrow = 4, ncol = 2,
                              byrow = TRUE), 
                       matrix(data = i, nrow = 4, ncol = 2, byrow = TRUE))
    MC4_bi_ranmean[[length(MC4_bi_ranmean)+1]] <- data
    MC4_bi_means[[length(MC4_bi_means)+1]] <- mean
  }
  base::print(i)
}

#calculate various polarization scores for simulated distributions

MC4_bi_ranvar_diff <- pblapply(MC4_bi_ranvar,
                               function(x) diff.bi(x = x, y = 2, k = 4),
                               cl = cores)
MC4_bi_ranvar_var <- pblapply(MC4_bi_ranvar, Euclidean, cl = cores)
MC4_bi_ranvar_mardia <- pblapply(MC4_bi_ranvar,
                                 function(x) mardia(x = x, plot = FALSE),
                                 cl = cores)
MC4_bi_ranvar_CPC <- pblapply(MC4_bi_ranvar,
                              function(x) CPC(x, 4, "kmeans", adjust = TRUE,
                                              nstart = 30), cl = cores)

MC4_bi_ranmean_diff <- pblapply(MC4_bi_ranmean,
                                function(x) diff.bi(x = x, y = 2, k = 4),
                                cl = cores)
MC4_bi_ranmean_var <- pblapply(MC4_bi_ranmean, Euclidean, cl = cores)
MC4_bi_ranmean_mardia <- pblapply(MC4_bi_ranmean,
                                  function(x) mardia(x = x, plot = FALSE),
                                  cl = cores)
MC4_bi_ranmean_CPC <- pblapply(MC4_bi_ranmean,
                               function(x) CPC(x, 4, "kmeans", adjust = TRUE,
                                               nstart = 30), cl = cores)

#extract kurtosis element from mardia object

MC4_bi_ranvar_kurt <- list()
MC4_bi_ranmean_kurt <- list()

for (i in 1:4000) {
  MC4_bi_ranvar_kurt[[length(MC4_bi_ranvar_kurt) + 1]] <- MC4_bi_ranvar_mardia[[i]]$b2p
  MC4_bi_ranmean_kurt[[length(MC4_bi_ranmean_kurt) + 1]] <- MC4_bi_ranmean_mardia[[i]]$b2p
}

#combine random variance polarization results into data frames, holding means constant within data frame for accurate rank calculations

MC4_bi_ranvar1 <- as.data.frame(matrix(data = c(unlist(MC4_bi_variances[1:1000]),
                                                unlist(MC4_bi_ranvar_diff[1:1000]),
                                                unlist(MC4_bi_ranvar_var[1:1000]),
                                                unlist(MC4_bi_ranvar_kurt[1:1000]),
                                                unlist(MC4_bi_ranvar_CPC[1:1000])),
                                       nrow = 1000))
colnames(MC4_bi_ranvar1) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC4_bi_ranvar2 <- as.data.frame(matrix(data = c(unlist(MC4_bi_variances[1001:2000]),
                                                unlist(MC4_bi_ranvar_diff[1001:2000]),
                                                unlist(MC4_bi_ranvar_var[1001:2000]),
                                                unlist(MC4_bi_ranvar_kurt[1001:2000]),
                                                unlist(MC4_bi_ranvar_CPC[1001:2000])),
                                       nrow = 1000))
colnames(MC4_bi_ranvar2) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC4_bi_ranvar3 <- as.data.frame(matrix(data = c(unlist(MC4_bi_variances[2001:3000]),
                                                unlist(MC4_bi_ranvar_diff[2001:3000]),
                                                unlist(MC4_bi_ranvar_var[2001:3000]),
                                                unlist(MC4_bi_ranvar_kurt[2001:3000]),
                                                unlist(MC4_bi_ranvar_CPC[2001:3000])),
                                       nrow = 1000))
colnames(MC4_bi_ranvar3) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

MC4_bi_ranvar4 <- as.data.frame(matrix(data = c(unlist(MC4_bi_variances[3001:4000]),
                                                unlist(MC4_bi_ranvar_diff[3001:4000]),
                                                unlist(MC4_bi_ranvar_var[3001:4000]),
                                                unlist(MC4_bi_ranvar_kurt[3001:4000]),
                                                unlist(MC4_bi_ranvar_CPC[3001:4000])),
                                       nrow = 1000))
colnames(MC4_bi_ranvar4) <- c("ranvar", "difference", "variance", "kurtosis", "CPC")

#combine random mean polarization results into data frames, holding variance constant within data frame for accurate rank calculations

MC4_bi_ranmean1 <- as.data.frame(matrix(data = c(unlist(MC4_bi_means[1:1000]),
                                                 unlist(MC4_bi_ranmean_diff[1:1000]),
                                                 unlist(MC4_bi_ranmean_var[1:1000]),
                                                 unlist(MC4_bi_ranmean_kurt[1:1000]),
                                                 unlist(MC4_bi_ranmean_CPC[1:1000])),
                                        nrow = 1000))
colnames(MC4_bi_ranmean1) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

MC4_bi_ranmean2 <- as.data.frame(matrix(data = c(unlist(MC4_bi_means[1001:2000]),
                                                 unlist(MC4_bi_ranmean_diff[1001:2000]),
                                                 unlist(MC4_bi_ranmean_var[1001:2000]),
                                                 unlist(MC4_bi_ranmean_kurt[1001:2000]),
                                                 unlist(MC4_bi_ranmean_CPC[1001:2000])),
                                        nrow = 1000))
colnames(MC4_bi_ranmean2) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

MC4_bi_ranmean3 <- as.data.frame(matrix(data = c(unlist(MC4_bi_means[2001:3000]),
                                                 unlist(MC4_bi_ranmean_diff[2001:3000]),
                                                 unlist(MC4_bi_ranmean_var[2001:3000]),
                                                 unlist(MC4_bi_ranmean_kurt[2001:3000]),
                                                 unlist(MC4_bi_ranmean_CPC[2001:3000])),
                                        nrow = 1000))
colnames(MC4_bi_ranmean3) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

MC4_bi_ranmean4 <- as.data.frame(matrix(data = c(unlist(MC4_bi_means[3001:4000]),
                                                 unlist(MC4_bi_ranmean_diff[3001:4000]),
                                                 unlist(MC4_bi_ranmean_var[3001:4000]),
                                                 unlist(MC4_bi_ranmean_kurt[3001:4000]),
                                                 unlist(MC4_bi_ranmean_CPC[3001:4000])),
                                        nrow = 1000))
colnames(MC4_bi_ranmean4) <- c("ranmean", "difference", "variance", "kurtosis",
                               "CPC")

#merge data sets together for plotting

MC4_bi_ranvar_all <- rbind(MC4_bi_ranvar1, MC4_bi_ranvar2, MC4_bi_ranvar3,
                           MC4_bi_ranvar4)
MC4_bi_ranmean_all <- rbind(MC4_bi_ranmean1, MC4_bi_ranmean2, MC4_bi_ranmean3,
                            MC4_bi_ranmean4)

#add means and standard deviation labels and scale measures

ranvar_labels <- c(rep(c("means = (-5, 5)", "means = (-4, 4)", "means = (-3, 3)",
                         "means = (-2, 2)"), each = 1000))
ranmean_labels <- c(rep(c("sd = 0.5", "sd = 1", "sd = 1.5", "sd = 2"), each = 1000))

MC4_bi_ranvar_all <- cbind(MC4_bi_ranvar_all, label = ranvar_labels)
MC4_bi_ranmean_all <- cbind(MC4_bi_ranmean_all, label = ranmean_labels)

MC4_bi_ranvar_all <- mutate(MC4_bi_ranvar_all,
                            difference = scales::rescale(difference, to = c(0, 1)),
                            variance = scales::rescale(variance, to = c(0, 1)),
                            kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                            CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

MC4_bi_ranmean_all <- mutate(MC4_bi_ranmean_all,
                             difference = scales::rescale(difference,
                                                          to = c(0, 1)),
                             variance = scales::rescale(variance, to = c(0, 1)),
                             kurtosis = scales::rescale(kurtosis, to = c(0, 1)),
                             CPC = scales::rescale(CPC, to = c(0, 1))) %>%
  pivot_longer(cols = 2:5, names_to = "Measure",
               values_to = "value")

#plot simulation results (Figure S10)

ggplot(MC4_bi_ranvar_all, aes(x = ranvar, y = value, linetype = Measure,
                              color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Standard Deviations", y = "Estimated Level of Polarization")

ggplot(MC4_bi_ranmean_all, aes(x = ranmean, y = value, linetype = Measure,
                               color = Measure)) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_linetype_manual(values = c("solid", "longdash", "dotdash", "dotted")) +
  scale_color_manual(values = c("black", "gray50", "gray50", "gray50")) +
  facet_wrap(~ label) +
  labs(x = "Component Means", y = "Estimated Level of Polarization")

#calculate polarization ranks for random variance distributions

MC4_bi_ranvar1 <- mutate(MC4_bi_ranvar1, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-5, 5)", 1000))
MC4_bi_ranvar2 <- mutate(MC4_bi_ranvar2, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-4, 4)", 1000))
MC4_bi_ranvar3 <- mutate(MC4_bi_ranvar3, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance),
                         kurtosis = rank(-kurtosis),
                         CPC = rank(CPC),
                         mean_label = rep("means = (-3, 3)", 1000))
MC4_bi_ranvar4 <- mutate(MC4_bi_ranvar4, rank_true = rank(-ranvar),
                         difference = rank(difference),
                         variance = rank(variance), kurtosis = rank(-kurtosis),
                         CPC = rank(CPC), mean_label = rep("means = (-2, 2)", 1000))

#calculate polarization ranks for random mean distributions

MC4_bi_ranmean1 <- mutate(MC4_bi_ranmean1, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 0.5", 1000))
MC4_bi_ranmean2 <- mutate(MC4_bi_ranmean2, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 1", 1000))
MC4_bi_ranmean3 <- mutate(MC4_bi_ranmean3, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 1.5", 1000))
MC4_bi_ranmean4 <- mutate(MC4_bi_ranmean4, rank_true = rank(ranmean),
                          difference = rank(difference),
                          variance = rank(variance),
                          kurtosis = rank(-kurtosis),
                          CPC = rank(CPC),
                          variance_label = rep("sd = 2", 1000))

#combine polarization rank data into single data set for plotting

MC4_bi_ranvar_plot <- rbind(MC4_bi_ranvar1, MC4_bi_ranvar2, MC4_bi_ranvar3,
                            MC4_bi_ranvar4)
MC4_bi_ranmean_plot <- rbind(MC4_bi_ranmean1, MC4_bi_ranmean2, MC4_bi_ranmean3,
                             MC4_bi_ranmean4)

#tidy data and create labels for plotting

MC4_bi_ranvar_plot <- pivot_longer(MC4_bi_ranvar_plot, cols = 2:5,
                                   names_to = "measure", values_to = "value")
MC4_bi_ranmean_plot <- pivot_longer(MC4_bi_ranmean_plot, cols = 2:5,
                                    names_to = "measure", values_to = "value")

MC4_bi_ranvar_plot$mean_label <- factor(MC4_bi_ranvar_plot$mean_label,
                                        levels = c("means = (-2, 2)",
                                                   "means = (-3, 3)",
                                                   "means = (-4, 4)",
                                                   "means = (-5, 5)"))
MC4_bi_ranmean_plot$variance_label <- factor(MC4_bi_ranmean_plot$variance_label,
                                             levels = c("sd = 2", "sd = 1.5",
                                                        "sd = 1", "sd = 0.5"))
MC4_bi_ranvar_plot$measure <- factor(MC4_bi_ranvar_plot$measure,
                                     levels = c("difference", "variance",
                                                "kurtosis", "CPC"))
MC4_bi_ranmean_plot$measure <- factor(MC4_bi_ranmean_plot$measure,
                                      levels = c("difference", "variance",
                                                 "kurtosis", "CPC"))

#plot results of bivariate Monte Carlo simulations (Figure S12)

ggplot(MC4_bi_ranvar_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ mean_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(MC4_bi_ranmean_plot, aes(x = rank_true)) +
  geom_point(aes(y = value), size = 0.5) +
  xlab("True Polarization Rank") +
  ylab("Estimated Polarization Rank") +
  facet_grid(measure ~ variance_label) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###################################################################################
#RMSE and MAE
###################################################################################

#apply dynamic labels and combine distance and concentration data frames

MC_uni_ranmean_plot$dynamic <- "distance"
MC_uni_ranvar_plot$dynamic <- "concentration"
MC_bi_ranmean_plot$dynamic <- "distance"
MC_bi_ranvar_plot$dynamic <- "concentration"

MC3_uni_ranmean_plot$dynamic <- "distance"
MC3_uni_ranvar_plot$dynamic <- "concentration"
MC3_bi_ranmean_plot$dynamic <- "distance"
MC3_bi_ranvar_plot$dynamic <- "concentration"

MC4_uni_ranmean_plot$dynamic <- "distance"
MC4_uni_ranvar_plot$dynamic <- "concentration"
MC4_bi_ranmean_plot$dynamic <- "distance"
MC4_bi_ranvar_plot$dynamic <- "concentration"

MC_uni_both <- rbind(MC_uni_ranmean_plot[,c(2, 4:6)],
                     MC_uni_ranvar_plot[,c(2, 4:6)])
MC_bi_both <- rbind(MC_bi_ranmean_plot[,c(2, 4:6)],
                    MC_bi_ranvar_plot[,c(2, 4:6)])

MC3_uni_both <- rbind(MC3_uni_ranmean_plot[,c(2, 4:6)],
                      MC3_uni_ranvar_plot[,c(2, 4:6)])
MC3_bi_both <- rbind(MC3_bi_ranmean_plot[,c(2, 4:6)],
                     MC3_bi_ranvar_plot[,c(2, 4:6)])

MC4_uni_both <- rbind(MC4_uni_ranmean_plot[,c(2, 4:6)],
                      MC4_uni_ranvar_plot[,c(2, 4:6)])
MC4_bi_both <- rbind(MC4_bi_ranmean_plot[,c(2, 4:6)],
                     MC4_bi_ranvar_plot[,c(2, 4:6)])

#create separate data frames of all observations and apply dynamic labels

MC_uni_all <- rbind(MC_uni_ranmean_plot[,c(2, 4:6)],
                    MC_uni_ranvar_plot[,c(2, 4:6)]) %>%
  mutate(dynamic = "all")
MC_bi_all <- rbind(MC_bi_ranmean_plot[,c(2, 4:6)],
                   MC_bi_ranvar_plot[,c(2, 4:6)]) %>%
  mutate(dynamic = "all")

MC3_uni_all <- rbind(MC3_uni_ranmean_plot[,c(2, 4:6)],
                     MC3_uni_ranvar_plot[,c(2, 4:6)]) %>%
  mutate(dynamic = "all")
MC3_bi_all <- rbind(MC3_bi_ranmean_plot[,c(2, 4:6)],
                    MC3_bi_ranvar_plot[,c(2, 4:6)]) %>%
  mutate(dynamic = "all")

MC4_uni_all <- rbind(MC4_uni_ranmean_plot[,c(2, 4:6)],
                     MC4_uni_ranvar_plot[,c(2, 4:6)]) %>%
  mutate(dynamic = "all")
MC4_bi_all <- rbind(MC4_bi_ranmean_plot[,c(2, 4:6)],
                    MC4_bi_ranvar_plot[,c(2, 4:6)]) %>%
  mutate(dynamic = "all")

#create unified data sets for both univariate and bivariate distributions and apply dimension labels

MC_uni <- rbind(MC_uni_both, MC_uni_all) %>%
  mutate(dimensions = "univariate")
MC_bi <- rbind(MC_bi_both, MC_bi_all) %>%
  mutate(dimensions = "bivariate")

MC3_uni <- rbind(MC3_uni_both, MC3_uni_all) %>%
  mutate(dimensions = "univariate")
MC3_bi <- rbind(MC3_bi_both, MC3_bi_all) %>%
  mutate(dimensions = "bivariate")

MC4_uni <- rbind(MC4_uni_both, MC4_uni_all) %>%
  mutate(dimensions = "univariate")
MC4_bi <- rbind(MC4_bi_both, MC4_bi_all) %>%
  mutate(dimensions = "bivariate")

#create unified data set of Monte Carlo predictions and rankings

MC_all <- rbind(MC_uni, MC_bi)
MC3_all <- rbind(MC3_uni, MC3_bi)
MC4_all <- rbind(MC4_uni, MC4_bi)

#calculate RMSE, MAE and percent predicted for univariate and bivariate distance, concentration, and total (Table 1)

MC_metrics <- group_by(MC_all, dimensions, dynamic, measure) %>%
  dplyr::summarize(RMSE = round(Metrics::rmse(rank_true, value), 2),
                   MAE = round(Metrics:: mae(rank_true, value), 2))

MC3_metrics <- group_by(MC3_all, dimensions, dynamic, measure) %>%
  dplyr::summarize(RMSE = round(Metrics::rmse(rank_true, value), 2),
                   MAE = round(Metrics:: mae(rank_true, value), 2))

MC4_metrics <- group_by(MC4_all, dimensions, dynamic, measure) %>%
  dplyr::summarize(RMSE = round(Metrics::rmse(rank_true, value), 2),
                   MAE = round(Metrics:: mae(rank_true, value), 2))
