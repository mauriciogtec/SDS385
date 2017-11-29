Rcpp::sourceCpp('Solutions08/gaussian_filter.cpp')

library(tidyverse)
library(RColorBrewer)
library(leaflet)
library(sp)

concentration <-  read_csv("https://raw.githubusercontent.com/jgscott/SDS385/master/data/co2.csv")

x <- data.matrix(concentration[ , c(2,3)])
z <- data.matrix(concentration[ , 4])

system.time({ # 6 minutes
  z2 <- gaussian_filter(x, z, 5, ncores = 8)
})
concentration$co2avgret_smooth <- z2[, 1]

plot(z2[ ,1], z[ ,1], xlab = "original", "averaged")

# myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 8))

# ggplot(concentration, aes(x = lon, y = lat, colour = co2avgret)) +
#   geom_point(size = 0.1) + theme_minimal() +
#   scale_colour_gradientn(colours = heat.colors(10))

n <- 100000
train_ratio <- 0.8
samp <- sample(nrow(concentration), n)

trainid <- samp[1:floor(train_ratio * n)]
xtrain <- data.matrix(concentration[trainid, c(2,3)])
ztrain <- data.matrix(concentration[trainid, 4])

b <- 10
z2 <- gaussian_filter(xtrain, ztrain, 1, ncores = 8)

df <- SpatialPointsDataFrame(
  data.frame(xtrain),
  data.frame(ztrain)
)

pal <- colorNumeric("RdYlBu", domain = df@data$co2avgret)
leaflet(df) %>% 
  addTiles() %>%
  addCircleMarkers(
    radius = 2, 
    stroke = FALSE, 
    color = ~pal(co2avgret), 
    label = ~as.character(co2avgret)
    ) %>% 
  addLegend(
    "bottomright",
    values = ~co2avgret,
    pal = pal)

# 
# # Number of folds for cross-validation
# nfolds <- 5
# N <- nrow(xtrain)
# 
# # Create 10 equal size folds as a list of validation indexes
# folds <- split(sample(N), cut(1:N, breaks = nfolds, labels = FALSE))
# 
# # Run cross-validation
# mse_data <- data.frame()
# for (k in 1:nfolds) {
#   idk <- folds[[k]]
#   fitk <- gaussian_filter(xtrain[-idk, ], y_train[-idk], lambda = lambda_seq)
#   mse_cv <- sapply(1:nlambda, function(i) {
#     mean((y_train[idk] - predict(fitk, X_train[idk, ], s = lambda_seq[i]))^2)
#   }) 
#   mse_data <- rbind(mse_data, data.frame(
#     fold = k,
#     lambda = lambda_seq,
#     mse_cv = mse_cv
#   ))
# }
