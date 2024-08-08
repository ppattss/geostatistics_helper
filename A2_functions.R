########################################################################################
########################################################################################
#
#     Project:  SPA553 Assessment 2 - function definitions
#     Author:   Peter Patterson
#     Date:     18/06/2024
#
#
#   = assignment operator creates local variable
#   <- assignment operator creates a global variable
#   'vg' = variogram object
#   'vgm' = variogram model object
#   object variables cannot begin with a number, eg not test$1, do test$a instead
#   functions need to be initialised before calling
#
#
#   debug_function()
#   eda_analysis(eda_input, eda_attribute)
#   new_variogram_model(df, type, lag, cutoff, nugget, psill, range, name)
#   ok_function()
#   sk_function()
#   uk_function()
#
########################################################################################
########################################################################################

print("Loading helper functions ...")


debug_function <- function() {
  print("Debug function output")

}

# ---------

summary_stats <- function(data) {
  summary_data <- numSummary(data, statistics=c("mean", "sd", "se(mean)", "var", "IQR", "quantiles", "CV", "skewness", "kurtosis"), quantiles=c(0,.025,.25,.5,.75,.975,1), type="2")
  return(summary_data$table)
}


# --------

create_map <- function(est_data, source_data, map_colour, map_limits, title) {
  this_map <- NULL
  this_map <- ggplot() +
              geom_stars(data = est_data, aes(fill = var1.pred, x = x, y = y)) +
              xlab(NULL) +
              ylab(NULL) +
              geom_sf(data = source_data,aes(size=pb), pch=21) +
              scale_fill_viridis_c(option = map_colour, limits = map_limits, na.value = "transparent", name="Pb (ppm)") +
              labs(x = "x", y = "y", title = title)
  return(this_map)
}


########################################################################################
#     Exploratory Data Analysis
########################################################################################

# numSummary(input_dataframe$pb, statistics=c("mean", "sd", "se(mean)", "var", "IQR", "quantiles", "CV", "skewness", "kurtosis"), quantiles=c(0,.025,.25,.5,.75,.975,1), type="2")

eda_analysis <-function(eda_input, eda_attribute = "Pb (ppm)", hist_bins) {

  print(summary(eda_input))

  eda_data <- NULL

  eda_data$variance <- round(var(eda_input$pb))
  print(paste("Sampled Data Variance = ", toString(eda_data$variance)))

  eda_data$mean <- mean(eda_input$pb)
  print(paste("Sampled Data Mean = ", toString(eda_data$mean)))

  eda_data$title <- "Spatial Distribution of Pb Samples"
  eda_data$color <- "H"
  #eda_data$label <- "Pb (ppm)"
  eda_data$min <- 0
  eda_data$max <- 120
  eda_data$limits <- c(eda_data$min, eda_data$max)

  eda_data$spatial_dist <- ggplot(data = eda_input) + geom_sf(aes(color = pb,size=pb)) + scale_color_viridis_c(option = eda_data$color , limits = eda_data$limits, name = eda_attribute) + labs(x = "x", y = "y", title= eda_data$title)

  eda_data$histogram <- with(eda_input, hist(pb, freq = TRUE, breaks=hist_bins, col="darkgray",xlab=eda_attribute, xlim=c(0, 120), main="Sampled Data Histogram", plot = TRUE))

  return (eda_data)
}


########################################################################################
#     Variogram Model & Plots
########################################################################################

#   new_variogram_model(input_df = , type, = ,direction = , lag = , nugget = , psill = , range = , cutoff = , name = "")
new_variogram_model <- function(data, type, lag, cutoff, nugget, psill, range, name) {
  #width <- lag

  this_model = NULL
  this_model$variogram = variogram(pb~1, data, width=lag, cutoff=cutoff)
  this_model$variogram_plot = plot(this_model$variogram, plot.numbers = T, pch = 19, col = "red", ylim=c(0, 1200))

  this_model$variogram_model = vgm(psill, model = type, range, nugget)
  this_model$variogram_model_plot = plot(this_model$variogram, pl = T, model = this_model$variogram_model, pch = 19, col = "blue", main = name,ylim=c(0, 1200))

  cat(name, "\n", "Type", "Lag", "Cutoff", "Nugget", "pSill", "Sill", "Range\n", type, lag, cutoff, nugget, psill, range, sep = "\t")
  plot(this_model$variogram_model_plot)

  return (this_model)
}


########################################################################################
#     Ordinary Kriging
########################################################################################

ordinary_kriging <- function(search, ok_input, this_variogram_model, this_grid, title, est_label, est_color, var_label, var_color, est_limits, ind_label) {

  ok_data <- NULL

  # ok_data$maxDistance <- 1200 # maximum distance to look for sampled data (usually a bit further then the range of variogram model)
  # ok_data$nMin <- 30 # minimum number of sampled data to calculate the estimate
  # ok_data$nMax <- 60 # maximum number of sampled data to calculate the estimate
  # ok_data$oMax <- 3 # maximum number of sampled data per octant to calculate the estimate

  #ok_data$krige <- krige(pb ~ 1, locations = ok_input, newdata = this_grid, model = this_variogram_model, maxdist = ok_data$maxDistance,nmin = ok_data$nMin, nmax= ok_data$nMax, omax=ok_data$oMax)

  ok_data$krige <- krige(pb ~ 1, locations = ok_input, newdata = this_grid, model = this_variogram_model, maxdist = search$maxDist,nmin = search$nMin, nmax= search$nMax, omax=search$oMax)

  #ok_data$krige_cv <- krige.cv(pb ~ 1, ok_input, model=this_variogram_model, maxdist = search$maxDist,nmin = search$nMin, nmax= search$nMax, omax=search$oMax)

  ok_data$krige_cv <- krige.cv(pb ~ 1, ok_input, model=this_variogram_model, nmax= search$nMax)

  #LOOCV nmax
  #cv_exp_ok_leave <- krige.cv(pb ~ 1, input_dataframe, model=vgm_exponential$variogram_model, nmax=32)
  #K-fold cross-validation (K=10)
  #cv_exp_ok_kfold <- krige.cv(pb ~ 1, input_dataframe, model=vgm_exponential$variogram_model, nmax=32, nfold=10)

  print(title)
  print(ok_data$krige)

  #OK estimates
  ok_data$estMin <- 19#min(ok_data$krige$var1.pred)
  ok_data$estMax <- 83#max(ok_data$krige$var1.pred)
  ok_data$estLimits <- c(ok_data$estMin, ok_data$estMax)
  ok_data$estMap <- ggplot() + geom_stars(data = ok_data$krige, aes(fill = var1.pred, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = ok_input,aes(size=pb), pch=21) + scale_fill_viridis_c(option = est_color, limits = est_limits, na.value = "transparent", name = est_label) + labs(x = "x", y = "y", title=title)

  #OK variance
  ok_data$varMin <- 408#min(ok_data$krige$var1.var)
  ok_data$varMax <- 1028#max(ok_data$krige$var1.var)
  ok_data$varLimits <- c(ok_data$varMin, ok_data$varMax)
  ok_data$varMap <- ggplot() + geom_stars(data = ok_data$krige, aes(fill = var1.var, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = ok_input, aes(size=pb), pch=21) + scale_fill_viridis_c(option = var_color, limits = c(157, 1191), na.value = "transparent", name = "Variance") + labs(x = "x", y = "y", title = var_label)

  #grid.arrange(ok_data$estMap, ok_data$varMap, nrow=1)
  plot(ok_data$estMap)

  ind_value <- 45.0
  ik_data <- ok_data$krige$var1.pred
  ik_results <- ok_data$krige$var1.pred
  for (x in 1:20) {
    #print(paste("X = ",x))
    for (y in 1:20) {
      #print(paste("X = ", x," Y = ", y))
      this_cell = ok_data$krige$var1.pred[x,y]
      #print(this_cell)
      #ik_data[x,y] <- ifelse(this_cell < ind_value, 1,0)
      #ik_results[x,y] <- ifelse(this_cell < 45.0, 1,0)
      #ik_results[x,y] <- ifelse(this_cell < 40, 1, ifelse(this_cell < 50, 0.5, 0))
      ik_results[x,y] <- ifelse(this_cell < 40, "Min - 40 ppm", ifelse(this_cell < 45, "40 - 45 ppm", ifelse(this_cell < 50, "45 - 50 ppm","50 - Max ppm")))
    }
  }
  ok_data$krige$var1.ind <- ik_results
  ind_values <- c("0" = "blue", "1" = "red")

  ok_data$indMap <- ggplot() +
    geom_stars(data = ok_data$krige, aes(fill = var1.ind, x = x, y = y)) +
    xlab(NULL) +
    ylab(NULL) +
    geom_sf(data = ok_input, pch=21) +
    scale_color_manual(values = ind_values) +
    labs(color = "Indicator", x="x", y="y", title=ind_label)

  #plot(ok_data$indMap)


  return (ok_data)
}


########################################################################################
#     Simple Kriging
########################################################################################

simple_kriging <- function(sk_mean_beta, sk_input, this_variogram_model, this_grid, title, est_label, est_color, var_label, var_color, est_limits, ind_label) {
  sk_data <- NULL
  # #input_mean <- mean(input_dataframe$pb)

  sk_data$krige <- krige(pb ~ 1, locations = sk_input, newdata = this_grid, model = this_variogram_model, beta= sk_mean_beta)

  sk_data$krige_cv <- krige.cv(pb ~ 1, sk_input, model=this_variogram_model, beta= sk_mean_beta)

  print(title)
  print(sk_data$krige)

  #SK estimates
  sk_data$estMin <- 19#min(sk_data$krige$var1.pred)
  sk_data$estMax <- 83#max(sk_data$krige$var1.pred)
  sk_data$estLimits <- c(sk_data$estMin, sk_data$estMax)

  sk_data$estMap <- ggplot() + geom_stars(data = sk_data$krige, aes(fill = var1.pred, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = sk_input,aes(size=pb), pch=21) + scale_fill_viridis_c(option = est_color, limits = est_limits, na.value = "transparent", name = est_label) + labs(x = "x", y = "y", title = title)

  # SK variance
  sk_data$varMin <- 408#min(sk_data$krige$var1.var)
  sk_data$varMax <- 1028#max(sk_data$krige$var1.var)
  sk_data$varLimits <- c(sk_data$varMin, sk_data$varMax)

  sk_data$varMap <- ggplot() + geom_stars(data = sk_data$krige, aes(fill = var1.var, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = sk_input,aes(size=pb), pch=21) + scale_fill_viridis_c(option = var_color, limits = c(157, 1191), na.value = "transparent", name = "Variance") + labs(x = "x", y = "y", title = var_label)

  plot(sk_data$estMap)

  ind_value <- 45.0
  ik_data <- sk_data$krige$var1.pred
  ik_results <- sk_data$krige$var1.pred
  for (x in 1:20) {
    #print(paste("X = ",x))
    for (y in 1:20) {
      #print(paste("X = ", x," Y = ", y))
      this_cell = ik_data[x,y]
      #print(this_cell)
      #ik_data[x,y] <- ifelse(this_cell < ind_value, 1,0)
      #ik_results[x,y] <- ifelse(this_cell < 45.0, 1,0)
      #ik_results[x,y] <- ifelse(this_cell < 40, 1, ifelse(this_cell < 50, 0.5, 0))
      ik_results[x,y] <- ifelse(this_cell < 40, "Min - 40 ppm", ifelse(this_cell < 45, "40 - 45 ppm", ifelse(this_cell < 50, "45 - 50 ppm","50 - Max ppm")))
    }
  }
  sk_data$krige$var1.ind <- ik_results
  ind_values <- c("0" = "blue", "1" = "red")

  sk_data$indMap <- ggplot() +
    geom_stars(data = sk_data$krige, aes(fill = var1.ind, x = x, y = y)) +
    xlab(NULL) +
    ylab(NULL) +
    geom_sf(data = sk_input, pch=21) +
    scale_color_manual(values = ind_values) +
    labs(color = "Indicator", x="x", y="y", title=ind_label)

  #plot(sk_data$indMap)

  return (sk_data)
}


########################################################################################
#     Universal Kriging
########################################################################################

universal_kriging <- function(uk_input, this_variogram_model, this_grid, title, est_label, est_color, var_label, var_color, est_limits, ind_label) {
  uk_data <- NULL

  #General formula for omnidirectional residuals UK
  uk_data$krige <- krige(pb~x + y, locations = uk_input, newdata = this_grid, model = this_variogram_model)

  # Cross-validation LOOCV
  uk_data$krige_cv <- krige.cv(pb ~ 1, uk_input, model=this_variogram_model)

  # UK estimate
  #uk_data$estMap <- ggplot() + geom_stars(data = uk_data$krige, aes(fill = var1.pred, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = uk_input,aes(size=pb), pch=21) + scale_fill_viridis_c(option = est_color, limits = est_limits, na.value = "transparent", name="Pb (ppm)") + labs(x = "x", y = "y", title = title)
  print("trying function")
  uk_data$estMap <- create_map(uk_data$krige, uk_input, est_color, est_limits, title)
  print("finsihed function")

  # UK variance
  uk_data$varMap <- ggplot() + geom_stars(data = uk_data$krige, aes(fill = var1.var, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = uk_input,aes(size=pb), pch=21) + scale_fill_viridis_c(option = var_color, limits = c(157, 1191), na.value = "transparent", name="Variance") + labs(x = "x", y = "y", title = var_label)

  print(title)
  print(uk_data$krige)
  plot(uk_data$estMap)
  print("mapplotted")

  ind_value <- 45.0
  ik_data <- uk_data$krige$var1.pred
  ik_results <- uk_data$krige$var1.pred
  for (x in 1:20) {
    #print(paste("X = ",x))
    for (y in 1:20) {
      #print(paste("X = ", x," Y = ", y))
      this_cell = ik_data[x,y]
      #print(this_cell)
      #ik_data[x,y] <- ifelse(this_cell < ind_value, 1,0)
      #ik_results[x,y] <- ifelse(this_cell < 45.0, 1,0)


      #ik_results[x,y] <- ifelse(this_cell < 40,3, ifelse(this_cell < 45, 2, ifelse(this_cell < 50, 1,0)))
      ik_results[x,y] <- ifelse(this_cell < 40, "Min - 40 ppm", ifelse(this_cell < 45, "40 - 45 ppm", ifelse(this_cell < 50, "45 - 50 ppm","50 - Max ppm")))

    }
  }
  uk_data$krige$var1.ind <- ik_results
  ind_values <- c("< 40 ppm" = "blue", "40 - 45 ppm" = "lightblue", "45 - 50 ppm", "salmon", "> 50 ppm", "red")

  uk_data$indMap <- ggplot() +
    geom_stars(data = uk_data$krige, aes(fill = var1.ind, x = x, y = y)) +
    xlab(NULL) +
    ylab(NULL) +
    geom_sf(data = uk_input, pch=21) +
    scale_color_manual(values = ind_values) +
    labs(color = "Indicator", x="x", y="y", title=ind_label)

  #plot(uk_data$indMap)

  return(uk_data)
}


########################################################################################
#     Indicator Kriging
########################################################################################

indicator_kriging <- function(ik_input) {
  ik_data <- NULL


  return(ik_data)
}


########################################################################################
#     Cross-Validation
########################################################################################

cross_validation <- function(cv_input) {
  cv_data <- NULL


  return(cv_data)
}



########################################################################################
#
########################################################################################


########################################################################################
#
########################################################################################

