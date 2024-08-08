########################################################################################
########################################################################################
#
#     Project:  SPA553 Assessment 2 - processing script
#     Author:   Peter Patterson
#     Date:     18/06/2024
#
#   = assignment operator creates local variable
#   <- assignment operator creates a global variable
#   <<- assignment operator forces variable to be global, even when declared in a function
#   object variables cannot begin with a number, eg can not object$1, must do object$a instead
#   functions need to be initialised before calling
#   dedicated config file to preload functions
#   eda_results["variance"] === eda_results$variance
#
#   control+click on a function call to lgo to its definition. (command+click on macOS)
#   access function documentation using ?FUNCTION, eg ?ggplot

#   'vg' = variogram object
#   'vgm' = variogram model object
#   'lm' = linear model (used to fit linear models, or carry out regression)
#   'var' = variance

########################################################################################
########################################################################################

#source("/Users/peterpatterson/Documents/r_data/Config.R")

on_load <- c(setwd("/Users/peterpatterson/Documents/r_data/"),
             source("A2_header.R"),

             input_path <- "/Users/peterpatterson/Documents/r_data/HM_Pb_R.csv",
             input_dataframe <- st_as_sf(as.data.frame(read.csv(input_path)), coords = c("x", "y"), remove = FALSE),

             grid_resolution <- 50,  #grid resolution value in metres
             kgrid <- st_as_stars(st_bbox(input_dataframe), dx = grid_resolution), #Create a grid resolution dx

             attribtue_label <- "Pb (ppm)",
             estimate_color_palette <- "C",
             variance_color_palette<- "viridis",

             estimate_min <- 17,
             estimate_max <- 85,
             est_limits <- c(estimate_min, estimate_max), #c(17,85)
             variance_min <- 157,
             variance_max <- 1191,
             var_limits <- c(variance_min,variance_max),  #c(408,1189). # c(157, 1191)
             ok_range <- c(22,88),
             sk_range <- c(22,88),
             uk_range <- c(16,89),
             hist_range <- c(15,90)
)



input_raw <- (read.csv(input_path))
force_df <- st_as_sf(as.data.frame(read.csv(input_path)), coords = c("x", "y"))
#test_df <- st_as_sf(read.csv(input_path), coords = c("x", "y"), remove=FALSE)
input_datafram.grid <- st_as_stars(st_bbox(input_dataframe), dx = grid_resolution) #Create a grid resolution dx



#print.data.frame(input_dataframe)
#print.summaryDefault(variogram_model_a)

#
# data(meuse)
# coordinates(meuse) = ~x+y
#
# data(meuse.grid)
#
# gridded(meuse.grid) = ~x+y
#
# sp_df <- (read.csv(input_path))
# data(sp_df)<-(read.csv(input_path))
#


########################################################################################
# Exploratory Data Analysis -   eda_analysis(eda_input, eda_attribute = "Pb (ppm)", hist_bins)


eda_results <- eda_analysis(input_dataframe, hist_bins = 24)

eda_results$spatial_dist

plot(eda_results$histogram, xlim = c(0,120), main = "Sampled Data", xlab = attribtue_label)

plot(eda_results$histogram, main = "Sampled Data", xlab = attribtue_label)


 this_figure <- c(
    # Layout to split the screen
    layout(mat = matrix(c(1,2),2,1, byrow=TRUE),  height = c(1,8)),

    par(mar=c(0, 4, 1.1, 2.1)),
    boxplot(input_dataframe$pb , horizontal=TRUE , ylim=c(0,120), xaxt="n" , col=rgb(0.8,0.8,0,0.5) , frame=F),

    par(mar=c(4, 4, 1.1, 2.1)),
    hist(input_dataframe$pb, freq = TRUE, breaks=24 , col="darkgray" , border="black" , main="Sampled Data" , xlab="Pb (ppm)", xlim=c(0,120), ylab ="Frequency"),

    abline(v = mean(input_dataframe$pb), col='red', lwd = 1),
    abline(v = median(input_dataframe$pb), col='blue', lwd = 1),

    this_quantile <- quantile(input_dataframe$pb, probs = c(0,0.25,0.5, 0.75, 1)),

    str_min <- paste("Min:\t\t\t", this_quantile[1]),
    str_q1 <- paste("1st Qu.:\t", this_quantile[2]),
    str_mean <- paste("Mean:\t\t", round(mean(input_dataframe$pb), 2)),
    str_median <- paste("Median:\t", this_quantile[3]),
    str_q3 <- paste("3rd Qu.:\t", this_quantile[4]),
    str_max <- paste("Max:\t\t\t", this_quantile[5])

)


legend("topright",
       legend = c(str_min, str_q1,str_mean, str_median, str_q3,str_max),
       col = c(NA, NA,"red", "blue", NA, NA),
       lwd = 2,
       lty = 1
)

########################################################################################
# Variogram Models (vgm) -    new_variogram_model(input_dataframe, name, type, lag, cutoff, nugget, psill, range)

# spherical_parameters <- data.frame(name = "Spherical (omnidirectional)",
#                                    type = "Sph",
#                                    lag = 100,
#                                    cutoff = 1000,
#                                    nugget = 300,
#                                    psill = 536,
#                                    range = 465,
#                                    direction = data.frame(azimuth = 0,
#                                                           dip = 0,
#                                                           tolerence = 180,
#                                                           bandwidth = 100
#                                                           )
#                                    )
# print(spherical_parameters)

vgm_spherical <-  new_variogram_model(data = input_dataframe,
                                      name = "Spherical (omnidirectional)",
                                      #direction="(omnidirectional)",
                                      type="Sph",
                                      lag=100,
                                      cutoff = 1000,
                                      nugget=300,
                                      psill=536,
                                      range = 465 )

vgm_exponential <- new_variogram_model(data = input_dataframe,
                                     name = "Exponential (omnidirectional)",
                                     #direction="(omnidirectional)",
                                     type="Exp",
                                     lag=100,
                                     cutoff = 1000,
                                     nugget=100,
                                     psill=736,
                                     range = 420 )

vgm_gaussian <-  new_variogram_model(data = input_dataframe,
                                     name = "Gaussian (omnidirectional)",
                                     #direction="(omnidirectional)",
                                     type="Gau",
                                     lag=100,
                                     cutoff = 1000,
                                     nugget=400,
                                     psill=436,
                                     range = 420 )

grid.arrange(vgm_spherical$variogram_model_plot, vgm_exponential$variogram_model_plot, vgm_gaussian$variogram_model_plot, nrow=1)

########################################################################################
# Ordinary Kirging Estimation - ok_function(search, input_dataframe, variogram_model, grid, title, est_label, est_color, var_label, var_color)

ok_search_parameters <- data.frame(maxDist = 1200,
                                   nMin = 30,
                                   nMax = 60,
                                   oMax = 3 )

sph_ok_results <- ordinary_kriging(ok_search_parameters, input_dataframe, vgm_spherical$variogram_model, kgrid, "Spherical OK Estimate", attribtue_label, estimate_color_palette, "Spherical OK Variance", variance_color_palette, uk_range, "Spherical OK Indicator Map")

exp_ok_results <- ordinary_kriging(ok_search_parameters, input_dataframe, vgm_exponential$variogram_model, kgrid, "Exponential OK Estimate", attribtue_label, estimate_color_palette, "Exponential OK Variance", variance_color_palette, uk_range, "Spherical OK Indicator Map")

gau_ok_results <- ordinary_kriging(ok_search_parameters, input_dataframe, vgm_gaussian$variogram_model, kgrid, "Gaussian OK Estimate", attribtue_label, estimate_color_palette, "Gaussian OK Variance", variance_color_palette, uk_range, "Spherical OK Indicator Map")

grid.arrange(sph_ok_results$estMap, exp_ok_results$estMap, gau_ok_results$estMap, nrow = 1)

#sph_ok_results$varMap

########################################################################################
# Simple Kriging - sk_function(sk_mean_beta, input_dataframe, variogram_model, grid, title, est_label, est_color, var_label, var_color, limits)

#sk_mean_beta <- eda_results$mean
sk_mean_beta <- round(eda_results$mean, digits = 5)

sph_sk_results <- simple_kriging(sk_mean_beta, input_dataframe, vgm_spherical$variogram_model, kgrid, "Spherical SK Estimate", attribtue_label, estimate_color_palette, "Spherical SK Variance", variance_color_palette, uk_range, "Spherical SK Indicator Map")

exp_sk_results <- simple_kriging(sk_mean_beta, input_dataframe, vgm_exponential$variogram_model, kgrid, "Exponential SK Estimate", attribtue_label, estimate_color_palette, "Exponential SK Variance", variance_color_palette, uk_range, "Exponential SK Indicator Map")

gau_sk_results <- simple_kriging(sk_mean_beta, input_dataframe, vgm_gaussian$variogram_model, kgrid, "Gaussian SK Estimate", attribtue_label, estimate_color_palette, "Gaussian SK Variance", variance_color_palette, uk_range, "Gaussian SK Indicator Map")

grid.arrange(sph_sk_results$estMap, exp_sk_results$estMap, gau_sk_results$estMap, nrow = 1)


########################################################################################
# Universal Kriging -


##Add x,y coordinates as covariates in the spatial dataset to be used as trend
input_dataframe$x = st_coordinates(input_dataframe)[,1]
input_dataframe$y = st_coordinates(input_dataframe)[,2]

#Get residuals from x,y trend to evaluate sill
uk_residuals <- lm(pb~x+y, input_dataframe)

uk_variance <- var(uk_residuals$residuals) #residuals variance as indicative sill to determine psill and nugget
paste("UK Variacnce: ", uk_variance)


# ----- apply UK krige

sph_uk_results <- universal_kriging(input_dataframe, vgm_spherical$variogram_model, kgrid, "Spherical UK Estimate", attribtue_label, estimate_color_palette, "Spherical UK Variance", variance_color_palette, uk_range, "Spherical UK Indicator Map")

exp_uk_results <- universal_kriging(input_dataframe, vgm_exponential$variogram_model, kgrid, "Universal Kriging Estimation (Exponential Model)", attribtue_label, estimate_color_palette, "Exponential UK Variance", variance_color_palette, uk_range, "Exponential UK Indicator Map")

gau_uk_results <- universal_kriging(input_dataframe, vgm_gaussian$variogram_model, kgrid, "Gaussian UK Estimate", attribtue_label, estimate_color_palette, "Gaussian UK Variance", variance_color_palette, uk_range, "Gaussian UK Indicator Map")


#exp_uk_results$estMap
grid.arrange(sph_sk_results$estMap, exp_sk_results$estMap, gau_sk_results$estMap, nrow = 1)


#
# # ------ plot results
# #
# # uk_min <-  min(uk_lead_sph$var1.pred)
# # uk_max <- max(uk_lead_sph$var1.pred)
# #
# # legend_limits <- c(uk_min,uk_max)
#
# #legend_limits <- c(0, 120)
#
# uk_data_sph <- NULL
#
# #General formula for omnidirectional residuals UK
# uk_data_sph$krige = krige(pb~x + y, locations = input_dataframe, newdata = kgrid, model = vgm_spherical$variogram_model)
#
# uk_data_sph$krige
#
# uk_data_sph$estMap <- ggplot() + geom_stars(data = uk_data_sph$krige, aes(fill = var1.pred, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = input_dataframe,aes(size=pb), pch=21) + scale_fill_viridis_c(option = "C", limits = uk_range, na.value = "transparent", name="Pb (ppm)") + labs(x = "x", y = "y", title = "Spherical UK Estimate")
#
# plot(uk_data_sph$estMap)
#
# #UK variance
# uk_data_sph$varMap <- ggplot() + geom_stars(data = uk_lead_sph, aes(fill = var1.var, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = input_dataframe,aes(size=pb), pch=21) + scale_fill_viridis_c(option = "viridis", limits = var_limits, na.value = "transparent", name="Variance") + labs(x = "x", y = "y", title = "Exponential UK Variance")
#
#
# # -------  exponential models.
# uk_data_exp <- NULL
#
# uk_data_exp$krige = krige(pb~x + y, locations = input_dataframe, newdata = kgrid, model = vgm_exponential$variogram_model)
# uk_data_exp$krige
#
#
# uk_data_exp$estMap <- ggplot() + geom_stars(data = uk_data_exp$krige, aes(fill = var1.pred, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = input_dataframe,aes(size=pb), pch=21) + scale_fill_viridis_c(option = "C", limits = uk_range, na.value = "transparent", name="Pb (ppm)") + labs(x = "x", y = "y", title = "Exponential UK Estimate")
#
# plot(uk_data_exp$estMap)
#
# #UK variance
# uk_data_exp$varMap <- ggplot() + geom_stars(data =  uk_data_exp$krige, aes(fill = var1.var, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = input_dataframe,aes(size=pb), pch=21) + scale_fill_viridis_c(option = "viridis", limits = var_limits, na.value = "transparent", name="Variance") + labs(x = "x", y = "y", title = "Exponential UK Variance")
#
# plot(uk_data_exp$varMap)
#
#
#
# # -------  gaussian models.
# uk_data_gau <- NULL
#
#
# uk_data_gau$krige = krige(pb~x + y, locations = input_dataframe, newdata = kgrid, model = vgm_gaussian$variogram_model)
# uk_data_gau$krige
#
#
#
# uk_data_gau$estMap <- ggplot() + geom_stars(data = uk_data_gau$krige, aes(fill = var1.pred, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = input_dataframe,aes(size=pb), pch=21) + scale_fill_viridis_c(option = "C", limits = uk_range, na.value = "transparent", name="Pb (ppm)") + labs(x = "x", y = "y", title = "Gaussian UK Estimate")
#
# plot(uk_data_gau$estMap)
#
# #UK variance
# uk_data_gau$varMap <- ggplot() + geom_stars(data = uk_data_gau$krige, aes(fill = var1.var, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = input_dataframe,aes(size=pb), pch=21) + scale_fill_viridis_c(option = "viridis", limits = var_limits, na.value = "transparent", name="Variance") + labs(x = "x", y = "y", title = "Gaussian UK Variance")
#
#
#
#
# grid.arrange(uk_data_sph$estMap,uk_data_exp$estMap, uk_data_gau$estMap, nrow=1 )
#

########################################################################################
#         Indicator Krigign
########################################################################################


ik_threshold <- 45

ik_cellsize <- 50
#ik_cells <- (1000 / ik_cellsize)^2

ik_input <- exp_uk_results$krige$var1.pred
ik_cells <- length(ik_input)
typeof(ik_input)
ik_input[1,1]

for (double in ik_input) {
  print(double)
  #print(ik_input[double])

}

ind_value <- 45
ik_data <- ik_input

for (x in 1:20) {
  #print(paste("X = ",x))

  for (y in 1:20) {
    print(paste("X = ", x," Y = ", y))
    this_cell = ik_input[x,y]

    ik_data[x,y] <- ifelse(this_cell < ind_value, 1,0)

  }
}

exp_uk_results$krige$ind <- ik_data

exp_uk_results$indMap <- ggplot() + geom_stars(data = exp_uk_results$krige, aes(fill = ind, x = x, y = y)) + xlab(NULL) + ylab(NULL) + geom_sf(data = input_dataframe,aes(size=pb), pch=21) + scale_fill_viridis_c(option = "H") + labs(x = "x", y = "y", title = "title")


exp_uk_results$indMap
#Create indicator variable based on a relevant threshold
ind_value <- 45

input_dataframe$exp_uk_pb_binary <- ifelse(exp_uk_results$krige$var1.pred < ind_value, 1,0)
input_dataframe



input_dataframe$lead_binary <- ifelse(input_dataframe$pb < ind_value, 1,0)



#Spatial distribution Indicator variable
ggplot(data = input_dataframe) +
  geom_sf(aes(color = factor(lead_binary)), size = 3) +
  scale_color_manual(values = c("0" = "blue", "1" = "red")) +
  labs(color = "Indicator", x="x", y="y", title="Indicator Kriging (Pb Threshold = 45 mg/kg)")




#Calculate the experimental indicator variogram

#This is easier to calculate directly in R

#Try to get the experimental variogram with the default parameters
#to then specify lag and cutoff if needed

#
# #Default experimental indicator variogram
# vi <- variogram(lead_binary ~ 1, location=spatial_lead)
# plot(vi, plot.numbers = T,pch=19,col="red")
#
# #Experimental indicator variogram with lag and cutoff
# lag =   100     # lag distance (m)
# cutoff =   1000  # maximum distance for plotting pairs of points (m)
#
# vi <- variogram(lead_binary ~ 1, location=spatial_lead,width=lag, cutoff=cutoff)
# plot(vi, plot.numbers = T,pch=19,col="red")
#
# #Manual variogram model fitting
#
# var(spatial_lead$lead_binary) #variance as indicative sill to determine psill and nugget
#
# psill =  535    #partial sill
# range =   465   # continuity range
# nugget =  300   #nugget effect
#
# #Remember to use only the "model" that applies
# model = "Sph"
# #model="Exp"
#
# vim <- vgm(psill=psill,model= model,range=range,nugget=nugget)
# plot(vi, pl = T, model = vim,pch=19,col="blue")
#
# #Specifying number of sampled data in local search neighborhood
# #This is important to retain nearby relevant samples for IK
# # You can test/use the following combinations:
# # max dist only
# # nmax only
# # maxdist and nmax
# # maxdist, nmax, nmin
# # max dist, omax (with omax normally equal to 2 or 3)
#
# maxdist =  1200   # maximum distance to look for sampled data (usually a bit further then the range of variogram model)
# nmin =  10      # minimum number of sampled data to calculate the estimate
# nmax =  30      #maximum number of sampled data to calculate the estimate
# omax =  30      # maximum number of sampled data per octant to calculate the estimate
#
#
# #General formula for IK (remove the search parameters as needed)
# ik_lead <- krige(input_dataframe$pb ~ 1, locations = input_dataframe,
#                  newdata = kgrid,
#                  model = vgm_exponential$variogram_model,
#                  maxdist = maxdist,nmin = nmin, nmax= nmax, omax=omax)
#
# ik_lead
#
# this_ik_lead <- as.list(ik_lead$var1.pred)
# ik_dataframe <- as.data.frame(this_ik_lead)
# spher_ik <- data.frame(`IK_est` = unlist(ik_dataframe, use.names = FALSE) )
#
# #IK estimates
#
# #If IK estimates negative or above 1, limit the predicted probabilities to the range:
# ik_lead$var1.pred <- pmin(1, ik_lead$var1.pred )
# ik_lead$var1.pred <- pmax(0, ik_lead$var1.pred)
#
# ik_lead
#
# #Probability map from Indicator Kriging
#
# ggplot() + geom_stars(data = ik_lead,
#                       aes(fill = var1.pred, x = x, y = y)) +
#   xlab("x") + ylab("y") +
#   coord_equal()+
#   scale_fill_viridis_c(option = "B", direction=-1,
#                        na.value = "transparent", name="Probability")+
#   ggtitle("Indicator Map for Pb Probability (res:50m)")


########################################################################################
########################################################################################


cv_results <- data.frame(exp_ok_cor = cor(exp_ok_results$krige_cv$observed, exp_ok_results$krige_cv$var1.pred),



)

########################################################################################
# Cross validation using LOOCV      exp_ok_results$krige_cv

plot(exp_ok_results$krige_cv$observed~exp_ok_results$krige_cv$var1.pred, main=" LOOCV OK Exponential",xlab="Sampled Pb (ppm) ", ylab="Predicted Pb (ppm)", pch=19)
abline(0,1,col="red")

cor(exp_ok_results$krige_cv$observed, exp_ok_results$krige_cv$var1.pred)



plot(exp_uk_results$krige_cv$observed~exp_uk_results$krige_cv$var1.pred, main=" LOOCV OK Exponential",xlab="Sampled Pb (ppm) ", ylab="Predicted Pb (ppm)", pch=19)
abline(0,1,col="red")

cor(exp_uk_results$krige_cv$observed, exp_uk_results$krige_cv$var1.pred)






########################################################################################
########################################################################################
#         Cross Validation
########################################################################################

#Leave-one-out cross validation (LOOCV)

cv_exp_ok_leave <- krige.cv(pb ~ 1, input_dataframe, model=vgm_exponential$variogram_model, nmax=32)

#cv_exp_ok_leaveB <- krige.cv(pb ~ 1, input_dataframe, model=vgm_exponential$variogram_model, nmax=32)
cv_exp_ok_leave


#K-fold cross-validation (K=10)

cv_exp_ok_kfold <- krige.cv(pb ~ 1, input_dataframe, model=vgm_exponential$variogram_model, nmax=32, nfold=10)
cv_exp_ok_kfold

# -----------------------------
#Comparison sampled dataset vs estimated values

plot(cv_exp_ok_leave$observed~cv_exp_ok_leave$var1.pred, main=" cv_exp_ok_leave Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")

cor(cv_exp_ok_leave$observed, cv_exp_ok_leave$var1.pred)


plot(cv_exp_ok_kfold$observed~cv_exp_ok_kfold$var1.pred, main=" cv_exp_ok_kfold Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")

cor(cv_exp_ok_kfold$observed, cv_exp_ok_kfold$var1.pred)


# -----------------------------
# Linear correlation using Lin's concordance coefficient
library(valmetrics) # Lin's concordance coefficient
# Values close to 1 indicate a strong concordance between predicted and observed values,
# Values near -1 indicate a strong discordance.
#Values close to 0 indicate no concordance.
#In a plot of predicted values versus observed values, an LCCC-value of 1 means that the all data points are on the 1.1-line.

lin_cor1<- lccc(cv_exp_ok_leave$observed,cv_exp_ok_leave$var1.pred)
lin_cor1

lin_cor2<- lccc(cv_exp_ok_kfold$observed,cv_exp_ok_kfold$var1.pred)
lin_cor2



#Comparing sampled dataset vs estimated values Using the histogram
hist_sampled <- ggplot(input_dataframe, aes(x = pb)) + geom_histogram(binwidth=25, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)",breaks=seq(0,750,100)) + ggtitle("Topsoil lead (Pb) - sampled")
hist_sampled

ok.df <- as.data.frame(exp_ok_results$krige)

hist_ok_omni <- ggplot(ok.df, aes(x = var1.pred)) + geom_histogram(binwidth=25, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)",breaks=seq(0,750,100)) + ggtitle("Topsoil lead (Pb) - Ok estimates omni")
hist_ok_omni



plot_cprob_sampled <- plot(ecdf(input_dataframe$pb),main="Pb sampled",xlab="Pb(ppm)",ylab="Cumulative Probability")
plot_cprob_omni <- plot(ecdf(ok.df$var1.pred),main="Pb OK omni",xlab="Pb(ppm)",ylab="Cumulative Probability")



# -----------------------------
#Cross-validation evaluation using the Residuals (estimated - observed)


#Residuals histogram (should be similar to a normal distribution)
hist_cv_ok_omni_leave <- ggplot(cv_exp_ok_leave, aes(x = residual)) + geom_histogram(binwidth=25, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("leave Residuals - Ok estimates omni")
hist_cv_ok_omni_leave

hist_cv_ok_omni_kfold <- ggplot(cv_exp_ok_kfold, aes(x = residual)) + geom_histogram(binwidth=25, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("kfold Residuals - Ok estimates dir")
hist_cv_ok_omni_kfold

# Plotting the residuals
residuals_cv_ok_leave_plot <- ggplot(data = cv_exp_ok_leave) + geom_sf(aes(color = residual), size=2) + scale_color_viridis_c(option = "H", name="residuals")
residuals_cv_ok_leave_plot

residuals_cv_ok_kfold_plot <- ggplot(data = cv_exp_ok_kfold) + geom_sf(aes(color = residual), size=2) + scale_color_viridis_c(option = "H", name="residuals")
residuals_cv_ok_kfold_plot



# -----------------------------
#Cross Validation statistics


#Bias or Mean error should be 0
#Comparison with the median helps assess magnitude of the errors

M1<- mean(cv_exp_ok_leave$residual)
#M2<- mean(cv_ok_ani$residual)

m1<- median(cv_exp_ok_leave$residual)
#m2<- median(cv_ok_ani$residual)

s1<- sd(cv_exp_ok_leave$residual)
#s2<- sd(cv_ok_ani$residual)

#Root mean squared error (RMSE) should be low

R1<- sqrt(mean(cv_exp_ok_leave$residual^2))
#R2<- sqrt(mean(cv_ok_ani$residual^2))
#rmse_


#Mean Squared Deviation Ratio (MSDR) should be 1

MR1<- mean((cv_exp_ok_leave$residual^2/cv_exp_ok_leave$var1.var))
#MR2<- mean((cv_ok_ani$residual^2/cv_ok_ani$var1.var))


########################################################################################
#         Uncertainty
########################################################################################

# cv_results <- data.frame(ok_sph_cor = cor(sph_ok_results$krige_cv$observed, sph_ok_results$krige_cv$var1.pred),
#                          ok_exp_cor = cor(exp_ok_results$krige_cv$observed, exp_ok_results$krige_cv$var1.pred),
#                          ok_gau_cor = cor(gau_ok_results$krige_cv$observed, gau_ok_results$krige_cv$var1.pred),
#                          sk_sph_cor = cor(sph_sk_results$krige_cv$observed, sph_sk_results$krige_cv$var1.pred),
#                          sk_exp_cor = cor(exp_sk_results$krige_cv$observed, exp_sk_results$krige_cv$var1.pred),
#                          sk_gau_cor = cor(gau_sk_results$krige_cv$observed, gau_sk_results$krige_cv$var1.pred),
#                          uk_sph_cor = cor(sph_uk_results$krige_cv$observed, sph_uk_results$krige_cv$var1.pred),
#                          uk_exp_cor = cor(exp_uk_results$krige_cv$observed, exp_uk_results$krige_cv$var1.pred),
#                          uk_gau_cor = cor(gau_uk_results$krige_cv$observed, gau_uk_results$krige_cv$var1.pred),
#
#                          ok_sph_lin = lccc(sph_ok_results$krige_cv$observed, sph_ok_results$krige_cv$var1.pred),
#                          ok_exp_lin = lccc(exp_ok_results$krige_cv$observed, exp_ok_results$krige_cv$var1.pred),
#                          ok_gau_lin = lccc(gau_ok_results$krige_cv$observed, gau_ok_results$krige_cv$var1.pred),
#                          sk_sph_lin = lccc(sph_sk_results$krige_cv$observed, sph_sk_results$krige_cv$var1.pred),
#                          sk_exp_lin = lccc(exp_sk_results$krige_cv$observed, exp_sk_results$krige_cv$var1.pred),
#                          sk_gau_lin = lccc(gau_sk_results$krige_cv$observed, gau_sk_results$krige_cv$var1.pred),
#                          sk_sph_lin = lccc(sph_uk_results$krige_cv$observed, sph_uk_results$krige_cv$var1.pred),
#                          sk_exp_lin = lccc(exp_uk_results$krige_cv$observed, exp_uk_results$krige_cv$var1.pred),
#                          sk_gau_lin = lccc(gau_uk_results$krige_cv$observed, gau_uk_results$krige_cv$var1.pred)
# )


cv_results222 <- data.frame( `cv_header` <- data.frame("cor","lccc","mean","median","sd","rmse","msdr"),
                             `ok_sph_cv` <- data.frame(cor(sph_ok_results$krige_cv$observed, sph_ok_results$krige_cv$var1.pred),
                                                   lccc(sph_ok_results$krige_cv$observed, sph_ok_results$krige_cv$var1.pred),
                                                   mean(sph_ok_results$krige_cv$residual),
                                                   median(sph_ok_results$krige_cv$residual),
                                                   sd(sph_ok_results$krige_cv$residual),
                                                   sqrt(mean(sph_ok_results$krige_cv$residual^2)),
                                                   mean((sph_ok_results$krige_cv$residual^2/sph_ok_results$krige_cv$var1.var))),
                           `ok_exp_cv` <- data.frame( cor(exp_ok_results$krige_cv$observed, exp_ok_results$krige_cv$var1.pred),
                                                   lccc(exp_ok_results$krige_cv$observed, exp_ok_results$krige_cv$var1.pred),
                                                   mean(exp_ok_results$krige_cv$residual),
                                                   median(exp_ok_results$krige_cv$residual),
                                                    sd(exp_ok_results$krige_cv$residual),
                                                    sqrt(mean(exp_ok_results$krige_cv$residual^2)),
                                                    mean((exp_ok_results$krige_cv$residual^2/exp_ok_results$krige_cv$var1.var))),
                           `ok_gau_cv` <- data.frame(cor(gau_ok_results$krige_cv$observed, gau_ok_results$krige_cv$var1.pred),
                                                    lccc(gau_ok_results$krige_cv$observed, gau_ok_results$krige_cv$var1.pred),
                                                    mean(gau_ok_results$krige_cv$residual),
                                                   median(gau_ok_results$krige_cv$residual),
                                                   sd(gau_ok_results$krige_cv$residual),
                                                    sqrt(mean(gau_ok_results$krige_cv$residual^2)),
                                                  mean((gau_ok_results$krige_cv$residual^2/gau_ok_results$krige_cv$var1.var)))
)



cv_results2 <- data.frame( ok_sph_cv <- data.frame(ok_sph_cor = cor(sph_ok_results$krige_cv$observed, sph_ok_results$krige_cv$var1.pred),
                                                   ok_sph_lin_cor = lccc(sph_ok_results$krige_cv$observed, sph_ok_results$krige_cv$var1.pred),
                                                   ok_sph_mean = mean(sph_ok_results$krige_cv$residual),
                                                   ok_sph_median = median(sph_ok_results$krige_cv$residual),
                                                   ok_sph_sd = sd(sph_ok_results$krige_cv$residual),
                                                   ok_sph_rmse = sqrt(mean(sph_ok_results$krige_cv$residual^2)),
                                                   ok_sph_msdr = mean((sph_ok_results$krige_cv$residual^2/sph_ok_results$krige_cv$var1.var))),
                  ok_exp_cv <- data.frame(ok_exp_cor = cor(exp_ok_results$krige_cv$observed, exp_ok_results$krige_cv$var1.pred),
                                          ok_exp_lin_cor = lccc(exp_ok_results$krige_cv$observed, exp_ok_results$krige_cv$var1.pred),
                                          ok_exp_mean = mean(exp_ok_results$krige_cv$residual),
                                          ok_exp_median = median(exp_ok_results$krige_cv$residual),
                                          ok_exp_sd = sd(exp_ok_results$krige_cv$residual),
                                          ok_exp_rmse = sqrt(mean(exp_ok_results$krige_cv$residual^2)),
                                          ok_exp_msdr = mean((exp_ok_results$krige_cv$residual^2/exp_ok_results$krige_cv$var1.var))),
                  ok_gau_cv <- data.frame(ok_gau_cor = cor(gau_ok_results$krige_cv$observed, gau_ok_results$krige_cv$var1.pred),
                                          ok_gau_lin_cor = lccc(gau_ok_results$krige_cv$observed, gau_ok_results$krige_cv$var1.pred),
                                          ok_gau_mean = mean(gau_ok_results$krige_cv$residual),
                                          ok_gau_median = median(gau_ok_results$krige_cv$residual),
                                          ok_gau_sd = sd(gau_ok_results$krige_cv$residual),
                                          ok_gau_rmse = sqrt(mean(gau_ok_results$krige_cv$residual^2)),
                                          ok_gau_msdr = mean((gau_ok_results$krige_cv$residual^2/gau_ok_results$krige_cv$var1.var))),

                  sk_sph_cv <- data.frame(sk_sph_cor = cor(sph_sk_results$krige_cv$observed, sph_sk_results$krige_cv$var1.pred),
                                          sk_sph_lin_cor = lccc(sph_sk_results$krige_cv$observed, sph_sk_results$krige_cv$var1.pred),
                                          sk_sph_mean = mean(sph_sk_results$krige_cv$residual),
                                          sk_sph_median = median(sph_sk_results$krige_cv$residual),
                                          sk_sph_sd = sd(sph_sk_results$krige_cv$residual),
                                          sk_sph_rmse = sqrt(mean(sph_sk_results$krige_cv$residual^2)),
                                          sk_sph_msdr = mean((sph_sk_results$krige_cv$residual^2/sph_sk_results$krige_cv$var1.var))),
                  sk_exp_cv <- data.frame(sk_exp_ccor = cor(exp_sk_results$krige_cv$observed, exp_sk_results$krige_cv$var1.pred),
                                          sk_exp_clin_cor = lccc(exp_sk_results$krige_cv$observed, exp_sk_results$krige_cv$var1.pred),
                                          sk_exp_mean = mean(exp_sk_results$krige_cv$residual),
                                          sk_exp_cmedian = median(exp_sk_results$krige_cv$residual),
                                          sk_exp_csd = sd(exp_sk_results$krige_cv$residual),
                                          sk_exp_crmse = sqrt(mean(exp_sk_results$krige_cv$residual^2)),
                                          sk_exp_cmsdr = mean((exp_sk_results$krige_cv$residual^2/exp_sk_results$krige_cv$var1.var))),
                  sk_gau_cv <- data.frame(sk_gau_cor = cor(gau_sk_results$krige_cv$observed, gau_sk_results$krige_cv$var1.pred),
                                          sk_gau_lin_cor = lccc(gau_sk_results$krige_cv$observed, gau_sk_results$krige_cv$var1.pred),
                                          sk_gau_mean = mean(gau_sk_results$krige_cv$residual),
                                          sk_gau_median = median(gau_sk_results$krige_cv$residual),
                                          sk_gau_sd = sd(gau_sk_results$krige_cv$residual),
                                          sk_gau_rmse = sqrt(mean(gau_sk_results$krige_cv$residual^2)),
                                          sk_gau_msdr = mean((gau_sk_results$krige_cv$residual^2/gau_sk_results$krige_cv$var1.var))),

                  uk_sph_cv <- data.frame(uk_sph_cor = cor(sph_uk_results$krige_cv$observed, sph_uk_results$krige_cv$var1.pred),
                                          uk_sph_lin_cor = lccc(sph_uk_results$krige_cv$observed, sph_uk_results$krige_cv$var1.pred),
                                          uk_sph_mean = mean(sph_uk_results$krige_cv$residual),
                                          uk_sph_median = median(sph_uk_results$krige_cv$residual),
                                          uk_sph_sd = sd(sph_uk_results$krige_cv$residual),
                                          uk_sph_rmse = sqrt(mean(sph_uk_results$krige_cv$residual^2)),
                                          uk_sph_msdr = mean((sph_uk_results$krige_cv$residual^2/sph_uk_results$krige_cv$var1.var))),
                  uk_exp_cv <- data.frame(uk_exp_cor = cor(exp_uk_results$krige_cv$observed, exp_uk_results$krige_cv$var1.pred),
                                          uk_exp_lin_cor = lccc(exp_uk_results$krige_cv$observed, exp_uk_results$krige_cv$var1.pred),
                                          uk_exp_mean = mean(exp_uk_results$krige_cv$residual),
                                          uk_exp_median = median(exp_uk_results$krige_cv$residual),
                                          uk_exp_sd = sd(exp_uk_results$krige_cv$residual),
                                          uk_exp_rmse = sqrt(mean(exp_uk_results$krige_cv$residual^2)),
                                          uk_exp_msdr = mean((exp_uk_results$krige_cv$residual^2/exp_uk_results$krige_cv$var1.var))),
                  uk_gau_cv <- data.frame(uk_gau_ccor = cor(gau_uk_results$krige_cv$observed, gau_uk_results$krige_cv$var1.pred),
                                          uk_gau_clin_cor = lccc(gau_uk_results$krige_cv$observed, gau_uk_results$krige_cv$var1.pred),
                                          uk_gau_cmean = mean(gau_uk_results$krige_cv$residual),
                                          uk_gau_cmedian = median(gau_uk_results$krige_cv$residual),
                                          uk_gau_csd = sd(gau_uk_results$krige_cv$residual),
                                          uk_gau_crmse = sqrt(mean(gau_uk_results$krige_cv$residual^2)),
                                          uk_gau_cmsdr = mean((gau_uk_results$krige_cv$residual^2/gau_uk_results$krige_cv$var1.var)))
                  )

write.csv(cv_results2, "results/cross_validation_results.csv", row.names=TRUE, append=FALSE)



plot(sph_ok_results$krige_cv$observed~sph_ok_results$krige_cv$var1.pred, main=" sph_ok_results$krige_cv Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")

plot(exp_ok_results$krige_cv$observed~exp_ok_results$krige_cv$var1.pred, main=" exp_ok_results$krige_cv Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")

plot(gau_ok_results$krige_cv$observed~gau_ok_results$krige_cv$var1.pred, main=" gau_ok_results$krige_cv Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")



plot(sph_sk_results$krige_cv$observed~sph_sk_results$krige_cv$var1.pred, main=" sph_sk_results$krige_cv Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")

plot(exp_sk_results$krige_cv$observed~exp_sk_results$krige_cv$var1.pred, main=" exp_sk_results$krige_cv Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")

plot(gau_sk_results$krige_cv$observed~gau_sk_results$krige_cv$var1.pred, main=" gau_sk_results$krige_cv Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")



plot(sph_uk_results$krige_cv$observed~sph_uk_results$krige_cv$var1.pred, main=" sph_uk_results$krige_cv Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")

plot(exp_uk_results$krige_cv$observed~exp_uk_results$krige_cv$var1.pred, main=" exp_uk_results$krige_cv Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")

plot(gau_uk_results$krige_cv$observed~gau_uk_results$krige_cv$var1.pred, main=" gau_uk_results$krige_cv Scatterplot Example",xlab="True (sampled) value ", ylab="Predicted ", pch=19)
abline(0,1,col="red")



sph_ok_results$krige_cv$observed~sph_ok_results$krige_cv$var1.pred


# ----------  PREPARE RESIDUALS RESULTS

results_residuals <- NULL

#Residuals histogram (should be similar to a normal distribution)
results_residuals$sph_ok_hist <- ggplot(sph_ok_results$krige_cv, aes(x = residual)) + geom_histogram(binwidth=5, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Spherical OK Residuals")

results_residuals$exp_ok_hist <- ggplot(exp_ok_results$krige_cv, aes(x = residual)) + geom_histogram(binwidth=5, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Exponential OK Residuals")

results_residuals$gau_ok_hist <- ggplot(gau_ok_results$krige_cv, aes(x = residual)) + geom_histogram(binwidth=5, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Gaussian OK Residuals")


results_residuals$sph_sk_hist <- ggplot(sph_sk_results$krige_cv, aes(x = residual)) + geom_histogram(binwidth=5, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Spherical SK Residuals")

results_residuals$exp_sk_hist <- ggplot(exp_sk_results$krige_cv, aes(x = residual)) + geom_histogram(binwidth=5, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Exponential SK Residuals")

results_residuals$gau_sk_hist <- ggplot(gau_sk_results$krige_cv, aes(x = residual)) + geom_histogram(binwidth=5, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Gaussian SK Residuals")


results_residuals$sph_uk_hist <- ggplot(sph_uk_results$krige_cv, aes(x = residual)) + geom_histogram(binwidth=5, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Spherical UK Residuals")

results_residuals$exp_uk_hist <- ggplot(exp_uk_results$krige_cv, aes(x = residual)) + geom_histogram(binwidth=5, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Exponential UK Residuals")

results_residuals$gau_uk_hist <- ggplot(gau_uk_results$krige_cv, aes(x = residual)) + geom_histogram(binwidth=5, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Gaussian UK Residuals")


grid.arrange(results_residuals$sph_ok_hist,
             results_residuals$sph_sk_hist,
             results_residuals$sph_uk_hist,
             results_residuals$exp_ok_hist,
             results_residuals$exp_sk_hist,
             results_residuals$exp_uk_hist,
             results_residuals$gau_ok_hist,
             results_residuals$gau_sk_hist,
             results_residuals$gau_uk_hist,
             nrow=3,ncol=3 )


res_limits <- c(-47,98)
# Plotting the residuals
results_residuals$sph_ok_plot <- ggplot(data = sph_ok_results$krige_cv) + geom_sf(aes(color = residual, size=residual)) + scale_color_viridis_c(option = "H", name="residuals", limits=res_limits) + ggtitle("Spherical OK Residuals")

results_residuals$exp_ok_plot <- ggplot(data = exp_ok_results$krige_cv) + geom_sf(aes(color = residual, size=residual)) + scale_color_viridis_c(option = "H", name="residuals", limits=res_limits) + ggtitle("Exponential OK Residuals")

results_residuals$gau_ok_plot <- ggplot(data = gau_ok_results$krige_cv) + geom_sf(aes(color = residual, size=residual))+ scale_color_viridis_c(option = "H", name="residuals", limits=res_limits) + ggtitle("Gaussian OK Residuals")


results_residuals$sph_sk_plot <- ggplot(data = sph_sk_results$krige_cv) + geom_sf(aes(color = residual, size=residual)) + scale_color_viridis_c(option = "H", name="residuals", limits=res_limits) + ggtitle("Spherical SK Residuals")

results_residuals$exp_sk_plot <- ggplot(data = exp_sk_results$krige_cv) + geom_sf(aes(color = residual, size=residual)) + scale_color_viridis_c(option = "H", name="residuals", limits=res_limits) + ggtitle("Exponential SK Residuals")

results_residuals$gau_sk_plot <- ggplot(data = gau_sk_results$krige_cv) + geom_sf(aes(color = residual, size=residual))+ scale_color_viridis_c(option = "H", name="residuals", limits=res_limits) + ggtitle("Gaussian SK Residuals")


results_residuals$sph_uk_plot <- ggplot(data = sph_uk_results$krige_cv) + geom_sf(aes(color = residual, size=residual)) + scale_color_viridis_c(option = "H", name="residuals", limits=res_limits) + ggtitle("Spherical UK Residuals")

results_residuals$exp_uk_plot <- ggplot(data = exp_uk_results$krige_cv) + geom_sf(aes(color = residual, size=residual)) + scale_color_viridis_c(option = "H", name="residuals", limits=res_limits) + ggtitle("Exponential UK Residuals")

results_residuals$gau_uk_plot <- ggplot(data = gau_uk_results$krige_cv) +geom_sf(aes(color = residual, size=residual)) + scale_color_viridis_c(option = "H", name="residuals", limits=res_limits) + ggtitle("Gaussian UK Residuals")


grid.arrange(results_residuals$sph_ok_plot,
             results_residuals$sph_sk_plot,
             results_residuals$sph_uk_plot,
             results_residuals$exp_ok_plot,
             results_residuals$exp_sk_plot,
             results_residuals$exp_uk_plot,
             results_residuals$gau_ok_plot,
             results_residuals$gau_sk_plot,
             results_residuals$gau_uk_plot,
             nrow=3,ncol=3 )

########################################################################################
#.  Estimation results, statistics and histograms

numSummary(input_raw[,"pb", drop=FALSE], statistics=c("mean", "sd", "se(mean)", "var", "IQR", "quantiles", "CV", "skewness", "kurtosis"), quantiles=c(0,.25,.5,.75,1), type="2")



# ----------  COLLECT ESTIMATION RESULTS

results_estimation <- c(data.frame(`OK_SPH` = unlist((as.data.frame(sph_ok_results$krige$var1.pred)), use.names = FALSE)),
                        data.frame(`OK_EXP` = unlist((as.data.frame(exp_ok_results$krige$var1.pred)), use.names = FALSE)),
                        data.frame(`OK_GAU` = unlist((as.data.frame(gau_ok_results$krige$var1.pred)), use.names = FALSE)),
                        data.frame(`SK_SPH` = unlist((as.data.frame(sph_sk_results$krige$var1.pred)), use.names = FALSE)),
                        data.frame(`SK_EXP` = unlist((as.data.frame(exp_sk_results$krige$var1.pred)), use.names = FALSE)),
                        data.frame(`SK_GAU` = unlist((as.data.frame(gau_sk_results$krige$var1.pred)), use.names = FALSE)),
                        data.frame(`UK_SPH` = unlist((as.data.frame(sph_uk_results$krige$var1.pred)), use.names = FALSE)),
                        data.frame(`UK_EXP` = unlist((as.data.frame(exp_uk_results$krige$var1.pred)), use.names = FALSE)),
                        data.frame(`UK_GAU` = unlist((as.data.frame(gau_uk_results$krige$var1.pred)), use.names = FALSE)) )

#summary.data.frame(results_estimation)
write.csv(results_estimation, "results/results_estimation.csv", row.names=TRUE, append=FALSE)
#write.csv(summary.data.frame(results_estimation), "results/results_estimation_stats.csv", append=FALSE)

results_estimation_summary <- numSummary(results_estimation, statistics=c("mean", "sd", "se(mean)", "var", "IQR", "quantiles", "CV", "skewness", "kurtosis"), quantiles=c(0,.025,.25,.5,.75,.975,1), type="2")
results_estimation_summary
write.csv(results_estimation_summary$table, "results/results_estimation_summary.csv", row.names=TRUE, append=FALSE)


#temp_est_summ <- summary_stats(results_estimation)


# ----------  COLLECT VARIANCE RESULTS

results_variance <- c(data.frame(`OK_SPH` = unlist((as.data.frame(sph_ok_results$krige$var1.var)), use.names = TRUE)),
                        data.frame(`OK_EXP` = unlist((as.data.frame(exp_ok_results$krige$var1.var)), use.names = TRUE)),
                        data.frame(`OK_GAU` = unlist((as.data.frame(gau_ok_results$krige$var1.var)), use.names = TRUE)),
                        data.frame(`SK_SPH` = unlist((as.data.frame(sph_sk_results$krige$var1.var)), use.names = TRUE)),
                        data.frame(`SK_EXP` = unlist((as.data.frame(exp_sk_results$krige$var1.var)), use.names = TRUE)),
                        data.frame(`SK_GAU` = unlist((as.data.frame(gau_sk_results$krige$var1.var)), use.names = TRUE)),
                        data.frame(`UK_SPH` = unlist((as.data.frame(sph_uk_results$krige$var1.var)), use.names = TRUE)),
                        data.frame(`UK_EXP` = unlist((as.data.frame(exp_uk_results$krige$var1.var)), use.names = TRUE)),
                        data.frame(`UK_GAU` = unlist((as.data.frame(gau_uk_results$krige$var1.var)), use.names = TRUE)) )

#summary.data.frame(results_variance)
write.csv(results_variance, "results/results_variance.csv", row.names=TRUE, append=FALSE)
#write.csv(summary.data.frame(results_variance), "results/results_variance_stats.csv", append=FALSE)


results_variance_summary <- numSummary(results_variance, statistics=c("mean", "sd", "se(mean)", "var", "IQR", "quantiles", "CV", "skewness", "kurtosis"), quantiles=c(0,.025,.25,.5,.75,.975,1), type="2")
results_variance_summary
write.csv(results_variance_summary$table, "results/results_variation_summary.csv", row.names=TRUE, append=FALSE)


results_variance_hist <- NULL
results_variance_hist$sph_ok_hist <- ggplot(sph_ok_results$krige_cv, aes(x = var1.var)) + geom_histogram(binwidth=25, alpha=.5, position="identity") + scale_x_continuous(name="Pb (ppm)") + ggtitle("Spherical OK Residuals")


results_variance_hist$sph_ok_hist
# -----


results_cv_variance <- c(data.frame(`OK_SPH` = unlist((as.data.frame(sph_ok_results$krige_cv$var1.var)), use.names = TRUE)),
                      data.frame(`OK_EXP` = unlist((as.data.frame(exp_ok_results$krige_cv$var1.var)), use.names = TRUE)),
                      data.frame(`OK_GAU` = unlist((as.data.frame(gau_ok_results$krige_cv$var1.var)), use.names = TRUE)),
                      data.frame(`SK_SPH` = unlist((as.data.frame(sph_sk_results$krige_cv$var1.var)), use.names = TRUE)),
                      data.frame(`SK_EXP` = unlist((as.data.frame(exp_sk_results$krige_cv$var1.var)), use.names = TRUE)),
                      data.frame(`SK_GAU` = unlist((as.data.frame(gau_sk_results$krige_cv$var1.var)), use.names = TRUE)),
                      data.frame(`UK_SPH` = unlist((as.data.frame(sph_uk_results$krige_cv$var1.var)), use.names = TRUE)),
                      data.frame(`UK_EXP` = unlist((as.data.frame(exp_uk_results$krige_cv$var1.var)), use.names = TRUE)),
                      data.frame(`UK_GAU` = unlist((as.data.frame(gau_uk_results$krige_cv$var1.var)), use.names = TRUE)) )

results_cv_variance_summary <- numSummary(results_cv_variance, statistics=c("mean", "sd", "se(mean)", "var", "IQR", "quantiles", "CV", "skewness", "kurtosis"), quantiles=c(0,.025,.25,.5,.75,.975,1), type="2")
results_cv_variance_summary
write.csv(results_cv_variance_summary$table, "results/results_cv_variation_summary.csv", row.names=TRUE, append=FALSE)


# ----------  COLLECT RESIDUALS RESULTS

results_residuals <- c(data.frame(`OK_SPH` = unlist((as.data.frame(sph_ok_results$krige_cv$residual)), use.names = TRUE)),
                      data.frame(`OK_EXP` = unlist((as.data.frame(exp_ok_results$krige_cv$residual)), use.names = TRUE)),
                      data.frame(`OK_GAU` = unlist((as.data.frame(gau_ok_results$krige_cv$residual)), use.names = TRUE)),
                      data.frame(`SK_SPH` = unlist((as.data.frame(sph_sk_results$krige_cv$residual)), use.names = TRUE)),
                      data.frame(`SK_EXP` = unlist((as.data.frame(exp_sk_results$krige_cv$residual)), use.names = TRUE)),
                      data.frame(`SK_GAU` = unlist((as.data.frame(gau_sk_results$krige_cv$residual)), use.names = TRUE)),
                      data.frame(`UK_SPH` = unlist((as.data.frame(sph_uk_results$krige_cv$residual)), use.names = TRUE)),
                      data.frame(`UK_EXP` = unlist((as.data.frame(exp_uk_results$krige_cv$residual)), use.names = TRUE)),
                      data.frame(`UK_GAU` = unlist((as.data.frame(gau_uk_results$krige_cv$residual)), use.names = TRUE)) )

#summary.data.frame(results_residuals)
write.csv(results_residuals, "results/results_residuals.csv", row.names=TRUE, append=FALSE)
#write.csv(summary.data.frame(results_residuals), "results/results_residuals_stats.csv", append=FALSE)

results_residuals_summary <- numSummary(results_residuals, statistics=c("mean", "sd", "se(mean)", "var", "IQR", "quantiles", "CV", "skewness", "kurtosis"), quantiles=c(0,.025,.25,.5,.75,.975,1), type="2")
results_residuals_summary
write.csv(results_residuals_summary$table, "results/results_residuals_summary.csv", row.names=TRUE, append=FALSE)



# ----------  GENERATE IMAGE GRIDS OF ESTIMATION MAPS FOR EACH MODEL

grid.arrange(vgm_spherical$variogram_model_plot, vgm_exponential$variogram_model_plot, vgm_gaussian$variogram_model_plot, nrow=1)

grid.arrange(sph_ok_results$estMap, exp_ok_results$estMap, gau_ok_results$estMap, nrow = 1)

grid.arrange(sph_sk_results$estMap, exp_sk_results$estMap, gau_sk_results$estMap, nrow = 1)

grid.arrange(sph_uk_results$estMap,exp_uk_results$estMap, gau_uk_results$estMap, nrow=1 )

grid.arrange(sph_ok_results$estMap, sph_sk_results$estMap, sph_uk_results$estMap, exp_ok_results$estMap,exp_sk_results$estMap,exp_uk_results$estMap, nrow=2,ncol=3)

# variance map grid
grid.arrange(sph_ok_results$varMap, sph_sk_results$varMap, sph_uk_results$varMap, exp_ok_results$varMap,exp_sk_results$varMap,exp_uk_results$varMap, nrow=2,ncol=3)




grid.arrange(sph_ok_results$estMap, sph_sk_results$estMap, sph_uk_results$estMap, exp_ok_results$estMap,exp_sk_results$estMap,exp_uk_results$estMap, gau_ok_results$estMap, gau_sk_results$estMap, gau_uk_results$estMap,  nrow=3,ncol=3)



# variance map grid
grid.arrange(sph_ok_results$varMap, sph_sk_results$varMap, sph_uk_results$varMap, exp_ok_results$varMap,exp_sk_results$varMap,exp_uk_results$varMap, gau_ok_results$varMap, gau_sk_results$varMap, gau_uk_results$varMap, nrow=3,ncol=3)



# indicator map grid
grid.arrange(sph_ok_results$indMap, sph_sk_results$indMap, sph_uk_results$indMap, exp_ok_results$indMap,exp_sk_results$indMap,exp_uk_results$indMap, gau_ok_results$indMap, gau_sk_results$indMap, gau_uk_results$indMap, nrow=3,ncol=3)

grid.arrange(sph_ok_results$indMap, sph_sk_results$indMap, sph_uk_results$indMap, exp_ok_results$indMap,exp_sk_results$indMap,exp_uk_results$indMap, nrow=2,ncol=3)


# ------------  GENERATE HISTOGRAMS FOR ESTIMATION RESULTS

results_histogram <- NULL

results_histogram$sampled_histogram <- with(input_dataframe, hist(pb, freq=TRUE, breaks=24, col="darkgray", xlim=hist_range,ylim=c(0,6),xlab=attribtue_label, main="Sampled Data"))

# Ordinary Kriging histograms
results_histogram$sph_ok_histogram <- with(results_estimation, hist(OK_SPH, freq=TRUE, breaks=12, col="darkgray", xlim=hist_range,xlab=attribtue_label,main="Spherical OK Estimate"))
results_histogram$exp_ok_histogram <- with(results_estimation, hist(OK_EXP, freq=TRUE, breaks=12, col="darkgray", xlim=hist_range,xlab=attribtue_label,main="Exponential OK Estimate"))
results_histogram$gau_ok_histogram <- with(results_estimation, hist(OK_GAU, freq=TRUE, breaks=12, col="darkgray", xlim=hist_range,xlab=attribtue_label,main="Gaussian OK Estimate"))


# Simple Kriging histograms
results_histogram$sph_sk_histogram <- with(results_estimation, hist(SK_SPH, freq=TRUE, breaks=12, col="darkgray", xlim=hist_range,xlab=attribtue_label,main="Spherical SK Estimate"))
results_histogram$exp_sk_histogram <- with(results_estimation, hist(SK_EXP, freq=TRUE, breaks=12, col="darkgray", xlim=hist_range,xlab=attribtue_label,main="Exponential SK Estimate"))
results_histogram$gau_sk_histogram <- with(results_estimation, hist(SK_GAU, freq=TRUE, breaks=12, col="darkgray", xlim=hist_range,xlab=attribtue_label,main="Gaussian SK Estimate"))


# Universal Kriging histograms
results_histogram$sph_uk_histogram <- with(results_estimation, hist(UK_SPH, freq=TRUE, breaks=12, col="darkgray", xlim=hist_range,xlab=attribtue_label,main="Spherical UK Estimate"))
results_histogram$exp_uk_histogram <- with(results_estimation, hist(UK_EXP, freq=TRUE, breaks=12, col="darkgray", xlim=hist_range,xlab=attribtue_label,main="Exponential UK Estimate"))
results_histogram$gau_uk_histogram <- with(results_estimation, hist(UK_GAU, freq=TRUE, breaks=12, col="darkgray", xlim=hist_range,xlab=attribtue_label,main="Gaussian UK Estimate"))




# ------------  GENERATE HISTOGRAMS FOR VARIANCE RESULTS
var_hist_range <- c(157,1191)

results_variance_histogram <- NULL

results_variance_histogram$sampled_histogram <- with(input_dataframe, hist(pb, freq=TRUE, breaks=24, col="darkgray", xlim=hist_range,ylim=c(0,6),xlab=attribtue_label, main="Sampled Data"))

# Ordinary Kriging histograms
results_variance_histogram$sph_ok_histogram <- with(results_variance, hist(OK_SPH, freq=TRUE, breaks=12, col="darkgray", xlim=var_hist_range,xlab=attribtue_label,main="Spherical OK Variance"))
results_variance_histogram$exp_ok_histogram <- with(results_variance, hist(OK_EXP, freq=TRUE, breaks=12, col="darkgray", xlim=var_hist_range,xlab=attribtue_label,main="Exponential OK Variance"))
results_variance_histogram$gau_ok_histogram <- with(results_variance, hist(OK_GAU, freq=TRUE, breaks=12, col="darkgray", xlim=var_hist_range,xlab=attribtue_label,main="Gaussian OK Variance"))


# Simple Kriging histograms
results_variance_histogram$sph_sk_histogram <- with(results_variance, hist(SK_SPH, freq=TRUE, breaks=12, col="darkgray", xlim=var_hist_range,xlab=attribtue_label,main="Spherical SK Variance"))
results_variance_histogram$exp_sk_histogram <- with(results_variance, hist(SK_EXP, freq=TRUE, breaks=12, col="darkgray", xlim=var_hist_range,xlab=attribtue_label,main="Exponential SK Variance"))
results_variance_histogram$gau_sk_histogram <- with(results_variance, hist(SK_GAU, freq=TRUE, breaks=6, col="darkgray", xlim=var_hist_range,xlab=attribtue_label,main="Gaussian SK Variance"))


# Universal Kriging histograms
results_variance_histogram$sph_uk_histogram <- with(results_variance, hist(UK_SPH, freq=TRUE, breaks=12, col="darkgray", xlim=var_hist_range,xlab=attribtue_label,main="Spherical UK Variance"))
results_variance_histogram$exp_uk_histogram <- with(results_variance, hist(UK_EXP, freq=TRUE, breaks=12, col="darkgray", xlim=var_hist_range,xlab=attribtue_label,main="Exponential UK Variance"))
results_variance_histogram$gau_uk_histogram <- with(results_variance, hist(UK_GAU, freq=TRUE, breaks=12, col="darkgray", xlim=var_hist_range,xlab=attribtue_label,main="Gaussian UK Variance"))


# ########################################################################################
# #         Output Results
# ########################################################################################
# # Compile/Mosaic output maps
#
cprob_sampled <- plot(ecdf(input_raw$pb),main="Sampled Data",xlab="Pb(ppm)",ylab="Cumulative Probability", xlim=hist_range)

cprob_ok_sph <- plot(ecdf(estimation_results$OK_SPH),main="Spherical OK Estimate",xlab=attribtue_label,ylab="Cumulative Probability", xlim=hist_range)
cprob_ok_exp <- plot(ecdf(estimation_results$OK_EXP),main="Exponential OK Estimate",xlab=attribtue_label,ylab="Cumulative Probability", xlim=hist_range)

cprob_sk_sph <- plot(ecdf(estimation_results$SK_SPH),main="Spherical SK Estimate",xlab=attribtue_label,ylab="Cumulative Probability", xlim=hist_range)
cprob_sk_exp <- plot(ecdf(estimation_results$SK_EXP),main="Exponential SK Estimate",xlab=attribtue_label,ylab="Cumulative Probability", xlim=hist_range)

cprob_uk_sph <- plot(ecdf(uk_data_sph$krige$var1.pred),main="Spherical UK Estimate",xlab=attribtue_label,ylab="Cumulative Probability", xlim=hist_range)
cprob_uk_exp <- plot(ecdf(uk_data_exp$krige$var1.pred),main="Exponential UK Estimate",xlab=attribtue_label,ylab="Cumulative Probability", xlim=hist_range)



#
# this_ecdf <-  ggplot(ok_pred_results, aes(OK_est)) + stat_ecdf(geom = "step")
# this_ecdf
#
# this_ecdf2 <-  ggplot(ok_pred_results, aes(OK_est)) + stat_ecdf(geom = "step")
# this_ecdf2
#
# this_ecdf3 <-  ggplot(sk_pred_results, aes(OK_est)) + stat_ecdf(geom = "step")
# this_ecdf3
#
# ecdf(input_raw$pb)
#
# ggplot(df, aes(height)) + stat_ecdf(geom = "step")
#
# #grid.arrange(ok_results$estMap, ok_results$varMap, sk_results$estMap, sk_results$varMap, uk_results$estMap, uk_results$varMap,  ncol=2, nrow=3)
#
# cat("Estimate Min/Max\n",
#     "OK_est_min",
#     "OK_est_max",
#     "SK_est_min",
#     "SK_est_max",
#     "\n",
#     min(ok_results$krige$var1.pred),
#     max(ok_results$krige$var1.pred),
#     min(sk_results$krige$var1.pred),
#     max(sk_results$krige$var1.pred),
#     sep = "\t")
#
# cat("Variance Min/Max\n",
#     "OK_var_min",
#     "OK_var_max",
#     "SK_var_min",
#     "SK_var_max",
#     "\n",
#     min(ok_results$krige$var1.var),
#     max(ok_results$krige$var1.var),
#     min(sk_results$krige$var1.var),
#     max(sk_results$krige$var1.var),
#     sep = "\t")
#
#
#
#
# grid.arrange(ok_results$estMap, ok_results$varMap, sk_results$estMap, sk_results$varMap,  ncol=2, nrow=2)
#



