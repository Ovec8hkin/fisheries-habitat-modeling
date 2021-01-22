library(tidyverse)  
library(ggplot2)
library(dplyr)
library(gridExtra)

plot_gam_spline <- function(spline){
  p = plot(spline)+l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) + l_points(alpha=0.2) + geom_hline(yintercept=0)
  
  return(p)
}

plot_repredicted_distributions <- function(data){
  catch_plot = ggplot(data,aes(x=easting,y=northing, weight=catch))+geom_bin2d(bins=100)+ggtitle("Catch")
  predicted_plot = ggplot(data,aes(x=easting,y=northing, weight=re_predicted))+geom_bin2d(bins=100)+ggtitle("Predicted")
  diff_plot = ggplot(data,aes(x=easting,y=northing, weight=diff))+geom_bin2d(bins=100)+ggtitle("Diff")+scale_fill_gradient2()
  return(grid.arrange(catch_plot, predicted_plot, diff_plot, nrow=1))
}

plot_repredicted_params <- function(data, x, binwidth){
  # Temperature, Depth, and Salinity Plot Comparisons
  return (
    ggplot(data)+
      geom_histogram(aes(x=x, weight=catch), binwidth=binwidth, fill='black', alpha=0.5)+
      geom_histogram(aes(x=x, weight=re_predicted), binwidth=binwidth, fill='red', alpha=0.3)
  )
}

plot_cumulative_catch <- function(df, title){
  summed_predictions = df %>% group_by(lat, lon) %>% summarise(Total=sum(predicted))
  summed_predictions = as.data.frame(summed_predictions)
  p = ggplot(summed_predictions)+
    geom_point(
      aes(x=lon, 
          y=lat, 
          colour=log(floor(Total)), 
          alpha=ifelse(log(floor(Total)) < 0, 0, 1)
      ), 
      size=0.1
    )+
    scale_color_gradient(low="white", high="red")+
    labs(x="Longitude", 
         y="Latitude", 
         title=title, 
         color="Catch"
    )+
    guides(alpha=FALSE)
  return(p)
}

plot_average_catch <- function(df, title){
  average_predictions = df %>% group_by(lat, lon) %>% summarise(Average=mean(predicted))
  average_predictions = as.data.frame(average_predictions)
  p = ggplot(average_predictions)+
    geom_point(
      aes(x=lon, 
          y=lat, 
          colour=log(floor(Average)), 
          alpha=ifelse(log(floor(Average)) < 0, 0, 1)
      ), 
      size=0.1
    )+
    labs(x="Longitude",
         y="Latitude",
         title=title,
         colour="Catch"
    )+
    guides(alpha=FALSE)
  return(p)
}

plot_yearly_catch <- function(df, title){
  yearly_average_predictions = df %>% 
    group_by(lat, lon, year) %>% 
    summarise(Total=sum(predicted))
  yearly_average_predictions = as.data.frame(yearly_average_predictions)
  p = ggplot(yearly_average_predictions)+
    geom_point(
      aes(x=lon,
          y=lat,
          colour=log(floor(Total)),
          alpha=ifelse(log(floor(Total)) < 0, 0, 1)
      ),
      size=0.1
    )+
    scale_color_gradient2()+
    labs(x="Longitude",
         y="Latitude",
         title=title,
         color="Catch"
    )+
    guides(alpha=FALSE)+facet_wrap(~year)
  return(p)
}

plot_monthly_catch <- function(df, title){
  monthly_average_predictions = df %>% 
    group_by(lat, lon, month) %>% 
    summarise(Average=mean(predicted))
  
  p = ggplot(monthly_average_predictions)+
    geom_point(
      aes(x=lon, 
          y=lat, 
          colour=log(floor(Average)), 
          alpha=ifelse(log(floor(Average)) < 0, 0, 1)
      ), 
      size=0.1
    )+
    labs(x="Longitude",
         y="Latitude",
         title=title,
         colour="Catch"
    )+
    scale_colour_gradient(low="white", high="red")+
    guides(alpha=FALSE)+facet_wrap(~month)
  
  return(p)
}