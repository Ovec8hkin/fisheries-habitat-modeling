source("src/gam_functions.R")
source("src/data.R")
source("src/plot_functions.R")

# Load ecomon catch data from CSV files
cod_catch_data = get_ecomon_catch_data("gadmor_100m3", "auxdata/catch-data/ecomon_cod.csv")
had_catch_data = get_ecomon_catch_data("melaeg_100m3", "auxdata/catch-data/ecomon_haddock.csv")
flo_catch_data = get_ecomon_catch_data("limfer_100m3", "auxdata/catch-data/ecomon_flounder.csv")
mac_catch_data = get_ecomon_catch_data("scosco_100m3", "auxdata/catch-data/ecomon_mackerel.csv")
but_catch_data = get_ecomon_catch_data("pepspp_100m3", "auxdata/catch-data/ecomon_butterfish.csv")


# Fit global models as negative binomial GAMs of depth, bottom/surface temperature, and bottom/surface salinity
model.formula = catch ~ s(depth) +s(sfc_temp) + s(btm_temp) + s(sfc_salt) + s(btm_salt)
cod.global.model = mgcv::gam(model.formula, data=cod_catch_data, family="nb")
had.global.model = mgcv::gam(model.formula, data=had_catch_data, family="nb")
flo.global.model = mgcv::gam(model.formula, data=flo_catch_data, family="nb")
mac.global.model = mgcv::gam(model.formula, data=mac_catch_data, family="nb")
but.global.model = mgcv::gam(model.formula, data=but_catch_data, family="nb")

# Use MuMIN to to test all other possible models
options(na.action = "na.fail")
cod.all.models = dredge_models(cod.global.model)
had.all.models = dredge_models(had.global.model)
flo.all.models = dredge_models(flo.global.model)
mac.all.models = dredge_models(mac.global.model)
but.all.models = dredge_models(but.global.model)

# Format all model tables for viewing convenience
cod.model_table = format_model_table(cod.all.models)
had.model_table = format_model_table(had.all.models)
flo.model_table = format_model_table(flo.all.models)
mac.model_table = format_model_table(mac.all.models)
but.model_table = format_model_table(but.all.models)

# Select "Best Model" and "depth+sfc_temp+sfc_salt model" for further analysis 
cod.best_model = eval(getCall(cod.all.models, 1))
cod.simple_model = eval(getCall(cod.all.models, 4))

had.best_model = eval(getCall(had.all.models, 1))
had.simple_model = eval(getCall(had.all.models, 5))

flo.best_model = eval(getCall(flo.all.models, 1))
flo.simple_model = eval(getCall(flo.all.models, 7))

mac.best_model = eval(getCall(mac.all.models, 1))
mac.simple_model = eval(getCall(mac.all.models, 9))

but.best_model = eval(getCall(but.all.models, 1))
but.simple_model = eval(getCall(but.all.models, 7))

# Create list of models for use in predictive loops
best.models.list = list(cod.best_model, had.best_model, flo.best_model, mac.best_model, but.best_model)
simple.models.list = list(cod.simple_model, had.simple_model, flo.simple_model, mac.simple_model, but.simple_model)

# Save fitted GAM objects to RData file for use in predictive workflows
save(list=ls(.GlobalEnv), file="fitted_gams.Rdata")




species.list = c("Atlantic Cod", "Haddock", "Yellowtail Flounder", "Atlantic Mackerel", "American Butterfish")

compute_habitat_stats <- function(data){
  
  zero_cutoff = 0.01
  
  pos_data = subset(data, predicted > zero_cutoff)
  
  quantiles = quantile(pos_data$predicted, c(0.25, 0.75))
  
  habitat.high = subset(pos_data, predicted > quantiles[2])
  habitat.med = subset(pos_data, quantiles[1] <= predicted & predicted <= quantiles[2])
  habitat.low = subset(pos_data, predicted < quantiles[1] & predicted > zero_cutoff)
  
  areasum = sum(data$area_km2)
  
  habitat.high.areasum = sum(habitat.high$area_km2)
  habitat.med.areasum = sum(habitat.med$area_km2)
  habitat.low.areasum = sum(habitat.low$area_km2)
  habitat.none.areasum = areasum - habitat.high.areasum - habitat.med.areasum - habitat.low.areasum
  
  habitat.high.perc = habitat.high.areasum/areasum
  habitat.med.perc = habitat.med.areasum/areasum
  habitat.low.perc = habitat.low.areasum/areasum
  habitat.none.perc = habitat.none.areasum/areasum
  
  population = sum(pos_data$predicted*pos_data$area_km2)
  
  weighted.east = weighted.mean(data$east, data$predicted)
  weighted.north = weighted.mean(data$north, data$predicted)
  
  data = c(population,
           areasum,
           habitat.high.areasum, 
           habitat.med.areasum,
           habitat.low.areasum,
           habitat.none.areasum,
           habitat.high.perc,
           habitat.med.perc,
           habitat.low.perc,
           habitat.none.perc,
           weighted.east,
           weighted.north,
           quantiles[[1]],
           quantiles[[2]]
           )
  
  return(data)
  
}
load_environmental_data <- function(y, m, s){
  
  month = str_pad(m, 2, pad="0")
  fname=paste("auxdata/environmental-data/environmental_data_", y, month, ".csv", sep="")
  prediction.data = read.csv(fname)
  colnames(prediction.data) = c("east", "north", "sfc_temp", "btm_temp", "sfc_salt", "btm_salt", "depth", "area")
  prediction.data$area_km2 = prediction.data$area/1000/1000
  
  sub = (prediction.data$sfc_salt > 30 & prediction.data$sfc_salt < 35)
  
  edata = prediction.data[sub,]
  
  len = length(edata[,1])
  edata$year = rep(y, len)
  edata$month = rep(m, len)
  edata$species = rep(s, len)
  
  return(edata)
  
}

habitat.stats = data.frame()
predictions = data.frame()
for(i in 3:3){
  s = species.list[i]
  model = best.models.list[[i]]
  print(s)
  for(y in 1978:2016){
    print(y)
    for(m in 1:12){
      print(m)
      
      edata = load_environmental_data(y, m, s)
      
      pred_data = predict_data(model, edata)
      a = compute_habitat_stats(pred_data)
      a = c(y, m, a, s)
      
      habitat.stats = rbind(habitat.stats, a)
      predictions = rbind(predictions, pred_data)
    }
  }
}

colnames(habitat.stats) <- c("year", 
                             "month",
                             "population",
                             "area", 
                             "high", 
                             "medium",
                             "low", 
                             "none",
                             "perc-h",
                             "perc-m",
                             "perc-l",
                             "perc-n",
                             "weighted.east",
                             "weighted.north",
                             "25%", 
                             "75%", 
                             "species")

habitat.stats[,1:16] = lapply(habitat.stats[,1:16], function(x) if(is.character(x)) as.numeric(x) else x)
habitat.stats[,3:16] = round(habitat.stats[,3:16], 3)
habitat.stats = habitat.stats[,c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,17,15,16)]


habitat.stats.long = reshape::melt(habitat.stats, 
                                   id.vars=c("year","species"), 
                                   measure.vars=c("high", "medium", "low", "none", 
                                                  "perc-h", "perc-m", "perc-l", "perc-n", 
                                                  "population", "weighted.east", "weighted.north")
                                   )

population.stats = subset(habitat.stats.long, variable %in% c("population"))
habitat.percent.stats = subset(habitat.stats.long, variable %in% c("perc-h", "perc-m", "perc-l", "perc-n"))
habitat.area.stats = subset(habitat.stats.long, variable %in% c("high", "medium", "low", "none"))

weighted.north.stats = subset(habitat.stats.long, variable %in% c("weighted.north"))

annual_population = population.stats %>% 
  group_by(Year, Species) %>% 
  summarise(Pop=sum(value))

annual_percentage = habitat.percent.stats %>% 
  group_by(Year, Species, variable) %>% 
  summarise(percent=mean(value))

annual_area = habitat.area.stats %>% 
  group_by(Year, Species, variable) %>% 
  summarise(area=sum(value))

annual_north = weighted.north.stats %>%
  group_by(year, species) %>%
  summarise(north=mean(value))

# Annual population plot (Line)
ggplot(annual_population, aes(x=Year, y=Pop, color=Species))+geom_point()+geom_line()+ylim(0, 5000000)

# Annual available area (km2) by type (Line)
ggplot(annual_area, aes(x=Year, y=area, color=Species))+
  geom_point()+geom_line()+facet_wrap(~variable)

# Annual area type percentage (Stacked Bar)
ggplot(annual_percentage)+
  geom_col(aes(x=Year, y=percent, fill=variable))+facet_wrap(~Species)


ggplot(prediction.data)+
  geom_point(aes(x=east, y=north, color=predicted), size=0.001)
  

ggplot(annual_north)+geom_point(aes(x=year, y=north, color=species))+geom_line(aes(x=year, y=north, color=species))+
  geom_line(data=cod.lm.predic, aes(x=year, y=pred), color="red")+
  geom_line(data=had.lm.predic, aes(x=year, y=pred), color="green")
  
  
cod.lm.predic = data.frame(year=c(1978:2016))
had.lm.predic = data.frame(year=c(1978:2016))

cod_annual_north = subset(annual_north, species=="Atlantic Cod")
cod_north_lm = lm(north~year, data=cod_annual_north)
cod_north_lm_pred = predict(cod_north_lm, cod.lm.predic)
cod.lm.predic$pred = cod_north_lm_pred

had_annual_north = subset(annual_north, species=="Haddock")
had_north_lm = lm(north~year, data=had_annual_north)
had_north_lm_pred = predict(had_north_lm, had.lm.predic)
had.lm.predic$pred = had_north_lm_pred
