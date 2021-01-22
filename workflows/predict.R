source("src/data.R")

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
  
  data = data.frame(population=population,
                    area=areasum, 
                    high=habitat.high.areasum, 
                    medium=habitat.med.areasum,
                    low=habitat.low.areasum, 
                    none=habitat.none.areasum,
                    perc.h=habitat.high.perc,
                    perc.m=habitat.med.perc,
                    perc.l=habitat.low.perc,
                    perc.n=habitat.none.perc,
                    weighted.east=weighted.east,
                    weighted.north=weighted.north,
                    twenty.five=quantiles[[1]], 
                    seventy.five=quantiles[[2]]
                    )
  
  return(data)
  
}
species.list = c("Atlantic Cod", "Haddock", "Yellowtail Flounder", "Atlantic Mackerel", "American Butterfish")

years=1978:2016
months=1:12
species=length(species.list)

model=best.models.list[[1]]

habitat.stats = data.frame()
geo.preds = data.frame()

for(i in 2:2){
  s = species.list[i]
  model = best.models.list[[i]]
  
  print(s)
  
  for(y in years){
    print(y)
    for(m in months){
      
      print(m)
      
      edata=load_environmental_data(y, m, s)
      prediction=predict_data(model, edata)
      hs = compute_habitat_stats(prediction)
      
      hs$year = y
      hs$month = m
      hs$species = s
      
      habitat.stats = dplyr::bind_rows(habitat.stats, hs)
      geo.preds = dplyr::bind_rows(geo.preds, prediction)
      
    }
  }
}

#write.csv(geo.preds, file="data/geo_predictions/had_preds.csv")
#write.csv(habitat.stats, file="data/habitat_stats/had_hab.csv")
