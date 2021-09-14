
library(phyloseq)
library(tidyverse)
library(mgcv)

nthreads <- 24

metadata <- readRDS("../Data/metadata_detrended.RDS")

responses <- c("modified_PC1", "modified_PC2", "modified_PC3", "modified_richness", "modified_shannon", "modified_simpson")

predictors_m1 <- labels(terms(as.formula(modified_PC1~
                                        s(month_int, k=7, bs="cc") + 
                                        s(hydro_year, k=8, m=2) + 
                                        s(rain_monthly) +
                                        s(tempmax_monthly) +
                                        s(rain_annual, k=9) +
                                        s(tempmax_annual, k=9))))

predictors_m2 <- labels(terms(as.formula(modified_PC1~
                                           s(month_int, k=7, bs="cc") + 
                                           s(hydro_year, k=8, m=2) + 
                                           s(rain_monthly) +
                                           s(tempmax_monthly) +
                                           s(rain_annual, k=9) +
                                           s(tempmax_annual, k=9) +
                                           s(time, by=grp, k=20, m=2) + grp +
                                           s(group_size) +
                                           s(area_tot_sqkm) +
                                           s(frac_uniq) +
                                           te(lat, lon, time) +
                                           s(diet_PC1) +
                                           s(diet_PC2) +
                                           s(diet_PC3) +
                                           s(diet_PC4) +
                                           s(diet_PC5) +
                                           s(diet_PC6) +
                                           s(diet_PC7) +
                                           s(diet_PC8) +
                                           s(diet_PC9) +
                                           s(diet_PC10) +
                                           s(diet_PC11) +
                                           s(diet_PC12) +
                                           s(diet_PC13))))

predictors_m3 <- labels(terms(as.formula(modified_PC1~
                                           s(month_int, k=7, bs="cc") + 
                                           s(hydro_year, k=8, m=2) + 
                                           s(rain_monthly) +
                                           s(tempmax_monthly) +
                                           s(rain_annual, k=9) +
                                           s(tempmax_annual, k=9) +
                                           s(time, by=grp, k=20, m=2) + grp +
                                           s(group_size) +
                                           s(area_tot_sqkm) +
                                           s(frac_uniq) +
                                           te(lat, lon, time) +
                                           s(diet_PC1) +
                                           s(diet_PC2) +
                                           s(diet_PC3) +
                                           s(diet_PC4) +
                                           s(diet_PC5) +
                                           s(diet_PC6) +
                                           s(diet_PC7) +
                                           s(diet_PC8) +
                                           s(diet_PC9) +
                                           s(diet_PC10) +
                                           s(diet_PC11) +
                                           s(diet_PC12) +
                                           s(diet_PC13) +
                                           s(time, by=host, k=60, m=1) + host +
                                           s(age) +
                                           s(sex, bs="re") +
                                           s(proportional_rank))))

start.time1 <- Sys.time()

for(y in responses) {
  print(y)
  
  # MODEL 1
  start.time2 <- Sys.time()
  m1 <- bam(reformulate(response=y, predictors_m1), data=metadata, nthreads=c(nthreads,1), discrete=T)
  saveRDS(m1, paste0(str_remove(y,"modified_"),"_m1.RDS"))
  end.time2 <- Sys.time()
  time.taken2 <- end.time2 - start.time2
  print(paste0("finishing model 1 took: ", round(time.taken2, 2), " hours"))
  
  # MODEL 2
  start.time3 <- Sys.time()
  m2 <- bam(reformulate(response=y, predictors_m2), data=metadata, nthreads=c(nthreads,1), discrete=T)
  saveRDS(m2, paste0(str_remove(y,"modified_"),"_m2.RDS"))
  end.time3 <- Sys.time()
  time.taken3 <- end.time3 - start.time3
  print(paste0("finishing model 2 took: ", round(time.taken3, 2), " hours"))
  
  # MODEL 3
  start.time4 <- Sys.time()
  m3 <- bam(reformulate(response=y, predictors_m3), data=metadata, nthreads=c(nthreads,1), discrete=T)
  saveRDS(m3, paste0(str_remove(y,"modified_"),"_m3.RDS"))
  end.time4 <- Sys.time()
  time.taken4 <- end.time4 - start.time4
  print(paste0("finishing model 3 took: ", round(time.taken4, 2), " hours"))
  
}

end.time1 <- Sys.time()
time.taken1 <- end.time1 - start.time1
print(paste0("finishing all models took: ", round(time.taken1, 2), " hours"))


