library(phyloseq)
library(tidyverse)
library(mgcv)

nthreads <- 24

# This is for host_randm1 (we ran these models for host_randm1-10) 

metadata <- readRDS("metadata_shuffled_host_labels.RDS")

responses <- c("modified_PC1", "modified_PC2", "modified_PC3")

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
                                           s(time, by=host_randm1, k=60, m=1) + host_randm1 +
                                           s(age) +
                                           s(sex, bs="re") +
                                           s(proportional_rank))))


for(y in responses) {
  print(y)
  # MODEL 3
  start.time4 <- Sys.time()
  m3 <- bam(reformulate(response=y, predictors_m3), data=metadata, nthreads=c(nthreads,1), discrete=T)
  saveRDS(m3, paste0(str_remove(y,"modified_"),"_m3_host_randm1.RDS"))
  end.time4 <- Sys.time()
  time.taken4 <- end.time4 - start.time4
  print(paste0("finishing model 3 took: ", round(time.taken4, 2), " hours"))
