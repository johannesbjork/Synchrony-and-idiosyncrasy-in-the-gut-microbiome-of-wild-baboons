
library(phyloseq)
library(tidyverse)
library(mgcv)

nthreads <- 24

metadata <- readRDS("../Data/metadata_detrended.RDS")
families <- read.csv("../Data/family_order.csv")$taxon

predictors <- labels(terms(as.formula(modified_PC1~
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


predictors_updated <- c( predictors[c(-7:-8,-26:-27)][1:6], paste0(predictors[7]," + ", predictors[8]), predictors[c(-7:-8,-26:-27)][7:23], paste0(predictors[26]," + ", predictors[27]), predictors[c(-7:-8,-26:-27)][24:26] )

start.time1 <- Sys.time()
for(y in families) {
  print(y)
  
  full_model <- bam(reformulate(response=y, predictors_updated), data=metadata, control=list(nthreads=nthreads))
  saveRDS(full_model, paste0(str_remove(y,"modified_"),"_full_model.RDS"))

  # Loop through predictors in a leave-one-out fashion
  for(p in predictors_updated) {
    print(p)
    if(p=="s(month_int, k = 7, bs = \"cc\")") {
      idx <- match(names(full_model$sp)[grep(x=names(full_model$sp), pattern="month_int")],names(full_model$sp))
      
      m <- bam(reformulate(response=y, setdiff(predictors_updated, p)), sp=full_model$sp[-idx], data=metadata, nthreads=c(nthreads,1), discrete=TRUE)
      write.csv(data.frame(dev_expl=summary(m, re.test=F)$dev.expl), paste0(str_remove(y,"modified_"),"_month_int.csv"))
      saveRDS(m, paste0(str_remove(y,"modified_"),"_month_int.RDS"))
    } else if(p=="s(hydro_year, k = 8, m = 2)") {
      idx <- match(names(full_model$sp)[grep(x=names(full_model$sp), pattern="hydro_year")],names(full_model$sp))
      m <- bam(reformulate(response=y, setdiff(predictors_updated, p)), sp=full_model$sp[-idx], data=metadata, nthreads=c(nthreads,1), discrete=TRUE)
      write.csv(data.frame(dev_expl=summary(m, re.test=F)$dev.expl), paste0(str_remove(y,"modified_"),"_hydro_year.csv"))
      saveRDS(m, paste0(str_remove(y,"modified_"),"_hydro_year.RDS"))
    } else if(p=="s(rain_annual, k = 9)" ) {
      idx <- match(names(full_model$sp)[grep(x=names(full_model$sp), pattern="rain_annual")],names(full_model$sp))
      m <- bam(reformulate(response=y, setdiff(predictors_updated, p)), sp=full_model$sp[-idx], data=metadata, nthreads=c(nthreads,1), discrete=TRUE)
      write.csv(data.frame(dev_expl=summary(m, re.test=F)$dev.expl), paste0(str_remove(y,"modified_"),"_rain_annual.csv"))
      saveRDS(m, paste0(str_remove(y,"modified_"),"_rain_annual.RDS"))
    } else if(p=="s(tempmax_annual, k = 9)" ) {
      idx <- match(names(full_model$sp)[grep(x=names(full_model$sp), pattern="tempmax_annual")],names(full_model$sp))
      m <- bam(reformulate(response=y, setdiff(predictors_updated, p)), sp=full_model$sp[-idx], data=metadata, nthreads=c(nthreads,1), discrete=TRUE)
      write.csv(data.frame(dev_expl=summary(m, re.test=F)$dev.expl), paste0(str_remove(y,"modified_"),"_tempmax_annual.csv"))
      saveRDS(m, paste0(str_remove(y,"modified_"),"_tempmax_annual.RDS"))
    } else if(p=="s(time, by = grp, k = 20, m = 2) + grp") {
      idx <- match(names(full_model$sp)[grep(x=names(full_model$sp), pattern="grp")],names(full_model$sp))
      m <- bam(reformulate(response=y, setdiff(predictors_updated, p)), sp=full_model$sp[-first(idx):-last(idx)], data=metadata, nthreads=c(nthreads,1), discrete=TRUE)
      write.csv(data.frame(dev_expl=summary(m, re.test=F)$dev.expl), paste0(str_remove(y,"modified_"),"_grp_time.csv"))
      saveRDS(m, paste0(str_remove(y,"modified_"),"_grp_time.RDS"))
    } else if(p=="te(lat, lon, time)") {
      idx <- match(names(full_model$sp)[grep(x=names(full_model$sp), pattern="lat")],names(full_model$sp))
      m <- bam(reformulate(response=y, setdiff(predictors_updated, p)), sp=full_model$sp[-first(idx):-last(idx)], data=metadata, nthreads=c(nthreads,1), discrete=TRUE)
      write.csv(data.frame(dev_expl=summary(m, re.test=F)$dev.expl), paste0(str_remove(y,"modified_"),"_lon_lat_time.csv"))
      saveRDS(m, paste0(str_remove(y,"modified_"),"_lon_lat_time.RDS"))
    } else if(p=="s(time, by = host, k = 60, m = 1) + host") {
      idx <- match(names(full_model$sp)[grep(x=names(full_model$sp), pattern="host")],names(full_model$sp))
      m <- bam(reformulate(response=y, setdiff(predictors_updated, p)), sp=full_model$sp[-first(idx):-last(idx)], data=metadata, nthreads=c(nthreads,1), discrete=TRUE)
      write.csv(data.frame(dev_expl=summary(m, re.test=F)$dev.expl), paste0(str_remove(y,"modified_"),"_host_time.csv"))
      saveRDS(m, paste0(str_remove(y,"modified_"),"_host_time.RDS"))
    } else if(p=="s(sex, bs = \"re\")") {
      idx <- match(names(full_model$sp)[grep(x=names(full_model$sp), pattern="sex")],names(full_model$sp))
      m <- bam(reformulate(response=y, setdiff(predictors_updated, p)), sp=full_model$sp[-idx], data=metadata, nthreads=c(nthreads,1), discrete=TRUE)
      write.csv(data.frame(dev_expl=summary(m, re.test=F)$dev.expl), paste0(str_remove(y,"modified_"),"_sex.csv"))
      saveRDS(m, paste0(str_remove(y,"modified_"),"_sex.RDS"))
      
    } else {
    
      idx <- match(p, names(full_model$sp))
      m <- bam(reformulate(response=y, setdiff(predictors_updated, p)), sp=full_model$sp[-idx], data=metadata, nthreads=c(nthreads,1), discrete=TRUE)
      write.csv(data.frame(dev_expl=summary(m, re.test=F)$dev.expl), paste0(str_remove(y,"modified_"), "_", str_remove_all(string=p, pattern=paste0(c("s\\(","\\)"),collapse="|")),".csv"))
      saveRDS(m, paste0(str_remove(y,"modified_"), "_", str_remove_all(string=p, pattern=paste0(c("s\\(","\\)"),collapse="|")),".RDS"))
    }
    
  }
}
end.time1 <- Sys.time()
time.taken1 <- end.time1 - start.time1
print(paste0("finishing all models took: ", round(time.taken1, 2), " hours"))