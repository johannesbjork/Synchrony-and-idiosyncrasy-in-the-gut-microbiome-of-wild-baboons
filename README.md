
# Index
  * *aitchison_pca_based_analyses.Rmd:* Aitchison PCA analysis
    * *Figs. 2A and 3D; figs S3, S4, S10, and S12*  
  * *taxon_specific_linear_models_and_smooths.Rmd:* Taxon-specific linear models and smooths
    * *Fig. 2C; figs. S7, S8, and S9*
  * *neutral_models.Rmd:* Sloan's Neutral model 
    * *fig. S16*   
  * *autocorrelation_analysis.Rmd:* Autocorrelation analysis
    * *Figs. 3A, B, C; figs. S11*   
  * *taxon_log_ratio_analysis.Rmd:* Host and social group log-ratios
    * *Fig. 3E; figs. S14, S15, S20, and S21*
  * *grooming_analysis.Rmd:* Social connectedness - microbiome similarity analysis
    * *figs. S23 and S24*  
  
 **All of the above analyses apart from *grooming_analysis.Rmd:* use the file ps.RDS in the OSF Data/ repository**
 ***grooming_analysis.Rmd:* uses within_grp_groom.RDS in the OSF Data/ repository** 
  
  * GAMs/
    * hierarchical/
      * *gams_vi_pcs_and_alphadiv.R:* Code for fitting **model P**, **model P+G** and **model P+G+H** to microbiome PC1-3 and 3 alpha diversity metrics
      * *gams_vi_phylum.R:* Code for fitting **model P**, **model P+G** and **model P+G+H** to the clr-transformed relative abundance of 12 phyla  
      * *gams_vi_family.R* Code for fitting **model P**, **model P+G** and **model P+G+H** to the clr-transformed relative abundance of 34 families
    * variable_importance/ 
      * *gams_all_models_PCs_div.R* Code for fitting **model P+G+H** to microbiome PC1-3 and 3 alpha diversity metrics, successively removing one predictor variable at the time, keeping the model otherwise intact
      * *gams_all_models_phylum.R* Code for fitting **model P+G+H** to the clr-transformed relative abundance of 12 phyla, successively removing one predictor variable at the time, keeping the model otherwise intact
      * *gams_all_models_family.R* Code for fitting **model P+G+H** to the clr-transformed relative abundance of 34 families, successively removing one predictor variable at the time, keeping the model otherwise intact
    * santiy_checks/
      * *randomize_host_and_grp_labels.Rmd:* Code to randomize social group membership and host identity
      * *gams_randm_grp_labels.R:* Code to run model P+G+H with randomized social group labels 
      * *gams_randm_grp_labels.R:* Code to run model P+G+H with randomized host labels 
      * *simulation.Rmd:* Code to run a simulation investigating the effect of model complexity on deviance explained
 
**All of the above GAMs use the file metadata_detrended.RDS in the OSF Data/ repository**  

